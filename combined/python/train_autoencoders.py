#!/usr/bin/env python3
# ============================================================================
# Step_025: Dimensionality Reduction using Autoencoders (Automated)
# ============================================================================
# This script performs dimensionality reduction using:
#   1. scVI for RNA data
#   2. PeakVI for ATAC data
#   3. MultiVI for joint RNA+ATAC data
#
# Training is done on TRAINING data only, then projected to all cells.
# This prevents data leakage from validation/test sets.
#
# Input: MuData file (.h5mu) with RNA and ATAC modalities and data_split column
# Output: Latent representations for train/validation/test splits
#
# Usage:
#   python Step_025.Dimensionality_Reduction_Autoencoders_Automated.py
#
# Note: Requires config.R to be parsed for paths. This script reads config
#       values from environment variables or uses defaults.
# ============================================================================

import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import muon as mu
import os
import gc
import sys
import argparse
import torch

# ============================================================================
# CONFIGURATION
# ============================================================================
# These can be overridden via command line arguments or environment variables

def check_gpu_available():
    """Check if GPU is available and working."""
    if not torch.cuda.is_available():
        return False, "CUDA not available"
    
    try:
        # Try to create a small tensor on GPU to verify it works
        device = torch.device('cuda')
        test_tensor = torch.zeros(1, device=device)
        del test_tensor
        torch.cuda.empty_cache()
        return True, f"GPU available: {torch.cuda.get_device_name(0)}"
    except Exception as e:
        return False, f"GPU test failed: {str(e)}"

def get_config():
    """Parse command line arguments or use defaults."""
    parser = argparse.ArgumentParser(description='Autoencoder dimensionality reduction')
    
    parser.add_argument('--input', type=str, 
                        help='Path to input .h5mu file')
    parser.add_argument('--output-dir', type=str,
                        help='Output directory for latent representations')
    parser.add_argument('--sample-name', type=str, default='sample',
                        help='Sample name for file naming')
    parser.add_argument('--n-latent', type=int, default=30,
                        help='Number of latent dimensions')
    parser.add_argument('--max-epochs', type=int, default=500,
                        help='Maximum training epochs')
    parser.add_argument('--n-hvg', type=int, default=3000,
                        help='Number of highly variable genes for RNA')
    parser.add_argument('--atac-percentile', type=float, default=5,
                        help='Percentile cutoff for ATAC peak selection (q05 = 5)')
    parser.add_argument('--use-cpu', action='store_true',
                        help='Force CPU training (slower but more compatible)')
    
    args = parser.parse_args()
    
    # If arguments not provided, try to read from environment or use defaults
    config = {
        'input_file': args.input or os.environ.get('INPUT_MUDATA', None),
        'output_dir': args.output_dir or os.environ.get('OUTPUT_DIR', './latent_output'),
        'sample_name': args.sample_name or os.environ.get('SAMPLE_NAME', 'sample'),
        'n_latent_dims': args.n_latent,
        'max_epochs': args.max_epochs,
        'n_hvg': args.n_hvg,
        'atac_percentile': args.atac_percentile,
        'use_cpu': args.use_cpu,
    }
    
    return config


def prepare_adata_for_projection(adata):
    """
    Manually injects internal scVI columns for projection.
    This tricks the model into accepting a new AnnData object without re-training.
    CRITICAL for MultiVI projection to work correctly!
    """
    # 1. _indices: Unique integer ID for each cell
    adata.obs['_indices'] = np.arange(len(adata)).astype(int)
    
    # 2. _scvi_batch: All 0s (assuming single batch training)
    adata.obs['_scvi_batch'] = np.zeros(len(adata), dtype=int)
    
    # 3. _scvi_labels: All 0s (assuming unsupervised training)
    adata.obs['_scvi_labels'] = np.zeros(len(adata), dtype=int)
    
    return adata


def main():
    config = get_config()
    
    # Validate input
    if config['input_file'] is None:
        print("ERROR: No input file specified!")
        print("Usage: python Step_045.Dimensionality_Reduction_Autoencoders_Automated.py --input <file.h5mu> --output-dir <dir>")
        sys.exit(1)
    
    input_file = config['input_file']
    output_dir = config['output_dir']
    sample_name = config['sample_name']
    n_latent_dims = config['n_latent_dims']
    max_epochs = config['max_epochs']
    n_hvg = config['n_hvg']
    atac_percentile = config['atac_percentile']
    use_cpu = config['use_cpu']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Check GPU availability
    gpu_available, gpu_message = check_gpu_available()
    if use_cpu:
        use_gpu = False
        print("\n⚠ Forcing CPU mode (--use-cpu flag set)")
    elif gpu_available:
        use_gpu = True
        print(f"\n✓ {gpu_message}")
    else:
        use_gpu = False
        print(f"\n⚠ {gpu_message}")
        print("  Falling back to CPU training (slower but compatible)")

    def try_train_model(train_func, *args, **kwargs):
        """Try to train or project with GPU, fallback to CPU on CUDA error."""
        try:
            return train_func(*args, **kwargs, accelerator='gpu' if use_gpu else 'cpu')
        except RuntimeError as e:
            if 'CUDA error' in str(e) or 'cuda' in str(e).lower():
                print("\n[AutoFallback] CUDA error detected. Retrying on CPU...\n")
                return train_func(*args, **kwargs, accelerator='cpu')
            else:
                raise
    
    print("=" * 80)
    print("Dimensionality Reduction using Autoencoders")
    print("=" * 80)
    print(f"\\nConfiguration:")
    print(f"  Input file: {input_file}")
    print(f"  Output directory: {output_dir}")
    print(f"  Sample name: {sample_name}")
    print(f"  Latent dimensions: {n_latent_dims}")
    print(f"  Max epochs: {max_epochs}")
    print(f"  HVG count: {n_hvg}")
    print(f"  ATAC percentile cutoff: q{atac_percentile}")
    print(f"  Device: {'GPU' if use_gpu else 'CPU'}")
    
    # =========================================================================
    # 0. LOAD DATA
    # =========================================================================
    print("\n" + "=" * 80)
    print("Loading MuData...")
    try:
        mdata = mu.read_h5mu(input_file)
        print(f"✓ Loaded MuData with {mdata.n_obs} cells")
    except Exception as e:
        print(f"✗ Error loading data: {e}")
        sys.exit(1)
    
    # Extract RNA and ATAC modalities
    rna_data = mdata.mod['rna'].copy()
    atac_data = mdata.mod['atac'].copy()
    
    # Verify data_split column exists
    if "data_split" not in rna_data.obs.columns:
        print("✗ Error: data_split column not found in RNA modality!")
        sys.exit(1)
    
    print("\nData split distribution:")
    print(rna_data.obs["data_split"].value_counts())
    
    # Create split masks
    is_train = rna_data.obs["data_split"] == "train"
    is_val = rna_data.obs["data_split"] == "validation"
    is_test = rna_data.obs["data_split"] == "test"
    
    print(f"\nSplit sizes: Train={is_train.sum()}, Val={is_val.sum()}, Test={is_test.sum()}")
    
    # =========================================================================
    # 1. FEATURE SELECTION (FROM TRAINING DATA ONLY)
    # =========================================================================
    print("\n" + "=" * 80)
    print("Performing Feature Selection on TRAINING data only...")
    print("=" * 80)
    
    # RNA: Highly Variable Genes
    print("\nSelecting Highly Variable Genes for RNA...")
    try:
        rna_train_temp = rna_data[is_train].copy()
        
        sc.pp.highly_variable_genes(
            rna_train_temp,
            flavor="seurat_v3",
            n_top_genes=n_hvg,
            subset=False
        )
        
        rna_genes_to_keep = rna_train_temp.var_names[rna_train_temp.var['highly_variable']].tolist()
        print(f"✓ Selected {len(rna_genes_to_keep)} highly variable genes for RNA")
        
        del rna_train_temp
        gc.collect()
    except Exception as e:
        print(f"✗ Error in RNA feature selection: {e}")
        sys.exit(1)
    
    # ATAC: Percentile-based selection
    print(f"\nSelecting informative peaks for ATAC (q{atac_percentile} method)...")
    try:
        atac_train_temp = atac_data[is_train].copy()
        
        if hasattr(atac_train_temp.X, "todense"):
            atac_dense = atac_train_temp.X.todense()
        else:
            atac_dense = atac_train_temp.X
        
        sum_across_cells = np.array(atac_dense.sum(axis=0)).flatten()
        threshold = np.percentile(sum_across_cells, atac_percentile)
        peaks_to_keep = sum_across_cells >= threshold
        atac_peaks_to_keep = atac_train_temp.var_names[peaks_to_keep].tolist()
        
        print(f"✓ Q{atac_percentile} threshold: {threshold:.2f}")
        print(f"✓ Selected {len(atac_peaks_to_keep)} informative peaks for ATAC")
        
        del atac_train_temp, atac_dense, sum_across_cells, peaks_to_keep
        gc.collect()
    except Exception as e:
        print(f"✗ Error in ATAC feature selection: {e}")
        sys.exit(1)
    
    # Save feature lists for reproducibility
    with open(os.path.join(output_dir, "selected_rna_genes.txt"), 'w') as f:
        f.write('\n'.join(rna_genes_to_keep))
    with open(os.path.join(output_dir, "selected_atac_peaks.txt"), 'w') as f:
        f.write('\n'.join(atac_peaks_to_keep))
    print("\n✓ Saved feature lists to output directory")
    
    # Apply feature selection to full datasets
    print("\nApplying feature selection to full datasets...")
    rna_data = rna_data[:, rna_genes_to_keep].copy()
    atac_data = atac_data[:, atac_peaks_to_keep].copy()
    print(f"✓ RNA data: {rna_data.shape}")
    print(f"✓ ATAC data: {atac_data.shape}")
    
    # =========================================================================
    # 2. scVI FOR RNA
    # =========================================================================
    print("\n" + "=" * 80)
    print("=== scVI: RNA Dimensionality Reduction ===")
    print("=" * 80)
    
    try:
        rna_train = rna_data[is_train].copy()
        
        print(f"\nSetting up scVI on {rna_train.n_obs} training cells...")
        scvi.model.SCVI.setup_anndata(rna_train)
        
        print("Training scVI model...")
        scvi_model = scvi.model.SCVI(rna_train, n_latent=n_latent_dims)
        try_train_model(scvi_model.train, max_epochs=max_epochs)
        print("✓ Training complete")
        
        # Save model
        scvi_model_path = os.path.join(output_dir, "scvi_model")
        scvi_model.save(scvi_model_path, overwrite=True)
        print(f"✓ Saved scVI model to {scvi_model_path}")
        
        # Project to ALL cells
        print("Projecting ALL cells to scVI latent space...")
        latent_rna_all = scvi_model.get_latent_representation(rna_data)
        latent_rna_all_df = pd.DataFrame(
            latent_rna_all,
            index=rna_data.obs_names,
            columns=[f'scVI_RNA_{i}' for i in range(n_latent_dims)]
        )
        latent_rna_all_df['data_split'] = rna_data.obs['data_split'].values
        
        # Save
        latent_rna_all_df.to_csv(os.path.join(output_dir, "latent_scvi_rna_all.csv"))
        
        latent_rna_train = latent_rna_all_df[latent_rna_all_df['data_split'] == 'train']
        latent_rna_val = latent_rna_all_df[latent_rna_all_df['data_split'] == 'validation']
        latent_rna_test = latent_rna_all_df[latent_rna_all_df['data_split'] == 'test']
        
        latent_rna_train.to_csv(os.path.join(output_dir, "latent_scvi_rna_train.csv"))
        latent_rna_val.to_csv(os.path.join(output_dir, "latent_scvi_rna_val.csv"))
        latent_rna_test.to_csv(os.path.join(output_dir, "latent_scvi_rna_test.csv"))
        
        print(f"✓ scVI RNA complete: Train={len(latent_rna_train)}, Val={len(latent_rna_val)}, Test={len(latent_rna_test)}")
        
        del rna_train, scvi_model
        gc.collect()
        
    except Exception as e:
        print(f"✗ Error in scVI processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # =========================================================================
    # 3. PeakVI FOR ATAC
    # =========================================================================
    print("\n" + "=" * 80)
    print("=== PeakVI: ATAC Dimensionality Reduction ===")
    print("=" * 80)
    
    try:
        atac_train = atac_data[is_train].copy()
        
        print(f"\nSetting up PeakVI on {atac_train.n_obs} training cells...")
        scvi.model.PEAKVI.setup_anndata(atac_train)
        
        print("Training PeakVI model...")
        peakvi_model = scvi.model.PEAKVI(atac_train, n_latent=n_latent_dims)
        try_train_model(peakvi_model.train, max_epochs=max_epochs)
        print("✓ Training complete")
        
        # Save model
        peakvi_model_path = os.path.join(output_dir, "peakvi_model")
        peakvi_model.save(peakvi_model_path, overwrite=True)
        print(f"✓ Saved PeakVI model to {peakvi_model_path}")
        
        # Project to ALL cells
        print("Projecting ALL cells to PeakVI latent space...")
        latent_atac_all = peakvi_model.get_latent_representation(atac_data)
        latent_atac_all_df = pd.DataFrame(
            latent_atac_all,
            index=atac_data.obs_names,
            columns=[f'PeakVI_ATAC_{i}' for i in range(n_latent_dims)]
        )
        latent_atac_all_df['data_split'] = atac_data.obs['data_split'].values
        
        # Save
        latent_atac_all_df.to_csv(os.path.join(output_dir, "latent_peakvi_atac_all.csv"))
        
        latent_atac_train = latent_atac_all_df[latent_atac_all_df['data_split'] == 'train']
        latent_atac_val = latent_atac_all_df[latent_atac_all_df['data_split'] == 'validation']
        latent_atac_test = latent_atac_all_df[latent_atac_all_df['data_split'] == 'test']
        
        latent_atac_train.to_csv(os.path.join(output_dir, "latent_peakvi_atac_train.csv"))
        latent_atac_val.to_csv(os.path.join(output_dir, "latent_peakvi_atac_val.csv"))
        latent_atac_test.to_csv(os.path.join(output_dir, "latent_peakvi_atac_test.csv"))
        
        print(f"✓ PeakVI ATAC complete: Train={len(latent_atac_train)}, Val={len(latent_atac_val)}, Test={len(latent_atac_test)}")
        
        del atac_train, peakvi_model
        gc.collect()
        
    except Exception as e:
        print(f"✗ Error in PeakVI processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Free RNA and ATAC data before MultiVI
    del rna_data, atac_data
    gc.collect()
    
    # =========================================================================
    # 4. MultiVI (JOINT RNA+ATAC) 
    # =========================================================================
    print("\n" + "=" * 80)
    print("=== MultiVI: Joint RNA+ATAC Dimensionality Reduction ===")
    print("=" * 80)
    
    try:
        # Prepare MuData with selected features
        print("\nPreparing MuData for MultiVI...")
        mdata.mod['rna'] = mdata.mod['rna'][:, rna_genes_to_keep].copy()
        mdata.mod['atac'] = mdata.mod['atac'][:, atac_peaks_to_keep].copy()
        mdata.update()
        
        print(f"✓ RNA features: {mdata.mod['rna'].n_vars}")
        print(f"✓ ATAC features: {mdata.mod['atac'].n_vars}")
        
        # Ensure data_split column exists in global obs
        if "data_split" not in mdata.obs.columns:
            print("Adding data_split from RNA modality to global obs...")
            mdata.obs["data_split"] = mdata.mod["rna"].obs["data_split"].values
        
        # Create training subset
        is_train_multi = mdata.obs["data_split"] == "train"
        print(f"\nCreating training subset ({is_train_multi.sum()} cells)...")
        mdata_train = mdata[is_train_multi].copy()
        
        # Setup MultiVI on training data
        print("Setting up MultiVI on training data...")
        scvi.model.MULTIVI.setup_mudata(
            mdata_train,
            modalities={
                "rna_layer": "rna",
                "atac_layer": "atac",
            },
        )
        
        # Train model
        print(f"Training MultiVI on {mdata_train.n_obs} training cells...")
        multivi_model = scvi.model.MULTIVI(mdata_train, n_latent=n_latent_dims)
        try_train_model(multivi_model.train, max_epochs=max_epochs)
        print("✓ Training complete")
        
        # Save trained model
        multivi_model_path = os.path.join(output_dir, "multivi_model")
        print(f"\nSaving trained MultiVI model to {multivi_model_path}...")
        multivi_model.save(multivi_model_path, overwrite=True)
        print("✓ Model saved")
        
        del mdata_train
        gc.collect()
        
        # =====================================================================
        # Projection approach for MultiVI
        # =====================================================================
        print("\n--- Preparing data for MultiVI projection ---")
        print("Using organize_multiome_anndatas() + metadata patching...")
        
        # Step 1: Organize the full multiome data
        combined_data_full = scvi.data.organize_multiome_anndatas(mdata)
        combined_data_full.obs['data_split'] = mdata.obs['data_split'].values
        print(f"✓ Organized data: {combined_data_full.shape}")
        
        # Step 2: Patch metadata for projection (CRITICAL!)
        print("Patching metadata for projection...")
        combined_data_full = prepare_adata_for_projection(combined_data_full)
        print("✓ Metadata patched with _indices, _scvi_batch, _scvi_labels")
        
        # Step 3: Load model with patched adata and project
        print("\nLoading model with patched data and projecting ALL cells...")
        multivi_model_loaded = scvi.model.MULTIVI.load(multivi_model_path, adata=combined_data_full)

        def try_get_latent(model, adata):
            try:
                return model.get_latent_representation(adata)
            except RuntimeError as e:
                if 'CUDA error' in str(e) or 'cuda' in str(e).lower():
                    print("\n[AutoFallback] CUDA error detected during projection. Retrying on CPU...\n")
                    model.to_device('cpu')
                    return model.get_latent_representation(adata)
                else:
                    raise

        latent_multi_all = try_get_latent(multivi_model_loaded, combined_data_full)
        
        latent_multi_all_df = pd.DataFrame(
            latent_multi_all,
            index=combined_data_full.obs_names,
            columns=[f'MultiVI_{i}' for i in range(n_latent_dims)]
        )
        latent_multi_all_df['data_split'] = combined_data_full.obs['data_split'].values
        
        print("✓ Projection complete")
        
        # Save
        print("Saving MultiVI latent representations...")
        latent_multi_all_df.to_csv(os.path.join(output_dir, "latent_multivi_all.csv"))
        
        latent_multi_train = latent_multi_all_df[latent_multi_all_df['data_split'] == 'train']
        latent_multi_val = latent_multi_all_df[latent_multi_all_df['data_split'] == 'validation']
        latent_multi_test = latent_multi_all_df[latent_multi_all_df['data_split'] == 'test']
        
        latent_multi_train.to_csv(os.path.join(output_dir, "latent_multivi_train.csv"))
        latent_multi_val.to_csv(os.path.join(output_dir, "latent_multivi_val.csv"))
        latent_multi_test.to_csv(os.path.join(output_dir, "latent_multivi_test.csv"))
        
        print(f"✓ MultiVI complete: Train={len(latent_multi_train)}, Val={len(latent_multi_val)}, Test={len(latent_multi_test)}")
        
    except Exception as e:
        print(f"✗ Error in MultiVI processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # =========================================================================
    # 5. UMAP VISUALIZATION OF LATENT SPACES
    # =========================================================================
    print("\n" + "=" * 80)
    print("=== UMAP Visualization of Latent Spaces ===")
    print("=" * 80)
    
    try:
        import matplotlib.pyplot as plt
        from umap import UMAP
        
        # Create plots directory
        plots_dir = os.path.join(output_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        # Color palette for splits
        split_colors = {'train': '#1f77b4', 'validation': '#ff7f0e', 'test': '#2ca02c'}
        
        def plot_umap(latent_df, title, filename, n_latent):
            """Generate UMAP plot for a latent space."""
            print(f"\nGenerating UMAP for {title}...")
            
            # Extract latent dimensions (exclude data_split column)
            latent_cols = [col for col in latent_df.columns if col != 'data_split']
            X = latent_df[latent_cols].values
            
            # Compute UMAP
            umap_model = UMAP(n_components=2, random_state=42, n_neighbors=30, min_dist=0.3)
            umap_coords = umap_model.fit_transform(X)
            
            # Create plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Plot each split separately
            for split in ['train', 'validation', 'test']:
                mask = latent_df['data_split'] == split
                if mask.sum() > 0:
                    ax.scatter(
                        umap_coords[mask, 0],
                        umap_coords[mask, 1],
                        c=split_colors[split],
                        label=f'{split} (n={mask.sum()})',
                        alpha=0.5,
                        s=5
                    )
            
            ax.set_xlabel('UMAP 1', fontsize=12)
            ax.set_ylabel('UMAP 2', fontsize=12)
            ax.set_title(f'{title}\n({n_latent} latent dimensions)', fontsize=14)
            ax.legend(loc='upper right', markerscale=3)
            
            # Remove spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, filename), dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"✓ Saved: {filename}")
            
            return umap_coords
        
        # Plot scVI (RNA)
        umap_scvi = plot_umap(
            latent_rna_all_df, 
            "scVI Latent Space (RNA only)", 
            "umap_scvi_rna.png",
            n_latent_dims
        )
        
        # Plot PeakVI (ATAC)
        umap_peakvi = plot_umap(
            latent_atac_all_df, 
            "PeakVI Latent Space (ATAC only)", 
            "umap_peakvi_atac.png",
            n_latent_dims
        )
        
        # Plot MultiVI (joint)
        umap_multivi = plot_umap(
            latent_multi_all_df, 
            "MultiVI Latent Space (RNA + ATAC joint)", 
            "umap_multivi_joint.png",
            n_latent_dims
        )
        
        # Create combined comparison figure
        print("\nGenerating combined comparison plot...")
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        umaps = [
            (umap_scvi, latent_rna_all_df, "scVI (RNA)"),
            (umap_peakvi, latent_atac_all_df, "PeakVI (ATAC)"),
            (umap_multivi, latent_multi_all_df, "MultiVI (Joint)")
        ]
        
        for ax, (umap_coords, latent_df, title) in zip(axes, umaps):
            for split in ['train', 'validation', 'test']:
                mask = latent_df['data_split'] == split
                if mask.sum() > 0:
                    ax.scatter(
                        umap_coords[mask, 0],
                        umap_coords[mask, 1],
                        c=split_colors[split],
                        label=f'{split} (n={mask.sum()})',
                        alpha=0.5,
                        s=3
                    )
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.set_title(title, fontsize=12)
            ax.legend(loc='upper right', markerscale=3, fontsize=8)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        plt.suptitle(f'Comparison of Latent Space Embeddings ({sample_name})', fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "umap_comparison_all.png"), dpi=150, bbox_inches='tight')
        plt.close()
        
        print("✓ Saved: umap_comparison_all.png")
        
        # Save UMAP coordinates
        umap_scvi_df = pd.DataFrame(umap_scvi, index=latent_rna_all_df.index, columns=['UMAP1', 'UMAP2'])
        umap_scvi_df['data_split'] = latent_rna_all_df['data_split'].values
        umap_scvi_df.to_csv(os.path.join(plots_dir, "umap_coords_scvi.csv"))
        
        umap_peakvi_df = pd.DataFrame(umap_peakvi, index=latent_atac_all_df.index, columns=['UMAP1', 'UMAP2'])
        umap_peakvi_df['data_split'] = latent_atac_all_df['data_split'].values
        umap_peakvi_df.to_csv(os.path.join(plots_dir, "umap_coords_peakvi.csv"))
        
        umap_multivi_df = pd.DataFrame(umap_multivi, index=latent_multi_all_df.index, columns=['UMAP1', 'UMAP2'])
        umap_multivi_df['data_split'] = latent_multi_all_df['data_split'].values
        umap_multivi_df.to_csv(os.path.join(plots_dir, "umap_coords_multivi.csv"))
        
        print("\n✓ UMAP coordinates saved to plots directory")
        
    except ImportError as e:
        print(f"\nWarning: Could not generate UMAP plots. Missing package: {e}")
        print("Install with: pip install umap-learn matplotlib")
    except Exception as e:
        print(f"\nWarning: Error generating UMAP plots: {e}")
        import traceback
        traceback.print_exc()
        print("Continuing without UMAP visualization...")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 80)
    print("=== DIMENSIONALITY REDUCTION COMPLETE ===")
    print("=" * 80)
    print(f"\nAll latent representations saved to: {output_dir}")
    print("\nFiles created:")
    print("  scVI (RNA):")
    print(f"    - latent_scvi_rna_all.csv ({len(latent_rna_all_df)} cells)")
    print(f"    - latent_scvi_rna_train.csv")
    print(f"    - latent_scvi_rna_val.csv")
    print(f"    - latent_scvi_rna_test.csv")
    print("\n  PeakVI (ATAC):")
    print(f"    - latent_peakvi_atac_all.csv ({len(latent_atac_all_df)} cells)")
    print(f"    - latent_peakvi_atac_train.csv")
    print(f"    - latent_peakvi_atac_val.csv")
    print(f"    - latent_peakvi_atac_test.csv")
    print("\n  MultiVI (joint):")
    print(f"    - latent_multivi_all.csv ({len(latent_multi_all_df)} cells)")
    print(f"    - latent_multivi_train.csv")
    print(f"    - latent_multivi_val.csv")
    print(f"    - latent_multivi_test.csv")
    print("\n  Saved models:")
    print(f"    - {os.path.join(output_dir, 'scvi_model')}/")
    print(f"    - {os.path.join(output_dir, 'peakvi_model')}/")
    print(f"    - {os.path.join(output_dir, 'multivi_model')}/")
    print("\n  Feature lists:")
    print(f"    - selected_rna_genes.txt ({len(rna_genes_to_keep)} genes)")
    print(f"    - selected_atac_peaks.txt ({len(atac_peaks_to_keep)} peaks)")
    print("\n  UMAP plots:")
    print(f"    - plots/umap_scvi_rna.png")
    print(f"    - plots/umap_peakvi_atac.png")
    print(f"    - plots/umap_multivi_joint.png")
    print(f"    - plots/umap_comparison_all.png")
    print(f"    - plots/umap_coords_*.csv (coordinates)")
    print("\n" + "=" * 80)
    print("✓ All embeddings ready for downstream analysis!")
    print("=" * 80)


if __name__ == "__main__":
    main()
