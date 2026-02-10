# Pre-trained Models

This directory contains pre-trained models for reproducing published results.

## Directory Structure

```
models/
├── E7.5_rep1/
│   ├── combined/
│   │   ├── pca_lsi/
│   │   │   ├── HVG/
│   │   │   │   ├── linear_models/     # OLS, Ridge, Lasso, ElasticNet, RF
│   │   │   │   └── neural_network/    # Deep NN models
│   │   │   └── Random_genes/
│   │   └── wnn/
│   ├── rna_only/
│   └── atac_only/
├── E7.5_rep2/
│   └── ...
└── T_Cells/
    └── ...
```

## Using Pre-trained Models

To use these models for inference on new data:

```r
# Load a pre-trained model
model <- readRDS("models/E7.5_rep1/combined/pca_lsi/HVG/linear_models/elasticnet_models.rds")

# For neural networks (requires TensorFlow/Keras)
library(keras)
nn_model <- load_model_hdf5("models/E7.5_rep1/combined/pca_lsi/HVG/neural_network/model.h5")
```

## Model Formats

| Model Type | File Format | R Package |
|------------|-------------|-----------|
| Linear models | `.rds` | glmnet |
| Random Forest | `.rds` | ranger |
| Neural Network | `.h5` / `.keras` | keras/tensorflow |

## Notes

- Models were trained using the configuration files in `data/target_genes/`
- Ensure you use the same target gene lists when evaluating
- Neural network models require compatible TensorFlow version (see requirements.txt)
