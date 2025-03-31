# Coupled Tensor-Tensor Completion Method (CTTC)

## Author
- **Maryam Bagherian**

## Date
- **March 31, 2025**

## Description
This repository implements the **Coupled Tensor-Tensor Completion (CTTC)** method for predicting drug-target interactions. The code uses coupled matrix-tensor and tensor-tensor completion methods to predict missing entries in tensor data. This approach is part of the work to improve drug-target interaction prediction using multi-dimensional tensor data.

## References

1. Bagherian, M., et al., "Coupled tensor-tensor completion method", under revision. 
2. Bagherian, M., et al., "Coupled matrix–matrix and coupled tensor–matrix completion methods for predicting drug–target interactions."  
   Briefings in Bioinformatics 22, no. 2 (2021): 2161-2171, doi: [10.1093/bib/bbaa025](https://doi.org/10.1093/bib/bbaa025).
3. Bagherian, M., "Tensor denoising via dual Schatten norms." Optimization Letters 18, no. 5 (2024): 1285-1301, doi: [10.1007/s11590-023-02068-8](https://doi.org/10.1007/s11590-023-02068-8)
   
Some MATLAB commandas have been used from the following refrence: 

4. Chen, Yi-Lei, Chiou-Ting Hsu, and Hong-Yuan Mark Liao. "Simultaneous tensor decomposition and completion using factor priors."  
   IEEE Transactions on Pattern Analysis and Machine Intelligence 36, no. 3 (2013): 577-591.

## Code Description

This code implements the Coupled Tensor-Tensor Completion (CTTC) method for tensor completion. It provides flexible configuration options, including parameters for scaling, regularization, and iteration settings. The code is designed to handle missing data by randomly selecting a missing rate for tensor completion, enabling the prediction of missing values based on available data.

The CTTC algorithm leverages the Coupled Matrix-Matrix Completion (CMMC) function to learn local distance functions for each tensor mode, utilizing the CoupledTensorCompletion function. This process allows the model to adaptively learn the appropriate distance metrics for each mode of the tensor. Additionally, a Mahalanobis distance function is learned through the core tensor, Z, enabling more accurate predictions of missing entries by capturing the underlying relationships in the data.

The algorithm evaluates its performance by comparing the predicted values against the ground truth, allowing for optimization and fine-tuning of the method. The method aims to minimize the reconstruction error, improving the accuracy of the tensor completion and enhancing the model's predictive capabilities.

### Functionality
- **CTTC.m**: Main function to perform the Coupled Tensor-Tensor Completion method. It requires the input tensor, side information, and several other parameters.
- **Supporting Functions**: Additional functions are used for initialization, metric function construction, optimization, and more.

### Key Parameters
- `X`: Initial tensor data (drug-target interaction data).
- `Side`: Side information, including drug similarity, gene similarity, and cell similarity matrices.
- `scale`: Vector of scaling factors.
- `epsilon`: Vector of regularization parameters.
- `iter`: Number of iterations for the optimization process.
- `m_rate`: Missing rate (proportion of missing entries in the tensor).
  
### Running the Code
1. Load the required datasets into the MATLAB environment.
2. Set the parameters for the tensor completion process, such as `scale`, `epsilon`, `iter`, and `m_rate`.
3. Call the `CTTC.m` function with appropriate arguments.
4. The code will return the completed tensor and evaluate the performance based on the residual error.

## Requirements

- MATLAB (version R2019b or later)
- Required MATLAB toolboxes (if any) for optimization and matrix operations

## Data Subset
A small subset of the data used in this project is provided in the data folder for ease of access and testing purposes. This subset contains representative samples of the full dataset and can be used to quickly test the code without the need to download the entire dataset.

To access this subset, navigate to the data folder in the repository, where you'll find the following files:
- cell_sim.mat
- drug_sim_mat_dense_subset
- filtered_target_tensor_1_subset
- gene_sim_mat_normalized_subset

Please note that these files represent a small portion of the full dataset and may not reflect all the variations present in the complete dataset. For access to the full dataset, please refer to the instructions provided in the manuscript [1]. 

## For questions or issues related to the code, please contact Maryam Bagherian at [maryambagherian@isu.deu].

## Example Usage

```matlab
% Load sample data
X1 = load("path/to/filtered_target_tensor_1.mat");
X2 = load("path/to/filtered_target_tensor_2.mat");
X3 = load("path/to/filtered_target_tensor_3.mat");

% Combine tensors
X = cat(2, X1, X2, X3);

% Set parameters
Side = {gene_sim_matrix, drug_sim_matrix, cell_sim_matrix};  % Example side info
scale = [1, 1, 1];  % Scaling factors
epsilon = [0.1, 0.1, 0.1];  % Regularization parameters
iter = 50;  % Number of iterations
m_rate = 0.2;  % Missing data rate

% Call the CTTC function
completed_tensor = CTTC(X, Side, iter, m_rate);
``` 
## How to Cite
If you use this code in your research or projects, please cite the following paper:
```bibtex
@article{BagherianCTTC,
  author = {Maryam Bagherian, Albert Hung, Ivo Dinov, Joshua Welch},
  title = {Coupled tensor-tensor completion method},
  journal = {to be added later},
  volume = {},
  number = {},
  pages = {},
  year = {2025},
  doi = {}
}

For additional works or specific references related to this code, you may also refer to the refrences above. 
```





