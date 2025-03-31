# Coupled Tensor-Tensor Completion (CTTC) for Drug-Target Interaction Prediction

## Author
- **Maryam Bagherian**

## Date
- **March 31, 2025**

## Description
This repository implements the **Coupled Tensor-Tensor Completion (CTTC)** method for predicting drug-target interactions. The code uses coupled matrix-tensor and tensor-tensor completion methods to predict missing entries in tensor data. This approach is part of the work to improve drug-target interaction prediction using multi-dimensional tensor data.

## References

1. Bagherian, M., et al., "Coupled tensor-tensor completion methods", under revision. 
2. Bagherian, M., et al., "Coupled matrix–matrix and coupled tensor–matrix completion methods for predicting drug–target interactions."  
   Briefings in Bioinformatics 22, no. 2 (2021): 2161-2171, doi: [10.1093/bib/bbaa025](https://doi.org/10.1093/bib/bbaa025).
3. Bagherian, M., "Tensor denoising via dual Schatten norms." Optimization Letters 18, no. 5 (2024): 1285-1301.
4. Chen, Yi-Lei, Chiou-Ting Hsu, and Hong-Yuan Mark Liao. "Simultaneous tensor decomposition and completion using factor priors."  
   IEEE Transactions on Pattern Analysis and Machine Intelligence 36, no. 3 (2013): 577-591.

## Code Description

This code implements the CTTC method for tensor completion. It includes various parameters for scaling, regularization, and iteration settings. The code also handles missing data by randomly selecting a missing rate for tensor completion and applies the CTTC algorithm to predict missing values.

### Functionality
- **CTTC.m**: Main function to perform the Coupled Tensor-Tensor Completion method. It requires the input tensor, side information, and several other parameters.
- **Supporting Functions**: Additional functions are used for initialization, graph construction, optimization, and more.

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



For questions or issues related to the code, please contact Maryam Bagherian at [maryambagherian@isu.deu].
