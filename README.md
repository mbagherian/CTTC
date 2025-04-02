# Coupled Tensor-Tensor Completion Method (CTTC)

## Author
- **Maryam Bagherian**

## Date
- **March 31, 2025**

## Description
This repository implements the **Coupled Tensor-Tensor Completion (CTTC)** method for predicting drug-target interactions. The code uses coupled matrix-tensor and tensor-tensor completion methods to predict missing entries in tensor data. This approach is part of the work to improve drug-target interaction prediction using multi-dimensional tensor data.

## References

1. Bagherian, M., et al., "Coupled tensor-tensor completion method", under revision. 
2. Bagherian, M., et al., "Coupled matrix–matrix and coupled tensor–matrix completion methods for predicting drug–target interactions." Briefings in Bioinformatics 22, no. 2 (2021): 2161-2171, doi: [10.1093/bib/bbaa025](https://doi.org/10.1093/bib/bbaa025).
3. Bagherian, M., "Tensor denoising via dual Schatten norms." Optimization Letters 18, no. 5 (2024): 1285-1301, doi: [10.1007/s11590-023-02068-8](https://doi.org/10.1007/s11590-023-02068-8)
   
Some MATLAB updates have been used from the following refrence: 

4. Chen, Yi-Lei, Chiou-Ting Hsu, and Hong-Yuan Mark Liao. "Simultaneous tensor decomposition and completion using factor priors."  
   IEEE Transactions on Pattern Analysis and Machine Intelligence 36, no. 3 (2013): 577-591.

## Code Description

This code implements the Coupled Tensor-Tensor Completion (CTTC) method for tensor completion. It provides flexible configuration options, including parameters for scaling, regularization, and iteration settings. The code is designed to handle missing data by randomly selecting a missing rate for tensor completion, enabling the prediction of missing values based on available data.

The CTTC algorithm generalizes the Coupled Matrix-Matrix Completion (CMMC) and Coupled Tensor-Matrix Completion (CTMC)functions to learn local distance functions for each tensor mode, utilizing the CoupledTensorCompletion function. This process allows the model to adaptively learn the appropriate distance metrics for each mode of the tensor. Additionally, a Mahalanobis distance function is learned through the core tensor, Z, enabling more accurate predictions of missing entries by capturing the underlying relationships in the data.

The algorithm evaluates its performance by comparing the predicted values against the ground truth, allowing for optimization and fine-tuning of the method. The method aims to minimize the reconstruction error, improving the accuracy of the tensor completion and enhancing the model's predictive capabilities.

### Functionality
- **CTTC.m**: Main function to perform the Coupled Tensor-Tensor Completion method. It requires the input tensor, side information, and several other parameters.
- **Supporting Functions**: Additional functions are used for initialization, Distance Metric Function construction, optimization, and more.

### Key Parameters
- `X`: An m-mode initial tensor data (e.g., cell-drug-target), could be 3D or higher.
- `Side`: A cell containing m side information tensors, including drug similarity, gene similarity, and cell similarity tensors, specific to this application. 
- `scale`: Vector of length m containing scaling factors fro each mode. 
- `epsilon`: Vector of length m as regularization parameters, optimzied depending on the side tensors. 
- `iter`: Number of iterations for the optimization process.
- `m_rate`: Missing rate (proportion of missing entries in the tensor).
  
### Running the Code
1. Load the required datasets into the MATLAB environment.
2. Set the parameters for the tensor completion process, such as `scale`, `epsilon`, `iter`, and `m_rate`.
3. Call the `CTTC.m` function with appropriate arguments.
4. The code will return the completed tensor and evaluate the performance based on the residual error.

## Requirements

- MATLAB (version R2019b or later)
- Required MATLAB toolboxes: MATLAB tensortoolbox https://www.tensortoolbox.org/ for tensorial operations. One may use MATLAB ADMM 
   https://www.mathworks.com/help/mpc/ref/admm.html for optimization puprposes as well. 

## Data Subset
A small subset of the data used in this project is provided in the data folder for ease of access and testing purposes. This subset contains representative samples of the full dataset and can be used to quickly test the code without the need to download the entire dataset.

To access this subset, navigate to the data folder in the repository, where you'll find the following files:
- cell_sim.mat
- drug_sim_mat_dense_subset
- filtered_target_tensor_1_subset
- gene_sim_mat_normalized_subset

Please note that these files represent a small portion of the full dataset (due to space limitation) and may not reflect all the variations present in the complete dataset. For more infromation on the full dataset and the side tensors please refer to the instructions provided in the manuscript [1]. 

## Example Usage

```matlab
% Load main tensor to be completed
X = load("path/to/incomplete_tensor.mat");

% Load side tensors
gene_sim_tensor = load("path/to/side_tensor_1.mat");
drug_sim_tensor = load("path/to/side_tensor_2.mat");
cell_sim_tensor = load("path/to/side_tensor_3.mat");

% Set parameters
Side = {gene_sim_tensor, drug_sim_tensor, cell_sim_tensor};  % Example side info
scale = [09, 0.8, 0.9];  % Scaling factors to be tuned depending on the density and sparisity of side infromation
epsilon = [0.1, 0.1, 0.1];  % Regularization parameters, to be optimzied by a grid search (or other well-known approaches)
iter = 50;  % Number of iterations- for this experimernt even 10 iterations suffice. 
m_rate = 0.8;  % Missing data rate defiens a Boolean/Logical mask for error-evaluation purposes

% Call the CTTC function
completed_tensor = CTTC(X, Side, iter, m_rate);
```

Based on the size of th einput main tensor and side infromation, the code initializes local distances `V_i` and form the Mahalanobis distance function using the local distances and the core tensor `Z`, from which it calculates the disatnce-metric-learnign based similarities and predicts the missing entries. 


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
```
For additional works or specific references related to this code, you may also refer to the refrences above.

For questions or issues related to the code, please contact Maryam Bagherian at [maryambagherian@isu.deu].






