import os
import sys
import numpy as np
import pandas as pd
import torch
from loguru import logger
import pickle
import random
from model_MLP import *

# Parameters
random_state = 42
np.random.seed(random_state)
random.seed(random_state)

# Check available device
device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
logger.info(f"Device: {device}")


###----------------------------------------------------------------------------------------------------------------------------------



#usage: python 2_train_MLP_nestedCV.py ctranspath
if len(sys.argv) != 4:
    raise ValueError("Usage: python 2_train_MLP_nestedCV.py <feature_type> <LOO_fold_N> <train/test>")
feature_type = sys.argv[1]
LOO_fold_N = int(sys.argv[2]) 
train_or_test = sys.argv[3]  # "train" or "test"


if feature_type not in ["uni","ctranspath", "musk", "resnet50"]:
    raise ValueError(f"Invalid feature type: {feature_type}")
    exit()
if train_or_test not in ["train", "test"]:
    raise ValueError(f"Invalid option for train/test: {train_or_test}")
    exit()

###----------------------------------------------------------------------------------------------------------------------------------
feature_type_dims = {"uni": 1536, "ctranspath": 768, "musk": 2048, "resnet50": 2048}

###----------------------------------------------------------------------------------------------------------------------------------
gene_file = f"01_ST_QC/01.gene_sd.filtered.csv"
gene_names = pd.read_csv(gene_file, index_col=0, header=None).index.values

# Model parameters
n_inputs = int(feature_type_dims[feature_type])
n_hiddens = int(feature_type_dims[feature_type])
num_outputs = len(gene_names)
dropout = 0.2
K_fold = 7
lr = 0.001 ###try a higher lr 0.001?
patience = 50 #?
num_epochs = 500



    
###----------------------------------------------------------------------------------------------------------------------------------


result_dir = f"01.nestedCV_MLP_K{K_fold}_{feature_type}_lr{lr}_patience{patience}/"

loss_fn = torch.nn.MSELoss()

class MaskedMSELoss(torch.nn.Module):
    def __init__(self):
        super(MaskedMSELoss, self).__init__()
        self.mse_loss = torch.nn.MSELoss(reduction='none')  # Use 'none' to calculate element-wise loss

    def forward(self, predictions, targets):
        mask = (targets != 0).float()  # Create a mask where labels are not 0
        loss = self.mse_loss(predictions, targets)  # Compute element-wise MSE loss
        masked_loss = loss * mask  # Apply the mask
        return masked_loss.sum() / mask.sum()  # Normalize by the number of non-zero labels

# Replace the loss function in your code
#try_loss_wfmse = False
#loss_fn = MaskedMSELoss()
#result_dir = f"01_all_cancers_MLP_{feature_type}_masked_mse/"
#if try_loss_wfmse:
#    from loss import *
#    loss_fn = weighted_focal_mse_loss
#    result_dir = f"01_all_cancers_MLP_{feature_type}_wfmse/"


if not os.path.exists(result_dir):
    os.makedirs(result_dir)

#clean log file
if os.path.exists(f"{result_dir}/loo{LOO_fold_N}.log.training"):
    os.remove(f"{result_dir}/loo{LOO_fold_N}.log.training")

logger.add(f"{result_dir}/loo{LOO_fold_N}.log.training", format="{time} {level} {message}", level="INFO")
logger.add(sys.stdout, format="{time} {level} {message}", level="INFO")

logger.info(f"Train using {feature_type} features, {K_fold} folds, {n_inputs} = inputs, {n_hiddens} = hiddens, \
dropout = {dropout}, {num_outputs} = outputs\n\
lr = {lr}, num_epochs = {num_epochs}, patience = {patience}\n\
output will be saved in {result_dir}")

###----------------------------------------------------------------------------------------------------------------------------------




gene_file = f"01_ST_QC/01.gene_sd.filtered.csv"
gene_names = pd.read_csv(gene_file, index_col=0, header=None).index.values
num_outputs = len(gene_names)
logger.info(f"Number of genes: {num_outputs}")

path2features = f"features/"
path2rna = f"01_ST_QC/"

#01_ST_QC/BH23.logST.filtered.pkl
#sample_QC_file = f"01_ST_QC/01.sample_summary.QC.csv"
#metadata_file = "Data_locations_v20250516.txt"



slide_train = np.loadtxt(f"01_slide_train.txt", dtype= "str")
logger.info(f"Slide train: {slide_train}")
slide_test = np.loadtxt(f"01_slide_test.txt", dtype= "str")

#slide_train names before "_" or "PHE"
slide_train_patient = [slide.split("_")[0] for slide in slide_train]
slide_names_patient_dict = dict(zip(slide_train, slide_train_patient))

slide_test_patient = [slide.split("_")[0] for slide in slide_test]
sample_names_patient = [slide_names_patient_dict[sample] for sample in slide_train]
sample_names_patient_unique = sorted(list(set(sample_names_patient)))


#LOO_fold_N must be less or equal to the number of training samples
if LOO_fold_N > len(sample_names_patient_unique) or LOO_fold_N < 1:
    raise ValueError(f"LOO_fold_N {LOO_fold_N} must be 1:{len(sample_names_patient_unique)}.\
                     Number of unique patients in training set: {len(sample_names_patient_unique)}")


logger.info(f"Slide train: {slide_train}.\nSlide train size: {len(slide_train)}.\n")



###----------------------------------------------------------------------------------------------------------------------------------
if train_or_test == "train":
    logger.info(f"Starting nested cross-validation with LOO fold {LOO_fold_N} on training samples: {slide_train}, N train patients: {len(sample_names_patient_unique)}")
    nested_cross_validate(
        sample_names=slide_train,
        slide_names_patient_dict =slide_names_patient_dict,
        loo_fold_k = LOO_fold_N,
        model_class=MLPRegression,   
        optimizer_class=torch.optim.Adam,
        gene_names=gene_names,
        criterion= loss_fn,    
        device=device,
        path2features=path2features,
        path2selected_barcode = "selected_barcode/",
        path2rna=path2rna,
        feature_type=feature_type,
        #batch_size= 32,
        n_inputs=n_inputs,
        n_hiddens=n_hiddens,
        n_outputs=num_outputs,
        dropout=dropout,
        lr=lr,
        num_epochs=num_epochs,
        patience=patience,
        k_folds= K_fold,
        result_dir=result_dir,
        random_state=random_state
    )
elif train_or_test == "test":
    ''' After all folds are done, evaluate the test set
    '''    
    logger.info(f"Evaluating test samples: {slide_test}, N train patients: {len(sample_names_patient_unique)}")
    evaluate_final_test(len(sample_names_patient_unique), slide_test, model_class = MLPRegression, \
                        result_dir = result_dir, output_dir = f"{result_dir}/testSet/", \
                        device = device, path2rna = path2rna, path2features = path2features, \
                        path2selected_barcode = "selected_barcode/", 
                        gene_names = gene_names, feature_type = feature_type, \
                        n_inputs = n_inputs, n_hiddens = n_hiddens, n_outputs = num_outputs, \
                        dropout = dropout, k_folds = K_fold)


exit()

# Perform cross-validation
cross_val_results = cross_validate(
    sample_names=slide_train,
    model_class=MLPRegression,   
    optimizer_class=torch.optim.Adam,
    gene_names=gene_names,
    criterion= loss_fn,    
    device=device,
    path2features=path2features,
    path2rna=path2rna,
    feature_type=feature_type,
    batch_size= 32,
    n_inputs=n_inputs,
    n_hiddens=n_hiddens,
    n_outputs=num_outputs,
    dropout=dropout,
    lr=lr,
    num_epochs=num_epochs,
    patience=patience,
    k_folds=K_fold,
    result_dir=result_dir,
    random_state=random_state
)

# Access cross-validation results
avg_valid_loss = cross_val_results["avg_valid_loss"]
logger.info(f"Average Test Loss: {avg_valid_loss:.4f}")


# Evaluate the models on the test set
test_performance = evaluate_test_samples(
    sample_names=slide_test,
    model_class=MLPRegression,
    result_dir=result_dir,
    device=device,
    path2rna=path2rna,
    path2features=path2features,
    gene_names=gene_names,
    feature_type=feature_type,
    n_inputs=n_inputs,
    n_hiddens=n_hiddens,
    n_outputs=num_outputs,
    dropout=dropout,
    k_folds=K_fold
)

# Log the overall performance
logger.info(f"Test Performance Across Samples and Folds: {test_performance}")
