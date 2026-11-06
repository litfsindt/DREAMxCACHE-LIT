# PGK2 DREAMxCACHE Challenge – LIT-Team-SXB

This repository contains the workflow developed by **LIT-Team-SXB** for the PGK2 DREAMxCACHE challenge.

## Repository structure

### 0. Dataset preparation

Files related to the preparation of the input dataset used for model training :

- `selecting_actives.py`: Python script used to filter the 7M putative active compounds according to the filtering rules provided by the challenge for the DEL dataset
- `filter_PAINS_reactive_groups.txt`: OpenEye filter file used to remove PAINS and reactive functional groups from the selected active compounds from `selecting_actives.py` to create `actives.smi`
- `selecting_inactives.sh`: Shell script used to launch the random inactive selection process, accounting for the size of each Parquet file and ensuring that selected compounds are not part of the 7M putative active set
- `sample_random_opendel.py`: Python script used to randomly select inactive compounds, with sampling weighted according to the size of each Parquet file

### 1. Model

Files related to the Chemprop model training and prediction : 

- `Chemprop_env.yml`: Conda environment used to build and run the Chemprop models
- `dataset.csv`: Input dataset containing the SMILES and corresponding labels used for model training
- `Models/`: Directory containing the five models used for the challenge
- `model.sh`: Chemprop command-line script used to train the models
- `predictions_on_test_set.csv`: Consensus prediction results obtained from the five models on the test set

### 2. GOLD_Docking

- GOLD configuration file (`.conf`) used to rescore the top 10% of compounds selected based on the model predictions

### 3. Key_Interactions

Files related to the detection of key protein–ligand interactions

- `search_IChem_HOH_3_TYR_242.py`: Python script used to detect the key hydrogen-bond interactions involving HOH 3 and TYR242
- `ichem_ifp.cpython-312-x86_64-linux-gnu.so`: IChem Python API library used by the script

### 4. SQM

- `cuby.in`: Input configuration file for the Cuby4 semi-empirical quantum mechanics (SQM) calculations
