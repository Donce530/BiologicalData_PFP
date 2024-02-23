# BiologicalData_PFP
This repo contains all the files and notebooks required to fulfil the requirements of the Protein Function Prediction Project for the Biological Data course at UNIPD.

To reproduce the results:

First run the pre-processing notebook on the files in the biological_data_pfp.zip folder provided in the project specification. This will produce 6 dataframes, `train` and `labels` for each of the 3 GO ontologies, within a dictionary to be used for supervised training.

Next, the training dictionary can then be used in all of the architecture notebooks (FNN, FNN3, CNN, ResNet1D) to produce training metrics for the model evaluation.

Finally, the best model is loaded into the BD_Final_Project.ipynb which performs predictions on the test set and out puts the predictions file for the CAFA evaluator. This notebook also contains the code and computations for graphs that are used in the final report.
