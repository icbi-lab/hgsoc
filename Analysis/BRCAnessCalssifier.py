# importing required modules
import pandas as pd
import scanpy as sc
import anndata 
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_validate
from sklearn.metrics import auc
from sklearn.metrics import roc_curve
import numpy as np
import joblib
from joblib import Parallel, delayed
from scipy import interp
from itertools import cycle
import matplotlib.pyplot as plt
import statistics
import random

# Preprocessing scRNAseq data
class Preprocessing:
    # Preprocesses single-cell RNA sequencing data using input files.
    def preprocess__sc(sc_h5, annotations):
        # Read the .h5ad file and annotations, merge them, and preprocess the data.
        adata = sc.read_h5ad(sc_h5)
        obs = pd.read_csv(annotations, sep='\t', index_col=0)
        adata.obs = pd.merge(adata.obs, obs, left_index=True, right_index=True, how='left')
        adata = adata[adata.obs['sample'].notna()]
        adata.obs['group'] = '1'
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        # Filter out specific cells based on sorting parameters.
        adata_cd45_neg = adata[adata.obs['sort_parameters'] == 'singlet, live, CD45-']
        adata_cancer_cells = adata_cd45_neg[adata_cd45_neg.obs['cell_type'] == 'Ovarian.cancer.cell']
        
        return adata, adata_cd45_neg, adata_cancer_cells
    
    # Preprocesses TCGA bulk RNAseq data using input files.
    def preprocess_tcga(tcga_counts, tcga_annotations):
        # Read the TCGA counts and annotations, merge them, and preprocess the data.
        adata_tcga = pd.read_csv(tcga_counts, sep='\t', index_col=0)
        adata_tcga = anndata.AnnData(X=adata_tcga.transpose())
        obs_df = pd.read_csv(tcga_annotations, index_col=0)
        adata_tcga.obs = pd.merge(adata_tcga.obs, obs_df, left_index=True, right_index=True, how='left')
        adata_tcga = adata_tcga[adata_tcga.obs['BRCAness'].notna()]
        adata_tcga.obs['groups'] = '2'
        sc.pp.normalize_total(adata_tcga, target_sum=1e6)
        sc.pp.log1p(adata_tcga)
        
        return adata_tcga
    
    # Combines and prepares a training dataset from single-cell and TCGA data.
    def get_training_dataset(adata_cancer, adata_tcga, tcga_annotations, coding_genes):
        # Filters and combines data to create a training set for later analysis.
        adata_cancer_adnexa = adata_cancer[adata_cancer.obs['tumor_supersite'] == 'Adnexa']
        obs_df = pd.read_csv(tcga_annotations, index_col=0)
        genes_cancer_cells = list(adata_cancer_adnexa.var_names)
        genes = list(set(genes_cancer_cells).intersection(list(adata_tcga.var_names)))
        tcga_expression = sc.get.obs_df(adata_tcga, keys=genes)
        df_training = pd.merge(tcga_expression, obs_df, left_index=True, right_index=True, how='right')
        genes.append('BRCAness')
        df_training = df_training[genes]
        filter_df = pd.read_csv(coding_genes, header=None)
        filter_f = list(filter_df[0])
        
        return df_training, filter_f

# Class for feature selection using recursive feature elimination.
class Selectors:
    # Uses recursive feature elimination for feature selection with a classifier.
    def recursive_feature_elemination(data, n_features, lables, filter, clf):
        # Iteratively fits a model and removes the least important feature until n_features remain.
        features = filter
        data1 = data.copy()

        while len(features) > n_features:
            data = data1[features]
            lab = data1[lables]
            data = data.fillna(0)
            x, y = data.values, lab.values.ravel()
            model = clf.fit(x, y)
            importance = model.feature_importances_
            imp = (importance)
            df_imp = pd.DataFrame({'Feature': features, 'ImportanceScore': imp})
            df_imp.sort_values(['ImportanceScore'], ascending=False, inplace=True)
            df_imp.reset_index(inplace=True)
            df_imp = df_imp[df_imp['ImportanceScore'] > 0]
            Idx_min = df_imp['ImportanceScore'].idxmin()
            df_imp = df_imp.drop(Idx_min)
            features = list(df_imp['Feature'])
        
        return features

# Class for feature selection across multiple classifiers.
class FeatureSelection:
    # Runs feature selection across multiple classifiers in parallel.
    def run_feature_selection_multi(threads, data, filter, nr_jobs, n_features, lables):
        # Applies feature selection across different classifiers and returns common features.
        f_filter = list(data.columns)
        filterpreliminary = filter
        filter = [i for i in filterpreliminary if i in f_filter]

        classifierList = [
            RandomForestClassifier(n_estimators=1000, n_jobs=threads),
            GradientBoostingClassifier(n_estimators=100),
            AdaBoostClassifier(n_estimators=100)
        ]

        FeaturesList = [Parallel(n_jobs=nr_jobs)(delayed(Selectors.recursive_feature_elemination)(clf=i, data=data, n_features=n_features, lables=lables, filter=filter) for i in classifierList)]

        FeaturesL = []
        for i in FeaturesList:
            for x in i:
                FeaturesL.extend(x)

        diff = list({x for x in FeaturesL if FeaturesL.count(x) > 1})

        return diff
    
    # Reduces the dataset to a signature of selected features.
    def reduce_to_signature(data, signature, lables):
        # Maps labels to numerical values and filters the dataset to the given signature.
        lab = data[lables]
        lab_encoding = {label: idx for idx, label in enumerate(lab.unique())}
        lab = lab.replace(lab_encoding)
        filterfinal = [i for i in signature if i in data.columns]
        datadf = data[filterfinal].fillna(0)
        return datadf, datadf.values, lab.values.ravel(), datadf.columns


class classification:
    
    # Method to perform classification using RandomForestClassifier
    def classify(x, y, k_fold, threads, output):
                # Extracting feature names from the input DataFrame
                features = x.columns
                
                # Convert DataFrame to numpy array for processing by scikit-learn
                x = x.values
                
                # Initialize the RandomForestClassifier with the specified number of trees and threads
                clf = RandomForestClassifier(n_estimators=300, n_jobs=threads)
                
                # Set up cross-validation with stratified k-folds to maintain the proportion of classes in each fold
                cv = StratifiedKFold(n_splits=k_fold)
                
                # Fit the model on the entire dataset to calculate feature importances
                model = clf.fit(x, y)
                
                # Perform cross-validation and return the scores for each fold
                scoreSkf = cross_val_score(model, x, y, cv=cv)

                # Calculate feature importances and store them in a DataFrame
                importance = model.feature_importances_
                df_imp = pd.DataFrame({'Feature': features, 'ImportanceScore': importance})
                
                # Sort the features by importance and reset the index
                df_imp.sort_values(['ImportanceScore'], ascending=False, inplace=True)
                df_imp.reset_index(inplace=True)
                
                # Save the feature importances to a CSV file
                df_imp.to_csv(f'{output}FeatureImportance.csv')

                # Prepare lists to store true positive rates and area under curve scores
                tprs = []
                aucs = []
                mean_fpr = np.linspace(0, 1, 100)  # Mean false positive rate

                # Create a plot for ROC curves
                fig, ax = plt.subplots()
                
                # Perform cross-validation and plot ROC curve for each fold
                for i, (train, test) in enumerate(cv.split(x, y)):
                    clf.fit(x[train], y[train])
                    viz = RocCurveDisplay.from_estimator(
                        clf,
                        x[test],
                        y[test],
                        name="ROC fold {}".format(i),
                        alpha=0.3,
                        lw=1,
                        ax=ax,
                    )
                    
                    # Interpolate true positive rates and store them along with AUC scores
                    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                    interp_tpr[0] = 0.0
                    tprs.append(interp_tpr)
                    aucs.append(viz.roc_auc)

                # Plot the chance line (diagonal)
                ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

                # Calculate the mean and standard deviation of the true positive rates
                mean_tpr = np.mean(tprs, axis=0)
                mean_tpr[-1] = 1.0
                mean_auc = auc(mean_fpr, mean_tpr)
                std_auc = np.std(aucs)
                
                # Plot the mean ROC curve
                ax.plot(
                    mean_fpr,
                    mean_tpr,
                    color="b",
                    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
                    lw=2,
                    alpha=0.8,
                )

                # Fill the area between the mean + std and mean - std of true positive rates
                std_tpr = np.std(tprs, axis=0)
                tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
                tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
                ax.fill_between(
                    mean_fpr,
                    tprs_lower,
                    tprs_upper,
                    color="grey",
                    alpha=0.2,
                    label=r"$\pm$ 1 std. dev.",
                )

                # Set plot limits and title, and show the legend
                ax.set(
                    xlim=[-0.05, 1.05],
                    ylim=[-0.05, 1.05],
                    title="Receiver operating characteristic example",
                )
                ax.legend(loc="lower right")
                
                # Save the ROC curve figure
                ax.figure.savefig(f'{output}RandomForrest_roc.svg', bbox_inches="tight")
                
                # Predict class labels for the input data using cross-validation
                y_pred = cross_val_predict(clf, x, y, cv=10)
                
                # Display and save the confusion matrix
                ConfusionMatrixDisplay.from_predictions(y, y_pred)
                plt.savefig(f'{output}RandomForrest_ConfusionMatrix.svg', bbox_inches="tight")
                plt.show()

                # Save the trained classifier to a file for later use
                file_save = f'{output}FinalizedClassifier.sav'
                feature_list = list(features)  # List of feature names
                model.feature_names = list(features)  # Add feature names to the model
                joblib.dump(model, file_save)
                
                # Return the cross-validation scores, trained model, and feature list
                return scoreSkf, model, feature_list


# List of seed values for random number generation to ensure reproducibility
seeds = list(range(0, 43))

# Iterate over each seed value
for i in seeds:
    # Set the seed for random number generation
    random.seed(i)

    # Define the paths to the data and annotation files
    tcga = 'featureCounts_for_DEseq2.tab'
    scRNA = 'GSE180661_matrix.h5'
    tcga_annotations = 'TCGA_Annotations.csv'
    scRNA_annotations = 'GSE180661_GEO_cells.tsv'
    coding_genes = 'Proteincoding_genes.tsv'

    # Preprocess single-cell RNA sequencing data
    adata, adata_cd45_neg, adata_cancer_cells = Preprocessing.preprocess_sc(sc_h5=scRNA, annotations=scRNA_annotations)

    # Preprocess TCGA bulk RNA sequencing data
    adata_tcga = Preprocessing.preprocess_tcga(tcga_counts=tcga, tcga_annotations=tcga_annotations)

    # Combine single-cell and TCGA data for training the classifier
    df_training, filter_coding = Preprocessing.get_training_dataset(adata_tcga=adata_tcga, adata_cancer=adata_cancer_cells, tcga_annotations=tcga_annotations, coding_genes=coding_genes)

    # Feature selection: identify the most informative genes
    signature = FeatureSelection.run_feature_selection_multi(threads=1, data=df_training, filter=filter_coding, nr_jobs=3, n_features=50, labels='BRCAness')

    # Reduce the dataset to the selected features (signature)
    datadf, data_x, lab_y, features_t = FeatureSelection.reduce_to_signature(data=df_training, labels='BRCAness', signature=signature)

    # Perform classification using the selected features and evaluate the model
    score, model_BRCAness, f_list = classification.classify(x=datadf, y=lab_y, k_fold=10, threads=5, output=f'BRCAness_Classifier/Seed_{i}')
