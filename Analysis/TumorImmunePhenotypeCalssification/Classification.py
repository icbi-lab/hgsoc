
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import auc
from sklearn.metrics import roc_curve
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import cross_val_predict


from sklearn.ensemble import RandomForestClassifier
import joblib

class Classify:

        def classify(x,y,k_fold,threads,output):
                features = x.columns
                #print(features)
                x = x.values
                clf = RandomForestClassifier(n_estimators=300,n_jobs=threads)
                cv = StratifiedKFold(n_splits=k_fold)
                #print(x)
                #print(y)
                model = clf.fit(x,y)
                scoreSkf = cross_val_score(model,x,y,cv=cv)

                importance = model.feature_importances_
                imp = (importance)
                df_imp = pd.DataFrame({'Feature': features, 'ImportanceScore':imp})
                df_imp.sort_values(['ImportanceScore'],ascending=False,inplace=True)
                df_imp.reset_index(inplace=True)
                df_imp.to_csv(f'{output}FeatureImportance.csv')
                # print(df_imp)



                tprs = []
                aucs = []
                mean_fpr = np.linspace(0, 1, 100)

                #fig, ax = plt.subplots()
                #for i, (train, test) in enumerate(cv.split(x, y)):
                #        clf.fit(x[train], y[train])
                #        viz = RocCurveDisplay.from_estimator(
                #                clf,
                #                x[test],
                #                y[test],
                #                name="ROC fold {}".format(i),
                #                alpha=0.3,
                #                lw=1,
                #                ax=ax,
                #        )
                #        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                #        interp_tpr[0] = 0.0
                #        tprs.append(interp_tpr)
                #        aucs.append(viz.roc_auc)

                #ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)


                #mean_tpr = np.mean(tprs, axis=0)
                #mean_tpr[-1] = 1.0
        
                #mean_auc = auc(mean_fpr, mean_tpr)
                #std_auc = np.std(aucs)
                #ax.plot(
                #        mean_fpr,
                #        mean_tpr,
                #        color="b",
                #        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
                #        lw=2,
                #        alpha=0.8,
                #)

                #std_tpr = np.std(tprs, axis=0)
                #tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
                #tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
                #ax.fill_between(
                #        mean_fpr,
                #        tprs_lower,
                #        tprs_upper,
                #        color="grey",
                #        alpha=0.2,
                #        label=r"$\pm$ 1 std. dev.",
                #)

                #ax.set(
                #        xlim=[-0.05, 1.05],
                #        ylim=[-0.05, 1.05],
                #        title="Receiver operating characteristic example",
                #)
                #ax.legend(loc="lower right")
                #ax.figure.savefig(f'{output}RandomForrest_roc.svg',bbox_inches = "tight")
                #plt.show()
                y_pred = cross_val_predict(clf,x,y,cv=3)
                ConfusionMatrixDisplay.from_predictions(y,y_pred)
                plt.show()


                file_save = f'{output}FinalizedClassifier.sav'
                print(list(features))
                model.feature_names = list(features)
                joblib.dump(model,file_save)
                return scoreSkf

                
        def classify_loo(x,y,k_fold,threads,output):

                #clf = RandomForestClassifier(n_estimators=300,n_jobs=threads)
                clf = LogisticRegression(random_state=0)
                kf = LeaveOneOut()

                x, y = x[y != 2], y[y != 2]


                #print(x)
                #print(len(x))
                #print(y)

                if kf.get_n_splits(x) == len(x):
                        print("They are the same length, splitting correct")
                else:
                        print("Something is wrong")
                all_y = []
                all_probs=[]
                it = 0
                for train, test in kf.split(x, y):
                        it = it+1
                        #print(it)
                        #print(x.iloc[train])
                        all_y.append(y[test])
                        all_probs.append(clf.fit(x.iloc[train], y[train]).predict_proba(x.iloc[test])[:,1])
                all_y = np.array(all_y)
                all_probs = np.array(all_probs)
                #print(all_y) #For validation 
                #print(all_probs) #For validation

                fpr, tpr, thresholds = roc_curve(all_y,all_probs)
                print(fpr, tpr, thresholds) #For validation
                roc_auc = auc(fpr, tpr)
                plt.figure(1, figsize=(12,6))
                plt.plot(fpr, tpr, lw=2, alpha=0.5, label='LOOCV ROC (AUC = %0.2f)' % (roc_auc))
                plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', label='Chance level', alpha=.8)
                plt.xlim([-0.05, 1.05])
                plt.ylim([-0.05, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Receiver operating characteristic example')
                plt.legend(loc="lower right")
                plt.grid()
                plt.savefig(f'{output}Roc_plot_loocv.png',bbox_inches = "tight")
                plt.close()

        
        def classify_single_var_loo(x,y,threads,output):

                #clf = RandomForestClassifier(n_estimators=300,n_jobs=threads)
                clf = LogisticRegression(random_state=0)
                kf = LeaveOneOut()

                var = list(x.columns)
                for i in var:
                        
                        X = pd.DataFrame()
                        X[f'{i}'] = x[f'{i}']
                        print(X)

                        X, y = X[y != 2], y[y != 2]

                        #print(x)
                        #print(len(x))
                        #print(y)

                        if kf.get_n_splits(X) == len(X):
                                print("They are the same length, splitting correct")
                        else:
                                print("Something is wrong")
                        all_y = []
                        all_probs=[]
                        it = 0
                        for train, test in kf.split(X, y):
                                it = it+1
                                #print(it)
                                #print(x.iloc[train])
                                all_y.append(y[test])
                                all_probs.append(clf.fit(X.iloc[train], y[train]).predict_proba(X.iloc[test])[:,1])
                        all_y = np.array(all_y)
                        all_probs = np.array(all_probs)
                        #print(all_y) #For validation 
                        #print(all_probs) #For validation

                        fpr, tpr, thresholds = roc_curve(all_y,all_probs)
                        print(fpr, tpr, thresholds) #For validation
                        roc_auc = auc(fpr, tpr)
                        plt.figure(1, figsize=(12,6))
                        plt.plot(fpr, tpr, lw=2, alpha=0.5, label='LOOCV ROC (AUC = %0.2f)' % (roc_auc))
                        plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', label='Chance level', alpha=.8)
                        plt.xlim([-0.05, 1.05])
                        plt.ylim([-0.05, 1.05])
                        plt.xlabel('False Positive Rate')
                        plt.ylabel('True Positive Rate')
                        plt.title(f'Receiver operating characteristic {i}')
                        plt.legend(loc="lower right")
                        plt.grid()
                        plt.savefig(f'{output}Roc_plot-{i}_loocv.png',bbox_inches = "tight")
                        plt.close()
                        
                