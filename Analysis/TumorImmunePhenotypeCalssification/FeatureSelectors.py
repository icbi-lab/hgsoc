import pandas as pd
import numpy as np

from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier

from joblib import Parallel, delayed

class Selectors:

    def recursive_feature_elemination(data,n_features,lables,filter,clf):
        
        features = filter
        data1 = data.copy()

        while len(features) > n_features:
            data = data1[features]
            lab = data1[lables]
            data = data.fillna(0)
            x ,y = data.values, lab.values.ravel()
            model = clf.fit(x,y) 
            importance = model.feature_importances_
            imp = (importance)
            df_imp = pd.DataFrame({'Feature': features, 'ImportanceScore':imp})
            df_imp.sort_values(['ImportanceScore'],ascending=False,inplace=True)
            df_imp.reset_index(inplace=True)
            df_imp = df_imp[(df_imp['ImportanceScore'] > 0)]
            Idx_min = df_imp['ImportanceScore'].idxmin()
            df_imp = df_imp.drop(Idx_min)
            features = list(df_imp['Feature'])
        
        return features

class FeatureSelection:

    def run_feature_selection_multi(threads,data,filter,nr_jobs,n_features,lables):

        f_filter = list(data.columns)
        filterpreliminary = filter
        filter = []
        for i in filterpreliminary:
            if i in f_filter:
                filter.append(i)
            else:
                continue

        classifierList = [
            RandomForestClassifier(n_estimators=1000, n_jobs=threads),
            #ExtraTreesClassifier(n_estimators=1000, n_jobs=threads),
            GradientBoostingClassifier(n_estimators=100),
            AdaBoostClassifier(n_estimators=100)
        ]

        FeaturesList = []
       
        FeaturesList.append(Parallel(n_jobs=nr_jobs)(delayed(Selectors.recursive_feature_elemination)(clf=i,data=data,n_features=n_features,lables=lables,filter=filter) for i in classifierList))
        FeaturesL = []
        for i in FeaturesList:
            for x in i:
                FeaturesL.extend(x)

        diff = list({x for x in FeaturesL if FeaturesL.count(x) > 1})

        return diff
    
    def run_feature_selection_single(threads,data,filter,n_features,lables):
        
        f_filter = list(data.columns)
        filterpreliminary = filter
        filter = []
        for i in filterpreliminary:
            if i in f_filter:
                filter.append(i)
            else:
                continue

        classifier = RandomForestClassifier(n_estimators=1000, n_jobs=threads)

      
        features = Selectors.recursive_feature_elemination(clf=classifier,data=data,n_features=n_features,lables=lables,filter=filter)

        return features