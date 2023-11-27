import argparse
import classifier_IO as IO
import FeatureSelectors as FS
import Classification
from sklearn.ensemble import RandomForestClassifier

# Passed Flags
parser = argparse.ArgumentParser()
parser.add_argument('--TrainingData','-TR')
parser.add_argument('--TestData','-TS')
parser.add_argument('--Output','-O')
parser.add_argument('--Lable_Column','-LC')
parser.add_argument('--Threads','-t')
parser.add_argument('--FilterFeatures','-ff')
parser.add_argument('--nFeatures','-nf')
parser.add_argument('--Nr_Kfold','-K')

args = parser.parse_args()
train_data_path = args.TrainingData  
test_data_path = args.TestData
out_path = args.Output
lables_collumn_trainig_data = args.Lable_Column
threads = int(args.Threads)
filter_file = args.FilterFeatures
nf = int(args.nFeatures)
kf = int(args.Nr_Kfold)

#filter = IO.ImportFilter.load_filter(filter_file=filter_file)

data = IO.InportData.load_training_data(file=train_data_path)

#signature = FS.FeatureSelection.run_feature_selection_multi(threads=threads,data=data,filter=filter,nr_jobs=3,n_features=nf,lables=lables_collumn_trainig_data)
with open('/home/raphael/Desktop/OV_analysis_new/INF_classification/Scripts/ICON7geneset.txt') as file:
    lines = file.readlines()
    signature = [line.rstrip() for line in lines]

datadf, data_x, lab_y , features_t = IO.DataHandling.reduce_to_signature(data=data,lables=lables_collumn_trainig_data,signature=signature)

cross_val_score = Classification.Classify.classify(x=datadf,y=lab_y,k_fold=kf,threads=threads,output=out_path)

with open(f'{out_path}Signature.txt','w') as signatuer_file:
    for item in signature:
        signatuer_file.write(f'{item}\n')
        
with open(f'{out_path}Metrics.txt','w') as metrics_file:
    metrics_file.write(f'{cross_val_score}')
