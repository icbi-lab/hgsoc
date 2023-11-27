from ast import parse
import pandas as pd 
import joblib
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-M','--model',required=True)
parser.add_argument('-d','--data',required=True)
parser.add_argument('-S','--signature',required=True)
parser.add_argument('-O','--output',required=True)

args = parser.parse_args()
clf = args.model 
data = args.data 
sig = args.signature
out = args.output

model = joblib.load(clf)

df = pd.read_csv(data,index_col=0,sep=',')
#df = df.drop(columns='BRCAness')
df = df.transpose()
print(df)

with open(sig) as file:
    lines = file.readlines()
    signature = [line.rstrip() for line in lines]
print(df)
df = df[signature]

#print(signature)

#np.array(df)
x = df
#f_names = model.feature_names
#print(f_names)
y_predict = model.predict(x)

brcaness = list(y_predict)

print(brcaness)

df['TinfStatus'] = brcaness
df['TinfStatus'] = df['TinfStatus'].replace({0:'Infiltrated',1:'Desert',2:'Excluded'})

df.to_csv(out,index_label=0)