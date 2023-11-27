import pandas as pd

class InportData:

    def load_training_data(file):

        training_data = pd.read_csv(file)

        return training_data
    
    def load_test_Data(file):

        if file != None:
            test_data = pd.read_csv(file)
        else:
            test_data = 'split_training_data'

        return test_data
    
class ExportData:

        def write_results_tables(table,file):

            table.to_csv(f'{file}')

class ImportFilter:

    def load_filter(filter_file):

        filter_df = pd.read_csv(filter_file,header=None)
        filter_f = list(filter_df[0])
        return filter_f


class   DataHandling:

    def reduce_to_signature(data,signature,lables):
        lab = data[lables]
        l = lab.unique()
        d = {}
        x = 0
        for i in l:
            d[i]=x
            x = x+1
        lab = lab.replace(d)
        filt = list(data.columns)
        filterfinal = []
        for i in signature: 
            if i in filt:
                filterfinal.append(i)
            else:
                continue
        datadf = data[filterfinal]
        Features = datadf.columns
        datadf = datadf.fillna(0)
        return datadf, datadf.values, lab.values.ravel(), Features


    def reduce_to_filter(data,filter,lables):

        lab = data[lables]
        l = lab.unique()
        d = {}
        x = 0
        for i in l:
            d[i]=x
            x = x+1
        lab = lab.replace(d)
        filt = list(data.columns)
        filterfinal = []
        for i in filter: 
            if i in filt:
                filterfinal.append(i)
            else:
                continue
        datadf = data[filterfinal]
        Features = datadf.columns
        datadf = datadf.fillna(0)
        return datadf, datadf.values, lab.values.ravel(), Features

