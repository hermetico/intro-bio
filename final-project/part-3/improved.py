
# coding: utf-8

# In[1]:

## Import libraries
import math
import matplotlib
import matplotlib.pyplot as plt
## remove this line when running script from terminal, keep it when running notebooks
get_ipython().magic(u'matplotlib inline')

import numpy as np
import pandas as pd

import os # to join pathrs, etc..
import pickle # to store the models

from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier


# In[2]:

results_folder = 'results'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)

peaks_path = "output"


# In[3]:

def read_model( path ):
    """Helper function which loads models from pickle python files"""
    path = os.path.join(models_folder,path + '.pck' )
    with open( path  , 'r') as f:
        model = pickle.load(f)
    return model

def read_files( paths ):
    # reads all the files
    raw_peaks = [pd.read_csv(path, sep="\t") for path in paths ]
    # combines them all in one dataframe
    df = pd.concat(raw_peaks).reset_index(drop=True)
    # converts peak_name column into integers
    peaks = [ int(peak.replace('p','')) for peak in df['peak_name']] 
    df['peak_name'] = peaks
    return df

def create_train_matrix( data , n_clusters):
    total_clusters = n_clusters
    d = { colname: [0] * total_clusters for colname in data['measurement_name'].unique() }
    matrix = pd.DataFrame(data=d)
    ## fill the matrix
    for name in data['measurement_name'].unique():
        patient = data[data['measurement_name'] == name]
        clusters = patient['cluster_id']
        #for cluster in clusters:
        matrix[name][clusters] = 1
    
    return matrix.transpose()

def write_csv( names, labels, filename):
    path = os.path.join(results_folder, filename +".csv")
    with open(path, "w") as f:
        f.write("file,candy\n")
        for name, label, in zip(names, labels):
            f.write("%s,%s\n"%(name, label))
            


# In[4]:

#training kmeans again with all the data files of 3 parts
peaks_files =  [
    os.path.join(peaks_path, file) 
    for file in os.listdir(peaks_path)
        if file.endswith('csv') 
] + [
    os.path.join("../part-2/output", file) 
    for file in  os.listdir("../part-2/output")
        if file.endswith('csv') 
]# + [
#    os.path.join("../part-1/output", file) 
#    for file in  os.listdir("../part-1/output")
#        if file.endswith('csv') 
#] 

def align_peaks( data, peaks=90):
    KM = KMeans(n_clusters=peaks, n_init=100)
    KM.fit(data[['t', 'r']])
    data['cluster_id'] = KM.labels_
    return KM, data  


print peaks_files


# In[7]:

all_df = read_files( peaks_files )
# for the number of clusters, we expect to have some in common and some not, for the different types of candies
# that is why we are using the average + 0.5
CLUSTERS = int(np.average([ all_df[all_df['measurement_name'] == name].shape[0] for name in all_df['measurement_name'].unique()]) * 1.5)
print "Guessing for  ", CLUSTERS, " clusters"


# In[8]:

k_means, aligned_all_df = align_peaks( all_df ,peaks = CLUSTERS)
KM = k_means


# In[9]:

# training random forests with all the data with the previous 2 parts, and testing with the data of the part 3
train_files = [ file for file in  os.listdir(peaks_path) if file.endswith('csv') ] 
#aligned_all_df[[aligned_all_df['measurement_name'] in train_files]]

# select only the elements that are in the list of files of the part 3
test_df = aligned_all_df[aligned_all_df['measurement_name'].isin(train_files)]

# now, select the elements that are not in that list
train_df = aligned_all_df[~aligned_all_df['measurement_name'].isin(train_files)]


# In[10]:

# create matrices for both
train_matrix_df = create_train_matrix(train_df, KM.n_clusters)
test_matrix_df = create_train_matrix(test_df, KM.n_clusters)


# In[12]:

# now, create the list of labels for the training
# the first part labels need some cleaning
#part_1_df = pd.read_csv("../part-1/files categories.txt", sep=",")
#part_1_df = part_1_df.apply(lambda x: x.str.strip())
#part_1_df = part_1_df.rename(columns=lambda x: x.strip())
#part_1_df['file'] = part_1_df['file'].apply(lambda x: x +'_ims.csv')
#part_1_df['candy'][part_1_df['candy'] == 'halls_citruzzz'] = 'citrus'

# load labels second part
part_2_df = pd.read_csv("../part-2/all_labels.txt", sep="\t")

# merge labels
#all_labels = pd.concat([part_1_df, part_2_df]).reset_index(drop=True)
all_labels = part_2_df
# set file column as index
labels_df = all_labels[['file', 'candy']].set_index('file')


# In[13]:

labels_df


# In[16]:

train_matrix_df


# In[17]:

# now train the random forest
RF = RandomForestClassifier()
RF = RF.fit(train_matrix_df.as_matrix(), labels_df['candy'])


# In[18]:

# apply the trained RF to the test data
prediction = RF.predict(test_matrix_df.as_matrix())
print prediction


# In[19]:

write_csv(test_matrix_df.index, prediction, "random_forests_prediction")


# In[ ]:



