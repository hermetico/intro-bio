
# coding: utf-8

# In[1]:

## Import libraries
import matplotlib
import matplotlib.pyplot as plt
## remove this line when running script from terminal, keep it when running notebooks
get_ipython().magic(u'matplotlib inline')

import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
import os # to join pathrs, etc..


# In[2]:

## Define filenames and sutff
peaks_path = 'output'
labels_path = 'all_labels.txt'

# reads labels into a df
train_labels_df = pd.read_csv(labels_path, sep="\t")
## take care because there are duplicates!!!!!!!!!
train_labels_df = train_labels_df.drop_duplicates(subset='file').reset_index(drop=True)
## adding a column for the labels with the id of the label
labels = train_labels_df['candy'].unique()
train_labels_df['candy_id'] = np.where(labels == train_labels_df['candy'][:,None])[-1]


# In[3]:

print "Total of %i files with peaks"% train_labels_df['file'].shape[0]

peaks_files = [ os.path.join(peaks_path, file) for file in train_labels_df['file'] ]


# In[4]:


results_folder = 'results'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)


# # Part 2
# 
# ## Step 7
# 
# Create a peak alignment matrix with the new data

# In[5]:

def read_files( paths ):
    # reads all the files
    raw_peaks = [pd.read_csv(path, sep="\t") for path in paths ]
    # combines them all in one dataframe
    df = pd.concat(raw_peaks).reset_index(drop=True)
    # converts peak_name column into integers
    peaks = [ int(peak.replace('p','')) for peak in df['peak_name']] 
    df['peak_name'] = peaks
    return df

def align_peaks( data, peaks=90):
    KM = KMeans(n_clusters=peaks, random_state=0)
    KM.fit(data[['t', 'r']])
    data['cluster_id'] = KM.labels_
    return data
    
def create_train_matrix( data ):
    total_clusters = data['cluster_id'].max() + 1
    d = { colname: [0] * total_clusters for colname in data['measurement_name'].unique() }
    matrix = pd.DataFrame(data=d)
    ## fill the matrix
    for name in data['measurement_name'].unique():
        patient = data[data['measurement_name'] == name]
        clusters = patient['cluster_id']
        #for cluster in clusters:
        matrix[name][clusters] = 1
    
    return matrix.transpose()

raw_df = read_files( peaks_files )

aligned_df = align_peaks( raw_df )
train_matrix_df = create_train_matrix(aligned_df)

print "Showing head of the matrix"
train_matrix_df.head()


# ## Step 8
# 
# Train random forest classifier

# In[6]:

rf = RandomForestClassifier(n_estimators=20)
rf = rf.fit(train_matrix_df.as_matrix(), train_labels_df['candy'])


# In[7]:

prediction = rf.predict(train_matrix_df.as_matrix())


# In[8]:

confusion_matrix(train_labels_df['candy'], prediction)


# ## Step 9
# 
# Implement 5-fold cross validation
# 

# In[9]:

def k_cross_fold_validation( data, labels, k=5):
    """Returns the predictions and the indexes of each test slice"""
    rf = RandomForestClassifier(n_estimators=25)
    predictions, indexes = [], []
    size = data.shape[0]
    step_size = size / 5
    for i in range(k):
        indexes.append(range(i*k, min(i*k + k, size)))
        xtrain = np.array(data)
        ytrain = np.array(labels)

        # picks test slices
        xtest = xtrain[i*k:i*k + k]
        ytest = ytrain[i*k:i*k + k]

        # removes test slices from the training sets
        xtrain = np.delete(xtrain, np.s_[i*k:i*k + k], axis=0)
        ytrain = np.delete(ytrain, np.s_[i*k:i*k + k], axis=0)

        
        rf = rf.fit(xtrain, ytrain)
        prediction = rf.predict(xtest)

        predictions.append(prediction.tolist())
        
    return np.array(predictions), np.array(indexes)

pred, indexes = k_cross_fold_validation(train_matrix_df.as_matrix(), train_labels_df['candy'])


# ## Step 10 
# Report the mean accuracy, sensitivity and specificity

# In[10]:

def report_mean( preds, trues):
    flatten_pred = np.hstack(preds)
    results = np.zeros(trues.unique().shape[0])
    for i, label in enumerate(trues.unique()):
        # indexes for such label
        class_trues =  np.argwhere(trues == label).ravel()
        # predicted label on such indexes
        predicted = flatten_pred[class_trues]

        results[i] = np.count_nonzero(predicted == label) / float(len(class_trues))
    return results.mean()

mean_accuracy = report_mean(pred, train_labels_df['candy'])
print "Mean Accuracy: ", mean_accuracy


# In[11]:

def report_sens_spec( pred, trues):
    preds = np.hstack(pred)
    classes = trues.unique()
    # report sensitivity for each class
    for label in classes:
        print label.capitalize()
        # indexes for such label
        class_trues =  np.argwhere(trues == label).ravel()
        class_falses =  np.argwhere(trues != label).ravel()

        pred_trues = np.argwhere(preds == label).ravel()
        pred_falses = np.argwhere(preds != label).ravel()

        TP = sum([ 1. for p in pred_trues if p in class_trues ])
        TN = sum([ 1. for p in pred_falses if p in class_falses ])
        FN = sum([ 1. for p in pred_falses if p in class_trues ])
        
        sen = TP / (TP + FN)
        spec = TN / (TN + FN)
        print "Sensitivity: ", sen
        print "Specificity: ", spec
    
report_sens_spec(pred, train_labels_df['candy'])    


# ## Step 11
# Extract the five most discriminating features (peaks) by using the Gini index

# In[21]:

train_matrix_df.shape[1]


# In[ ]:



