import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import statsmodels.api as sm
from scipy import stats
import seaborn as sns
import glob

def merge_cuffdiff(genefile):
    spath = os.getcwd()
    gene_id = pd.read_csv(genefile)
    for root, dirs, filenames in os.walk(spath):
        for f in filenames:
            if f.endswith('gene_exp.diff'):
                x = pd.read_table(os.path.join(root, f))
                x = x.rename(columns={ x.columns[9]: "logFC" })
                y = os.path.basename(root)
                x.columns = [str(col) + '_' + y for col in x.columns]
                x = x.rename(columns={ x.columns[0]: "gene_id" })
                gene_id = gene_id.merge((x), on = 'gene_id', how = 'left')
    return gene_id

def FPKM_restriction(c):
    vals = []
    cols = list(c.columns.values)
    for i,y in enumerate(cols):
        if 'value' in y:
            vals.append(y)
        for z in vals:
            if 'p_value' in z:
                vals.remove(z)
            if 'q_value' in z:
                vals.remove(z)
    for q in vals:
        c = c[c[q] >= 1]
        merged_table = c
    return merged_table

def split_gene_ID(x):
    new_IDs = x
    new = new_IDs["gene_id"].str.partition(".")
    new.columns = ['gene_id','.','decimal']
    new_IDs['gene_id'] = new['gene_id']
    return new_IDs

def merge_metrics(d):
    for file in glob.glob("METRICS/*.txt"):
        file = pd.read_table(file)
        file = file.rename(columns={ file.columns[0]: "gene_id" })
        d = d.merge((file), on = 'gene_id', how = 'left')
    return d

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
sg = pd.read_csv('sgtranscriptome.csv')
sg = sg.rename(index=str, columns={"test_id": "gene_id"})
sg = sg.dropna()
sg = sg.drop(labels='Mito', axis=1)
for file in glob.glob("METRICS/*.txt"):
    file = pd.read_table(file)
    file = file.rename(columns={ file.columns[0]: "gene_id" })
    sg = sg.merge((file), on = 'gene_id', how = 'left')
sg = sg.fillna(0)
metric_cols = sg.columns[13:-1]
X_metrics = sg[metric_cols]
y_FC = sg['Fold change']
y_Loc = sg['Localization']

X_train, X_test, y_train, y_test = train_test_split(X_metrics, y_Loc,
                                                   random_state = 0)
from sklearn.neighbors import KNeighborsClassifier
knn = KNeighborsClassifier(n_neighbors = 25)
knn.fit(X_train, y_train)
print('Accuracy of K-NN classifier on training set: {:.2f}'
     .format(knn.score(X_train, y_train)))
print('Accuracy of K-NN classifier on test set: {:.2f}'
     .format(knn.score(X_test, y_test)))

#kNN Regression
from sklearn.neighbors import KNeighborsRegressor

X_train, X_test, y_train, y_test = train_test_split(X_metrics, y_FC, random_state = 0)

knnreg = KNeighborsRegressor(n_neighbors = 5).fit(X_train, y_train)

print(knnreg.predict(X_test))
print('R-squared kNN test score: {:.3f}'
     .format(knnreg.score(X_test, y_test)))

X_train, X_test, y_train, y_test = train_test_split(X_metrics, y_FC,
                                                   random_state = 0)

linreg = LinearRegression().fit(X_train, y_train)
print('linear model intercept: {}'
     .format(linreg.intercept_))
print('linear model coeff:\n{}'
     .format(linreg.coef_))
print('R-squared score (training): {:.3f}'
     .format(linreg.score(X_train, y_train)))
print('R-squared score (test): {:.3f}'
     .format(linreg.score(X_test, y_test)))


est = sm.OLS(y_FC, X_metrics)
est2 = est.fit()
print(est2.summary())


sns.pairplot(X_metrics)
plt.show()
