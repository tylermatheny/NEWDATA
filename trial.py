import pandas as pd
import matplotlib as plt
import numpy as np
import os
import statsmodels.api as sm
from scipy import stats
import seaborn as sns
spath = "/Users/tyler-matheny/NEWDATA/"
ids = pd.read_table('/Users/tyler-matheny/NEWDATA/STRESS/gene_exp.diff')
gene_id = ids['test_id']
gene_id = gene_id.to_frame()
for root, dirs, filenames in os.walk(spath):
    for f in filenames:
        if f.endswith('gene_exp.diff'):
            x = pd.read_table(os.path.join(root, f))
            x = x.rename(columns={ x.columns[9]: "logFC" })
            y = os.path.basename(root)
            x.columns = [str(col) + '_' + y for col in x.columns]
            x = x.rename(columns={ x.columns[0]: "test_id" })
            gene_id = gene_id.merge((x), on = 'test_id', how = 'outer')
cols = list(gene_id.columns.values)
vals = []
for i,y in enumerate(cols):
    if 'value' in y:
        vals.append(y)
    for z in vals:
        if 'p_value' in z:
            vals.remove(z)
        if 'q_value' in z:
            vals.remove(z)
for q in vals:
    x = 'gene_id.' + str(q)
    gene_id = gene_id.drop(gene_id[eval(x) < 1].index)
    merged_table = gene_id
sns_plot = sns.pairplot(merged_table)
sns_plot.savefig('output.png')
