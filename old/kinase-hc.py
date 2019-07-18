#%% [markdown]
# ## Libraries
import pandas as pd
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
import plotly.plotly as py
import plotly.offline

#%% [markdown]
# ## preparing data to be processed
data_path = 'data/human_protein_distance_matrix.csv'
df = pd.read_csv(data_path, delimiter=',', index_col=0)

print(df.info())
print(df.shape)
print(df.head(8))

# converting redundant distance matrix to condensed matrix
condensed_matrix = ssd.squareform(df)


#%% [markdown]
# ## Hierarchical clustering of the condensed matrix
H = linkage(condensed_matrix, 'ward')
fig = plt.figure(figsize=(25, 10))

dn = dendrogram(H)