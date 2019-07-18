import pandas as pd
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
import plotly.plotly as py
import plotly.offline

data_path = 'data/phosphorylation_data.txt'
df = pd.read_csv(data_path, delimiter=',', index_col=0)

print(df.info())
print(df.shape)
print(df.head(8))
