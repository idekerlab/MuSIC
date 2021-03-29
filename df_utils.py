import numpy as np
import pandas as pd
import sys
import os
from sklearn.metrics.pairwise import manhattan_distances, euclidean_distances, cosine_similarity
from scipy.spatial.distance import canberra

def upper_tri_values(df):
    ''' Return array with values of upper triangle of the DataFrame.
    
    Args:
        df: Symmetric DataFrame
    
    Return:
        Numpy array
    '''
    m = df.values
    return m[np.triu_indices(df.shape[0], k=1)]

def znorm(df):
    ''' Z-transform within each column.
    '''
    norm_df = pd.DataFrame(index=df.index, columns=df.columns)
    for c in df.columns:
        value = df[c]
        norm_df[c] = (value - value.mean()) / value.std()
    return norm_df


def cosine_similarity_scaled(df):
    ''' Calculate Cosine similarity between each pair of rows in a DataFrame.
        Similarity scaled into [0, 1]
    '''
    sim = cosine_similarity(df)
    shift = sim.min()
    sim -= shift
    scale = sim.max()
    sim /= scale
    return pd.DataFrame(sim, index=df.index.values, columns=df.index.values)


def manhattan_similarity(df):
    ''' Calculate Manhattan similarity between each pair of rows in a DataFrame.
        Similarity scaled into [0, 1]
    '''
    # Get manhattan distance
    dist = manhattan_distances(df)
    # Convert distance to similarity by max-minus
    sim = dist.max() - dist
    # Scale into [0,1]
    sim /= sim.max()
    return pd.DataFrame(sim, index=df.index.values, columns=df.index.values)

def euclidean_similarity(df):
    ''' Calculate Euclidean similarity between each pair of rows in a DataFrame.
        Similarity scaled into [0, 1]
    '''
    # Get euclidean distance
    dist = euclidean_distances(df)
    # Convert distance to similarity by max-minus
    sim = dist.max() - dist
    # Scale into [0,1]
    sim /= sim.max()
    return pd.DataFrame(sim, index=df.index.values, columns=df.index.values)

def canberra_similarity(df):
    ''' Calculate Canberra similarity between each pair of rows in a DataFrame.
        Similarity scaled into [0, 1]
    '''
    index = df.index.values
    dist = pd.DataFrame(0, index=index, columns=index, dtype=float)
    for i in range(len(index)-1):
        a = df.loc[index[i]].values
        for j in range(i+1, len(index)):
            b = df.loc[index[j]].values
            d = canberra(a, b)
            dist.at[index[i], index[j]] = d
            dist.at[index[j], index[i]] = d
    dist = dist.values
    # Convert distance to similarity by max-minus
    sim = dist.max() - dist
    # Scale into [0,1]
    sim /= sim.max()
    return pd.DataFrame(sim, index=index, columns=index)

def pearson_scaled(df):
    ''' Calculate Pearson correlation between each pair of rows in a DataFrame.
        Correlation scaled into [0, 1]
    '''
    corr = df.T.corr(method='pearson')
    shift = corr.min().min()
    corr -= shift
    scale = corr.max().max()
    corr /= scale
    return corr

def spearman_scaled(df):
    ''' Calculate Spearman correlation between each pair of rows in a DataFrame.
        Correlation scaled into [0, 1]
    '''
    corr = df.T.corr(method='spearman')
    shift = corr.min().min()
    corr -= shift
    scale = corr.max().max()
    corr /= scale
    return corr

def kendall_scaled(df):
    ''' Calculate Kendall correlation between each pair of rows in a DataFrame.
        Correlation scaled into [0, 1]
    '''
    corr = df.T.corr(method='kendall')
    shift = corr.min().min()
    corr -= shift
    scale = corr.max().max()
    corr /= scale
    return corr

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    ''' Check if the given numpy matrix is symmetric or not.
    '''
    return np.allclose(a, a.T, rtol=rtol, atol=atol)
