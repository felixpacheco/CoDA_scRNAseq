#!/usr/bin/env python
# coding: utf-8

# Data analysis
import pandas as pd
import numpy as np

# Data visualization
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.cm as cm
import seaborn as sns

# Math
import math
from scipy import special
import scipy.stats as ss
from scipy.spatial.distance import pdist, cdist
from sklearn.metrics.pairwise import euclidean_distances
from scipy.stats import f_oneway
from scipy.stats import ttest_ind

print("loaded all dependencies")

def closure(composition, kappa):
    ''' Returns the composition closed to a kappa value '''
    counts = pd.DataFrame(index=composition.columns)
    for column in composition.T:
        counts[column] = kappa * composition.T[column]/np.sum(composition.T[column])
    return counts.T

def gmean(composition):
    ''' Returns de geometric mean of a composition '''
    gmean = ss.mstats.gmean(composition)

    return [100 * i / np.sum(gmean) for i in gmean]

def sample_center(composition):
    ''' Returns de sample center of a series of compositions'''

    return pow(composition/gmean(composition), 1./np.sqrt(totvar(composition)))

def var_matrix(composition):
    ''' Returns de variation matrix of a series of compositions'''

    reduc = np.array(composition)[:, :min(500, np.shape(composition)[1])]

    # New vectorized version. Faster than ketchup!
    var_matrix = np.var(np.log(reduc[:, :, None] * 1. / reduc[:, None]), axis=0)

    return var_matrix

def totvar(composition):
    ''' Returns the total variation of a series of compositions'''

    totvar = 1. / (2 * np.shape(var_matrix(composition))[0]) * np.sum(var_matrix(composition))

    return totvar


def aitchison_dist (x, y):
    '''Returns Aitchison distance'''
    D = len(x)
    sum_ij = 0
    for i in range(D):
        for j in range(D):
            sum_ij += (math.log(x[i]/x[j]) - math.log(y[i]/y[j])) ** 2
    return(np.sqrt (1/(2*D) * sum_ij) )


def _clr_internal(obj):
    return (np.log(obj.T) - np.mean(np.log(obj.T))).T


def _alr_internal(obj):
    return (np.log(obj.T/obj.T.loc[obj.columns[-1]])).T


def _ilr_internal(obj, psi):
    return pd.DataFrame(np.dot(_clr_internal(obj), psi.T), index=obj.index)


def sbp_basis(obj):
    ''' Define basis to use in IRL transformation '''
    dim = np.shape(obj)[1]
    psi = np.zeros([dim-1, dim])
    for i in range(dim-1):
        for j in range(dim):
            if j+1 <= dim-i-1:
                psi[i, j] = np.sqrt(1./((dim-i-1)*(dim-i)))
            elif j+1 == dim-i:
                psi[i, j] = -np.sqrt((dim-i-1)/(dim-i))
    check_basis(psi)

    return psi


def check_basis(psi):
    ''' Check if basis is orthonormal '''
    ident = np.matmul(psi, psi.T)

    if np.trace(ident) != np.shape(ident)[0]:
        raise AttributeError("Error: Basis is not normalized.")
    if np.abs(np.sum(ident-np.diag(np.diagonal(ident)))) > 1e-6:
        raise AttributeError("Error: Basis is not orthogonal.")



def bayesian_transform(composition, n_samples, kind, output):
    ''' A method to calculate the logratio transform  with Bayesian replacement'''
    logratio = pd.DataFrame(index=composition.columns)

    for column in composition.T:
        p_matrix = ss.dirichlet.rvs(composition.T[column]+0.5, n_samples)
        if kind == 'clr':
            c_matrix = _clr_internal(p_matrix)
        elif kind == 'alr':
            c_matrix = [np.log(i/i[-1]) for i in p_matrix]
        if output == 'mean':
            logratio[column] = [np.mean(i) for i in zip(*c_matrix)]
        elif output == 'std':
            logratio[column] = [np.std(i) for i in zip(*c_matrix)]
    return logratio.T

def clr(composition, bayesian=False, n_samples=5000):
    ''' Wrapper for CLR '''
    if not bayesian:
        return _clr_internal(composition)
    return bayesian_transform(n_samples, 'clr', 'mean')

def clr_std(composition, n_samples=5000):
    ''' Wrapper for CLR bayesian error estimate'''
    return bayesian_transform(n_samples, 'clr', 'std')

def alr(composition, part=None, bayesian=False, n_samples=5000):
    ''' Wrapper for ALR '''
    if part:
        parts = composition.index.tolist()
        parts.remove(part)
        composition.reindex(parts+[part])

    if not bayesian:
        return _alr_internal(composition)[composition.columns[:-1]]
    return bayesian_transform(n_samples, 'alr', 'mean')[composition.columns[:-1]]

def alr_std(composition, part=None, n_samples=5000):
    ''' Wrapper for ALR error estimate'''
    if part:
        parts = composition.index.tolist()
        parts.remove(part)
        composition.reindex(parts+[part])

    return bayesian_transform(n_samples, 'alr', 'std')[composition.columns[:-1]]


def ilr(composition, psi=None, bayesian=False, n_samples=5000):
    if not psi:
        psi = sbp_basis(composition)
    else:
        check_basis(psi)
    if not bayesian:
        return _ilr_internal(composition, psi)
    return np.dot(bayesian_transform(n_samples, 'clr', 'mean'), psi.T)


def zero_replacement(composition):
    n_samples = 5000
    counts = pd.DataFrame(index=composition.columns)

    for column in composition.T:
        p_matrix = ss.dirichlet.rvs(composition.T[column]+0.5, n_samples)
        counts[column] = [np.mean(i) for i in zip(*p_matrix)]

    return counts.T

def _svd(clr):
    scores, eig_val, loadings = np.linalg.svd(clr)
    scores = pd.DataFrame(scores.T[0:2, :], columns=clr.index, index=['pc1', 'pc2'])
    loadings = pd.DataFrame(np.inner(eig_val*np.identity(len(eig_val)),
                                     loadings.T[0:len(eig_val), 0:len(eig_val)])[0:2],
                            columns=clr.columns[0:len(eig_val)], index=['pc1', 'pc2'])
    return scores, eig_val, loadings

print("all functions loaded")  

# READ DATA
counts = pd.read_table("../data/bulk/counts_bulk_patient.tsv")

print("Data loaded")

counts.set_index('Ensembl_gene_id', inplace=True)
counts = counts[(counts.T != 0).any()]

print(counts.head(5))


counts = zero_replacement(counts)
print(counts.head(5))

clr_data = clr(counts)
data_array = np.array(clr_data)
print(clr_data)

clr_data.to_csv('clr.csv', index = True)

np.savetxt("clr_array.csv", data_array, delimiter=",")
