#!/usr/bin/python
import matplotlib.pyplot as plt
import logging
import numpy as np
from itertools import combinations
from scipy.io import loadmat

from src.optknock import OptKnock
from src.models import *
from src.analysis_toolbox import *
from src.html_writer import HtmlWriter


def main():
    #draw('wild_type', None, None)
    draw('rubisco', dimension=2)
    draw('SBPase', dimension=2)
    draw('SBPase_with_rbc', dimension=2)
    #draw('deoxyribose')
    #draw('MOG')

def draw(title, dimension=2):
    if dimension == 3:
        raise ValueError('currently only 2D heatmaps can be visualized')

    data = loadmat('res/data_%dD_%s.mat' % (dimension, title))

    fig = plt.figure()
    yield_data = data['yields']
    slope_data = data['slopes']
    carbon_sources = list(data['cs'].flat)
    single_ko_list = list(data['ko'].flat)

    if 'electrons' in carbon_sources:
        i = carbon_sources.index('electrons')
        yield_data = np.delete(yield_data, i, 1)
        slope_data = np.delete(slope_data, i, 1)
        carbon_sources.pop(i)
    
    # reshape the slope matrix into a block matrix, where 
    # the x and y of each block are the first and second knockouts
    # and all carbon-sources are in inside the square block

    # calculate the number of columns/rows in the square 2D matrix 
    # which contains the heatmap data.
    N_ko = len(single_ko_list)
    N_cs = len(carbon_sources)
    N_block = int(np.ceil(np.sqrt(N_cs)))
    N_full = N_ko * N_block
    mat = np.zeros((N_full, N_full))

    labels = [''] * N_full
    for i, ko in enumerate(single_ko_list):
        j = i * N_block
        labels[j] = ko

    for i in xrange(N_ko):
        for j in xrange(N_ko):
            for k in xrange(N_cs):
                x = int(i*N_block + np.floor(k/N_block))
                y = int(j*N_block + (k % N_block))
                mat[x,y] = slope_data[i + j*N_ko, k]
                
    mat[np.where(np.isnan(mat))] = 0
    plt.pcolor(mat, figure=fig)
    plt.colorbar()
    fig.savefig('res/heatmap_%dD_%s.svg' % (dimension, title))
    
if __name__ == "__main__":
    main()
