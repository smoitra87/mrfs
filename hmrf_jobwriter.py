#! /usr/bin/env python


from itertools import product, repeat,chain,combinations
import os
import sys
import pickle
from operator import itemgetter
from os.path import split, splitext
import glob
import numpy as np
from collections import defaultdict,Counter
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

def plot_matrix(adj):
    # Plot the adj matrix
    plt.ion()
    plt.imshow(adj, interpolation='nearest')
    raw_input('Continue?')
    plt.close()


def create_adj_matrix(nVisNodes, nHidNodes, options):
    nNodes = nVisNodes + nHidNodes;
    adj = np.zeros((nNodes, nNodes))

    if 'linear_vis' in options:
        for idx in range(nVisNodes-1):
            adj[idx,idx+1] = 1

    if 'linear_hid' in options:
        for idx in range(nVisNodes, nVisNodes + nHidNodes-1):
            adj[idx,idx+1] = 1

    if 'vis_hid' in options:
        for idx in range(nVisNodes):
            adj[idx,idx+nVisNodes] = 1

    adj += adj.T

    #plot_matrix(adj)
    return adj


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Create Regression jobs")
    args = parser.parse_args()

    # Create the infoStruct
    infoStruct = {}
    infoStruct['useMex'] = 1.
    infoStruct['seed'] = 42;
    infoStruct['hasHidden'] = 0.;
    infoStruct['lambdaNode'] = 0.1;
    infoStruct['lambdaEdge'] = 1.
    infoStruct['inferFunc'] = 'loopy';
    infoStruct['condInferFunc'] = 'loopy';
    infoStruct['nHidStates'] = 4.;

    # Create the adj matrix
    #adj_options = ['linear_vis','linear_hid','vis_hid']
    adj_options = ['linear_vis']
    nVisNodes = 99
    nHidNodes = nVisNodes if infoStruct['hasHidden'] else 0
    adj = create_adj_matrix(nVisNodes,nHidNodes,adj_options)
    infoStruct['adj'] = adj


    # Create the options
    options = {}
    options['LS'] = 0.
    options['TolFun']=1e-2;
    options['TolX'] =1e-2;
    options['Method'] ='lbfgs';
    options['Display'] = 'on';
    options['MaxIter'] = 400.;
    options['DerivativeCheck'] = 'off';
    infoStruct['options'] = options;

    sio.savemat('infoStruct.mat', {'infoStruct': infoStruct})

#    for drug,top in product(DRUGS,toplist):
#        #outstr = "('%s','%s',1,1,'%s')"%(structf,msaf,instancef)
#        outstr = 0
#
#        print >>sys.stdout, outstr
