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

def default_infoStruct():
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
    return infoStruct

lambda_list = [ 1e-2, 1e-1, 1.]
nHidStates_list = [2., 5., 10.]
#data_list = ['PF00240', 'PF00595', 'sim3']
data_list = ['PF00240']
arch_list = ['linvis-linhid', '3dvis-3dhid', 'linvis-3dhid']

msaf_dict = {'PF00240' :'PF00240_train.msa',
             'PF00595' : 'PF00595_train.msa',
             'sim3' : 'sim3.train.msa'}
adjf_dict = {
    'PF00240' : '1UBQ_adj.npy',
    'PF00595' : '1BE9_adj.npy'
}


data_nVisNodes = {
    'PF00240' : 69,
    'PF00595' : 81
}

def create_infoStruct(archtype, datakey):
    infoStruct = default_infoStruct()
    if archtype == 'linvis-linhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', \
                'linear_hid', 'vis_hid'])
    elif archtype == 'linvis-3dhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', 'vis_hid'])
        adj_3d = np.load(adjf_dict[datakey])
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_3d;
    elif archtype == '3dvis-3dhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        adj = create_adj_matrix(nVisNodes,nHidNodes,['vis_hid'])
        adj_3d = np.load(adjf_dict[datakey])
        adj[:nVisNodes,:nVisNodes] = adj_3d;
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_3d;
    else:
        raise ValueError('Unknown arch', archtype)

    infoStruct['adj'] = adj
    return infoStruct

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Create param jobs")
    parser.add_argument("--prefix", type=str, default="hmrf")
    parser.add_argument("--infofdir", type=str, help="Directory for infoStruct")
    args = parser.parse_args()

    if not args.infofdir:
        raise ValueError('infofdir missing')

    from itertools import product

    for idx, tup in enumerate(product(data_list, arch_list, lambda_list, \
                                      nHidStates_list)):
        datakey, archtype, lambdaVal, nHidStates = tup
        infoStruct = create_infoStruct(archtype, datakey)
        infoStruct['archtype'] = archtype
        infoStruct['datakey'] = datakey
        infoStruct['lambdaNode'] = lambdaVal
        infoStruct['lambdaEdge'] = lambdaVal
        infoStruct['nHidStates'] = nHidStates

        infof = "{}_{}_infoStruct.mat".format(args.prefix,idx)
        sio.savemat(os.path.join(args.infofdir,infof), {'infoStruct': infoStruct})

        msaf = msaf_dict[datakey]
        outf =  "{}_{}_params.mat".format(args.prefix, idx)

        outstr = "('%s','%s','%s')"%(infof, msaf, outf)
        print >>sys.stdout, outstr
