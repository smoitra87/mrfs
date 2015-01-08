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

    if '2-vis' in options:
        for idx in range(nVisNodes-2):
            adj[idx,idx+2] = 1

    if '3-vis' in options:
        for idx in range(nVisNodes-3):
            adj[idx,idx+3] = 1

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

lambda_list = [ 1.]
nHidStates_list = [2.]
#data_list = ['PF00240', 'PF00595', 'sim3']
data_list = ['PF00240', 'PF00595']
arch_list = ['linvis', '3dvis', 'l1vis', '12vis', '123vis', 'linhid', \
             '3dhid', 'l1hid', 'linvis-linhid', 'linvis-3dhid','l1vis-l1hid' ]

train_dict = {'PF00240' :'PF00240_train.msa',
             'PF00595' : 'PF00595_train.msa',
             'sim3' : 'sim3.train.msa'}
valid_dict = {'PF00240' :'PF00240_valid.msa',
             'PF00595' : 'PF00595_valid.msa'}
test_dict = {'PF00240' :'PF00240_test.msa',
             'PF00595' : 'PF00595_test.msa'}

datamux = {'train' : train_dict,
          'valid' : valid_dict,
          'test' : test_dict}

adjf_dict = {
    'PF00240' : '1UBQ_adj.npy',
    'PF00595' : '1BE9_adj.npy'
}

l1_adjf_dict = {
    'PF00240' : 'results/l1/PF00240.list',
    'PF00595' : 'results/l1/PF00595.list'
}

data_nVisNodes = {
    'PF00240' : 69,
    'PF00595' : 81
}

def load_l1(datakey) :
    listf = l1_adjf_dict[datakey]
    with open(listf) as fin:
        edges = [map(int,line.strip().split(',')) for line in fin]

    nRes = data_nVisNodes[datakey]
    adj = np.zeros((nRes, nRes))
    for edge in edges:
        node1, node2 = edge
        adj[node1, node2] = 1

    adj += adj.T
    return adj

def create_infoStruct(archtype, datakey):
    infoStruct = default_infoStruct()

    if archtype == 'linvis':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = 0.
        infoStruct['hasHidden'] = 0.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis'])
    elif archtype == '12vis':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = 0.
        infoStruct['hasHidden'] = 0.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', '2-vis'])
    elif archtype == '123vis':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = 0.
        infoStruct['hasHidden'] = 0.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', '2-vis',\
                                                     '3-vis'])
    elif archtype == '3dvis':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = 0.
        infoStruct['hasHidden'] = 0.
        adj = create_adj_matrix(nVisNodes,nHidNodes,[])
        adj_3d = np.load(adjf_dict[datakey])
        adj[:nVisNodes,:nVisNodes] = adj_3d;
    elif archtype == 'l1vis':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = 0.
        infoStruct['hasHidden'] = 0.
        adj = create_adj_matrix(nVisNodes,nHidNodes,[])
        adj_l1 = load_l1(datakey)
        adj[:nVisNodes,:nVisNodes] = adj_l1;
    elif archtype == 'linhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_hid', 'vis_hid'])
    elif archtype == '3dhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['vis_hid'])
        adj_3d = np.load(adjf_dict[datakey])
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_3d;
    elif archtype == 'l1hid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['vis_hid'])
        adj_l1 = load_l1(datakey)
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_l1;
    elif archtype == 'l1vis-l1hid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['vis_hid'])
        adj_l1 = load_l1(datakey)
        adj[:nVisNodes,:nVisNodes] = adj_l1;
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_l1;
    elif archtype == 'linvis-linhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', \
                'linear_hid', 'vis_hid'])
    elif archtype == 'linvis-3dhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', 'vis_hid'])
        adj_3d = np.load(adjf_dict[datakey])
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_3d;
    elif archtype == 'linvis-3dhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
        adj = create_adj_matrix(nVisNodes,nHidNodes,['linear_vis', 'vis_hid'])
        adj_l1 = load_l1(datakey)
        adj[nVisNodes:nVisNodes+nHidNodes,nVisNodes:nVisNodes+nHidNodes] = \
                adj_l1;
    elif archtype == '3dvis-3dhid':
        nVisNodes = data_nVisNodes[datakey]
        nHidNodes = nVisNodes
        infoStruct['hasHidden'] = 1.
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
    parser.add_argument("--infodir", type=str, help="Directory for infoStruct")
    parser.add_argument("--datadir", type=str, help="Directory for datafile")
    parser.add_argument("--rootdir", type=str, help="Directory on workhorse")
    parser.add_argument("--resultsdir", type=str, help="Directory on workhorse")
    parser.add_argument("--param", action='store_true', help="Learn params flag")
    parser.add_argument("--eval_pll", action='store_true', help="Evaluate Pll")
    parser.add_argument("--eval_ll", action='store_true', help="Evaluate Test ll")

    args = parser.parse_args()

    if not args.infodir:
        raise ValueError('infodir missing')

    if not args.rootdir:
        raise ValueError('rootdir missing')

    if not args.datadir:
        raise ValueError('datadir missing')

    if not args.resultsdir:
        raise ValueError('resultsdir missing')

    from itertools import product

    for idx, tup in enumerate(product(data_list, arch_list, lambda_list, \
                                      nHidStates_list)):
        datakey, archtype, lambdaVal, nHidStates = tup

        if args.param:
            infoStruct = create_infoStruct(archtype, datakey)
            infoStruct['archtype'] = archtype
            infoStruct['datakey'] = datakey
            infoStruct['lambdaNode'] = lambdaVal
            infoStruct['lambdaEdge'] = lambdaVal
            infoStruct['nHidStates'] = nHidStates

            infof = "{}_{}_infoStruct.mat".format(args.prefix,idx)
            sio.savemat(os.path.join(args.infodir,infof), {'infoStruct': infoStruct})

            msaf = train_dict[datakey]
            outf =  "{}_{}_params.mat".format(args.prefix, idx)

            infof = os.path.join(args.rootdir,args.infodir, infof)
            msaf = os.path.join(args.rootdir,args.datadir, msaf)
            outf = os.path.join(args.rootdir,args.resultsdir, outf)

            outstr = "('%s','%s','%s')"%(infof, msaf, outf)
            print >>sys.stdout, outstr

        if args.eval_pll or args.eval_ll:
            paramf =  "{}_{}_params.mat".format(args.prefix, idx)
            paramf = os.path.join(args.rootdir,args.resultsdir, paramf)

            for dtype in datamux.keys():
                data_dict = datamux[dtype]
                msaf = data_dict[datakey]
                msaf = os.path.join(args.rootdir,args.datadir, msaf)

                if args.eval_pll:
                    outf =  "{}_{}_{}_pll.mat".format(args.prefix, idx, dtype)
                else:
                    outf =  "{}_{}_{}_ll.mat".format(args.prefix, idx, dtype)

                outf = os.path.join(args.rootdir,args.resultsdir, outf)

                outstr = "('%s','%s','%s')"%(paramf, msaf, outf)
                print >>sys.stdout, outstr
