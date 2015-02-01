import os, sys
import numpy as np
import scipy.io as sio

import matplotlib.pyplot as plt
plt.ion()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--mrf_paramf", type=str, help="MRF param file" )
    parser.add_argument("--outf", type=str, help="Output file" )
    parser.add_argument("--posdef_eps", type=float, default=0.01,\
                        help="epsilon for making matrix positive definite" )


    args = parser.parse_args()


    struct = sio.loadmat(args.mrf_paramf, squeeze_me=True, \
                         struct_as_record=False)
    w, nodeMap, edgeMap = struct['w'], struct['nodeMap'], struct['edgeMap']
    nRes = struct['adj'].shape[0]

    rbm_weights = np.zeros((nRes * 21, nRes * 21))

    edgeEnds = struct['edgeStruct'].edgeEnds
    nEdges = edgeEnds.shape[0]

    for edgeIdx in range(nEdges):
        eMap = np.squeeze(edgeMap[..., edgeIdx])-1
        eWeights = w[eMap.flatten()].reshape([21,21])
        eWeights /= 2. # Divide by 2 to account for symmetry in RBM

        e1,e2 = edgeEnds[edgeIdx, :]
        e1, e2 = e1 - 1, e2 - 1

        rbm_weights[e1*21:(e1*21+21), e2*21:(e2*21+21)] = eWeights

    rbm_weights += rbm_weights.T

    for nodeIdx in range(nRes):
        nMap = np.squeeze(nodeMap[nodeIdx,:]) - 1
        nWeights = w[nMap]

        rbm_weights[range(nodeIdx*21,nodeIdx*21+21), \
                    range(nodeIdx*21,nodeIdx*21+21)] = nWeights


    struct['wrbm'] = rbm_weights

    min_eig = np.linalg.eigvals(rbm_weights).min()

    # Make the matrix positive definite
    if min_eig < 0:
        rbm_weights += np.eye(nRes*21) * (1 + args.posdef_eps) \
                * (-min_eig)

    sio.savemat(args.outf, struct)











