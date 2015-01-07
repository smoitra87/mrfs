import os
import sys
import Bio.PDB
import numpy as np

"""
Analyzes contact map for a number of different proteins
"""
import matplotlib.pyplot as plt


def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    try:
        res1coord = residue_one["CB"].coord
    except KeyError:
        res1coord = residue_one["CA"].coord
    try:
        res2coord = residue_two["CB"].coord
    except KeyError:
        res2coord = residue_two["CA"].coord

    diff_vector = res1coord - res2coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two):
    """Returns a matrix of C-alpha distances between two chains"""

    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer


def get_cmap(pdbpath, prot, pdbf):

    if args.prot not in ('ubq','pdz'): return ValueError('Unknown protein')
    pdb_code = '1BE9' if args.prot=='pr' else '1UBQ'
    pdb_filename = pdbf
    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    model = structure[0]

    if args.prot == 'ubq':
        chainA = [res for (resid,res) in enumerate(model["A"]) if \
                  5<=resid<=73 and res.get_resname() != 'HOH']
    if args.prot == 'pdz':
        chainA = [res for resid,res in enumerate(model["A"]) \
                  if res.get_resname() != 'HOH' and 313 <= res.id[1] <=391 ]

    dist_matrix = calc_dist_matrix(chainA, chainA)
    cmap8 = (dist_matrix < 8.0)
    cmap12 = (dist_matrix < 12.0)
    cmap8[np.arange(cmap8.shape[0]),np.arange(cmap8.shape[0])] = False
    cmap8 = np.triu(cmap8)
    cmap12[np.arange(cmap12.shape[0]),np.arange(cmap12.shape[0])] = False
    cmap12 = np.triu(cmap12)

    return cmap8,cmap12

pdz = [(313,318), (319,333), (334,339), (340,351),(353,391)]
baker = [(1,6), (8,22), (25,30), (31,42), (43,81)]

pdz = [range(s,e+1) for (s,e) in pdz]
baker = [range(s,e+1) for (s,e) in baker]
import itertools;
pdz = itertools.chain(*pdz)
baker = itertools.chain(*baker)

pdz = [i-313 for i in  pdz]
baker = [i-1 for i in  baker]

def get_adjmap(pdbpath, prot, pdbf):
    cmap8,cmap12 = get_cmap(pdbpath, prot, pdbf)

    if prot == 'pdz':
        adjMap = np.zeros((81,81))
        for i in range(81-1):
            adjMap[i,i+1] = 1

        x,y = zip(*itertools.product(pdz,pdz))
        x2,y2 = zip(*itertools.product(baker,baker))

        adjMap[x2,y2] = cmap8[x,y]
    elif prot == 'ubq':
        adjMap = cmap8
    else :
        raise ValueError('Unknown prot', prot)
    adjMap += adjMap.T
    return adjMap


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Calculates contact map")
    parser.add_argument('pdbf', type=str, help="Path to pdb file",
                        default='../data/1A30.pdb')
    parser.add_argument('--prot', type=str, default='ubq', help="ubq/pdz")
    parser.add_argument('--adjf', type=str)
    parser.add_argument('--display', action='store_true')
    args = parser.parse_args()

    adjMap = get_adjmap(args.pdbf, args.prot, args.pdbf)

    if args.display:
        plt.ion()
        plt.imshow(adjMap,interpolation='nearest')
        raw_input('Press any key..')

    np.save(args.adjf, adjMap)

