import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import scipy.io as sio

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--countf", help="Count file")
    parser.add_argument("--reps", help="Mat containg reps")
    parser.add_argument("--title", help="Title of plot")
    parser.add_argument("--xlabel")
    parser.add_argument("--ylabel")
    parser.add_argument("--bar_width", type=float, default=0.2, help="Thickness of bars")
    parser.add_argument("--opacity", type=float, default=1., help="Opacity")
    parser.add_argument("--savef", type=str, help="Name of file to save fig")
    parser.add_argument("--show", action='store_true', help="Show fig")
    args = parser.parse_args()

    x = sio.loadmat(args.reps)
    reps  = x['reps']
    pca = PCA(n_components=2)
    x_r = pca.fit(reps).transform(reps)

    plt.figure()
    with open(args.countf) as fin :
        y = [int(l.strip()) for l in fin]
    plt.scatter(x_r[:,0],x_r[:,1],c=y)

    if args.xlabel:
        plt.xlabel(args.xlabel)
    if args.ylabel:
        plt.ylabel(args.ylabel)
    if args.title:
        plt.title(args.title)

    if args.savef:
        plt.savefig(args.savef)

    if args.show:
        plt.show()
