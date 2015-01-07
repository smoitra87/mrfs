import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
import csv

def plot_coldict(d, dat, metkey):
    max_labels = max(len(k.keys()) for k in d) -1
    xticklabels = [d2['name']+":"+"|".join(k for k in sorted(d2.keys()) if k!='name') for d2 in d]
    labels = [[k for k in sorted(d2.keys()) if k!='name'] for d2 in d]
    means = [[d2[l]['mean'] for l in labs] for (labs,d2) in zip(labels,d)]
    stds = [[d2[l]['std'] for l in labs] for (labs,d2) in zip(labels,d)]

    means = [l+[0.]*(max_labels-len(l)) for l in means]
    stds = [l+[0.]*(max_labels-len(l)) for l in stds]
    means, stds = zip(*means), zip(*stds)

    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)

    N = len(labels)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars

    colors = 'rgbymc'
    rects = []

    for idx,means2 in enumerate(means):
        rect = ax.bar(ind+width*idx, means2, width, color = colors[idx], yerr=stds[idx],\
                error_kw={'ecolor':'gray', 'lw':2})
        rects.append(rect)

    metkey = 'neg-'+metkey if 'll' in metkey else metkey

    ax.set_ylabel(metkey)
    ax.set_title('HyperParam Scan: Dataset {}'.format(dat))
    ax.set_xticks(ind+width)
    ax.set_xticklabels( xticklabels )

    zed = [tick.label.set_fontsize(7) for tick in ax.xaxis.get_major_ticks()]

    ax.legend(rects, ["Lab{}".format(i+1) for i in range(max_labels)])

    plt.tight_layout()
    plt.savefig("{}_{}.png".format(dat,metkey), dpi=100)
    plt.close()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = 'analyzes collected results')
    parser.add_argument('--csvf', type=str, help='csv file')
    parser.add_argument('--best', type=int, default=0, help='Best result to consider')
    parser.add_argument('--metric', type=str, default='valid-ll', \
                        help='Metric')
    args = parser.parse_args()

    if not args.csvf:
        raise ValueError("Missing csvf")

    with open(args.csvf, 'rb') as fin:
        reader = csv.reader(fin)
        header = next(reader)
        print header
        rows = [row for row in reader]

    metric_idx = next(idx for (idx,e) in enumerate(header) if e == args.metric)
    dataset_idx = next(idx for (idx,e) in enumerate(header) if e == 'datakey')
    start_idx = next(idx for (idx,e) in enumerate(header) if e == 'archtype')
    end_idx = next(idx for (idx,e) in enumerate(header) if e == 'nHidStates')

    # Convert ll into negLL
    extract_col = lambda(colidx): [row[colidx] for row in rows]

    # Plot fitted Gaussians
    metric = extract_col(metric_idx)
    metric = np.asarray(map(float, metric))

    # Convert ll into negLL
    if 'll' in args.metric and 'nll' not in args.metric:
        metric = -metric;
    from operator import itemgetter
    sortidx = np.argsort(metric)
    np.sort(metric)

    rows = [rows[idx] for idx in sortidx]

    dataset = extract_col(dataset_idx)
    uniq_d = set(dataset)
    dat_dict  = {}
    for dat in uniq_d:
        datidx = [jdx for (jdx,e2) in enumerate(dataset) if e2 == dat]

        if args.best > 0:
            datidx = datidx[:args.best]
        dat_dict[dat] = datidx

    from collections import defaultdict
    coldict_dict = defaultdict(list);
    for dat in sorted(dat_dict.keys()):
        print
        print("Analyzing Dataset {} .....").format(dat)
        print

        for idx in range(start_idx, end_idx+1):
            col = extract_col(idx)
            uniq_e = set(col)

            coldict = {}
            coldict['name'] = header[idx]

            for e in uniq_e :
                colidx = [jdx for (jdx,e2) in enumerate(col) if e2 == e and jdx in dat_dict[dat]]
                coldict[e] = {}
                coldict[e]['mean'] = np.mean(metric[colidx])
                coldict[e]['std'] = np.std(metric[colidx])
            print coldict
            coldict_dict[dat].append(coldict)

    coldict_dict = dict(coldict_dict)
    for dat in coldict_dict:
        plot_coldict(coldict_dict[dat], dat, args.metric)









