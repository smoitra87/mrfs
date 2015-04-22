import os,sys
import glob
import scipy.io as sio
import numpy as np



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Creates CSV from mat files")
    parser.add_argument("--globstr",type=str, help='globstring')
    parser.add_argument("--outf",type=str, help='outf')
    args = parser.parse_args()

    files = glob.glob(args.globstr)

    pllfs = [f for f in files if '_pll' in f ]
    blosumfs = [f for f in files if '_blosum' in f ]
    llfs = [f for f in files if '_ll' in f ]
    paramfs = [f for f in files if '_params' in f ]

    mat = sio.loadmat(paramfs[0], squeeze_me=True)
    infoStruct = mat['infoStruct']
    headers = infoStruct.dtype.names

    remove_headers = ['hasHidden', 'seed', 'options', 'adj', 'useMex']
    headers = sorted([h for h in headers if h not in remove_headers])
    headers = ['jobkey', 'train-ll','valid-ll', 'test-ll'] + \
            ['train-pll', 'valid-pll', 'test-pll'] + \
            ['train-imperr', 'valid-imperr', 'test-imperr'] + \
            ['train-bl90', 'valid-bl90', 'test-bl90'] + \
            ['train-imperr-serr', 'valid-imperr-serr', 'test-imperr-serr'] +\
            ['train-bl90-serr', 'valid-bl90-serr', 'test-bl90-serr'] +\
            headers;

    records = []
    for paramf in paramfs:
        base_key = os.path.splitext(os.path.basename(paramf))[0]
        base_key = "_".join(base_key.split("_")[:-1])

        llf = [f for f in llfs if base_key in f]
        pllf = [f for f in pllfs if base_key in f]
        blosumf = [f for f in blosumfs if base_key in f]

        mat = sio.loadmat(paramf, squeeze_me=True)
        infoStruct = mat['infoStruct']

        metric = {}
        metric['jobkey'] = base_key

        # initialize metrics
        for t in ('train', 'valid', 'test'):
            metric['{}-ll'.format(t)] = None
            metric['{}-pll'.format(t)] = None
            metric['{}-imperr'.format(t)] = None
            metric['{}-bl90'.format(t)] = None
            metric['{}-imperr-serr'.format(t)] = None
            metric['{}-bl90-serr'.format(t)] = None


        for f in llf :
            for t in ('train', 'valid', 'test'):
                if t in f:
                    metric['{}-ll'.format(t)] = sio.loadmat(f, squeeze_me=True)['avgLL']

        for f in pllf :
            for t in ('train', 'valid', 'test'):
                if t in f:
                    metric['{}-pll'.format(t)] = sio.loadmat(f, squeeze_me=True)['pll']
                    metric['{}-imperr'.format(t)] = sio.loadmat(f, squeeze_me=True)['impErr']
                    imperr = sio.loadmat(f, squeeze_me=True)['imperr_raw']
                    metric['{}-imperr-serr'.format(t)] = np.std(imperr) / np.sqrt(len(imperr))

        for f in blosumf :
            for t in ('train', 'valid', 'test'):
                if t in f:
                    metric['{}-bl90'.format(t)] = sio.loadmat(f, squeeze_me=True)['impErr']
                    bl90 = sio.loadmat(f, squeeze_me=True)['imperr_raw']
                    metric['{}-bl90-serr'.format(t)] = np.std(bl90) / np.sqrt(len(imperr))

        record = [metric[h] if h in metric else str(infoStruct[h]) for h in headers]
        records.append(record)

    import csv

    if not args.outf :
        raise ValueError("outf missing")

    with open(args.outf, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for record in records:
            writer.writerow(record)


