import numpy as np
import scipy.io as sio


mrf = 'ARNDCQEGHILKMFPSTWYV-';
rbm = 'AVLIPFWMGSTCYNQDEKRH-';

P = []
for idx, aa in enumerate(mrf):
    P.append(next(idx2 for (idx2,aa2) in enumerate(rbm) if aa == aa2))
P = np.asarray(P) + 1;

P2 = []
for idx, aa in enumerate(rbm):
    P2.append(next(idx2 for (idx2,aa2) in enumerate(mrf) if aa == aa2))
P2 = np.asarray(P2) + 1;

sio.savemat('aa_mrf_to_rbm.mat', {'P':P2})

map_to_mrf = dict(enumerate(P2 - 1))
map_to_rbm =  dict(enumerate(P - 1))



