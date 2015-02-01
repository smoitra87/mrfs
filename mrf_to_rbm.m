function [L, minL, Pmat] = mrf_to_rbm(paramf, alpha, outf)

load(paramf);
nRes = size(adj,1);

mrf_weights = zeros(nRes * 21);

edgeEnds = edgeStruct.edgeEnds;
nEdges = size(edgeEnds,1);

for edgeIdx = 1:nEdges
    eMap = edgeMap(:,:,edgeIdx);
    eWeights = reshape(w(eMap(:)),[21,21]);
    eWeights = eWeights/2;
    e = edgeEnds(edgeIdx, :);
    e1 = e(1); e2 = e(2);
    
    mrf_weights(((e1-1)*21+1):(e1*21), ((e2-1)*21+1):(e2*21)) = eWeights;
end

mrf_weights = mrf_weights + mrf_weights';

for nodeIdx = 1:nRes
    nMap = nodeMap(nodeIdx,:);
    nWeights = w(nMap);
    
    mrf_weights(((nodeIdx-1)*21+1):(nodeIdx*21), ...
        ((nodeIdx-1)*21+1):(nodeIdx*21)) = diag(nWeights);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the MRF weight matrix to a RBM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Permute the cols to fix mrf to rbm aa map
load('aa_mrf_to_rbm.mat');
P = double(P);
nFeats = nRes * 21;
Pmap = zeros(nFeats,1);

for idx = 1:nRes
    Pmap(((idx-1)*21+1):(idx*21)) = 21*(idx-1) + P;
end

% for i = 1:length(Pmap)
%     idx = floor((i-1)/21);
%     rem = mod(i-1,21)+1;
%     Pmap(i) = idx*21 + P(rem) ;
% end


mrf_weights = mrf_weights(Pmap, Pmap);

%% Add the diagonal element to make is pos def

% fix the pos def
wrbm = mrf_weights;
min_eig = min(eig(wrbm));
if min_eig < 0
    wrbm = wrbm + eye(nRes*21) * (1 + alpha) ...
        * (-min_eig);
end

% Calculate min fill
Pcol = amd(wrbm);

% calculate row permuter
Pmat = eye(nFeats); 
Pmat = Pmat(Pcol,:); 

wrbm_minfill = Pmat*wrbm*Pmat';

% Perform chol decomps
L = chol(wrbm);
minL = chol(wrbm_minfill);


%% A bunch of assertions
assert(sum(sum((L'*L - wrbm).^2))<1e-10);
assert(sum(sum((minL'*minL - wrbm(Pcol,Pcol)).^2))<1e-10);

assert(sum(sum((minL'*minL - wrbm_minfill).^2))<1e-10);
assert(sum(sum(((Pmat'*wrbm_minfill*Pmat - wrbm).^2)))<1e-10);

temp = L'*L + eye(nRes*21)*(1 + alpha)*min_eig;
assert(sum(sum(((temp - mrf_weights).^2)))<1e-10);

save(outf, 'mrf_weights', 'L','minL','Pcol', 'Pmat', 'Pmap', 'alpha', ...
    'min_eig');

end

