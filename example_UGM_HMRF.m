%clear all
close all
addpath(genpath('.'))
seed = 42;
s = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(s);
msaf = 'clean_prXedit_1k_train.msa';


%% Load the MSA
[names,y] = textread(msaf, '%s %s');
if(isempty(y{1}))
  y = textread(msaf, '%s');
end;
y = converttonumericmsa(y);
y = y-min(min(y));
[nInstances,nVisNodes] = size(y);
nHidNodes = nVisNodes;
nNodes = nVisNodes + nHidNodes;
y = y+1;
y = [y zeros(nInstances, nHidNodes)];
y = int32(y);


%% Make edgeStruct
nStates = [ones(1, nVisNodes)*21 ones(1, nHidNodes)*2 ];

adj = zeros(nVisNodes + nHidNodes);
for i = 1:nVisNodes
    adj(i,i+nVisNodes) = 1;
end

for i = 1:nVisNodes-1
    adj(i, i+1) = 1;
end

for i = 1:nHidNodes-1
    adj(nVisNodes + i,nVisNodes+i+1) = 1;
end
adj = adj+adj';

edgeStruct = UGM_makeEdgeStruct(adj,nStates);
edgeStruct.useMex  = 1;

maxState = max(nStates);
nEdges = edgeStruct.nEdges;

%% Make Node and Edge Pots
[nodeMap,edgeMap] = UGM_make_hiddenMRFmaps(edgeStruct);

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = randn(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

%% Learn the params

% Evaluate NLL
nll = UGM_MRF_Hidden_NLL(w,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP, @UGM_Infer_LBP)

% Set up regularization parameters
lambdaVal = 1;
lambda = lambdaVal*ones(size(w));
reglaFunObj = @(w)penalizedL2(w,@UGM_MRF_Hidden_NLL,lambda,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP, @UGM_Infer_LBP);

% LBFGS to find the weights
display('Training...');
options.LS=0;
options.TolFun=1e-2;
options.TolX=1e-2;
options.Method='lbfgs';
options.Display='on';
options.MaxIter=400;
options.DerivativeCheck='off';

% Optimize
tic;
w = minFunc(reglaFunObj,w, options);
toc;

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);


%% Do decoding and inference

decode = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);

[nodeBel,edgeBel,logZ] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);

samples = UGM_Sample_Gibbs(nodePot,edgePot,edgeStruct,1000);
figure;
imagesc(samples')
title('Samples from MRF model');
fprintf('(paused)\n');


%% Do conditional decoding/inference/sampling in learned model

clamped = zeros(nNodes,1);
clamped(1:2:90) = 2;

condDecode = UGM_Decode_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Decode_LBP);
condNodeBel = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Infer_LBP);
condSamples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Sample_Gibbs);

figure;
imagesc(condSamples')
title('Conditional samples from MRF model');
fprintf('(paused)\n');







