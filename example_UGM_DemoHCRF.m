clear all
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
y = int32(y+1);
nStates = 21;

[nInstances,nNodes] = size(y);


%% Make edgeStruct
% load(adjf,'adjFinal','infoStruct');
% adj = adjFinal;

adj = zeros(nNodes);
for i = 1:nNodes-1
    adj(i,i+1) = 1;
end
adj = adj+adj';

edgeStruct = UGM_makeEdgeStruct(adj,nStates);
edgeStruct.useMex  = 1;
maxState = max(nStates);
nEdges = edgeStruct.nEdges;

%% Make Node and Edge Pots
% TODO - fix the node and the edge map to untie the paramms

% % Make nodeMap
% nodeMap = zeros(nNodes,maxState,'int32');
% nodeMap(:,1) = 1;
% 
% % Make edgeMap
% edgeMap = zeros(maxState,maxState,nEdges,'int32');
% 
% edgeMap(1,1,:) = 2;
% edgeMap(2,2,:) = 2;
ising= 0 ; tied = 0;
[nodeMap,edgeMap] = UGM_makeMRFmaps(nNodes,edgeStruct,ising,tied);

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

% Example of making potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

%% Learn the params

% Compute sufficient statistics
suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);

% Evaluate NLL
nll = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain)

% Set up regularization parameters
lambdaVal = 1;
lambda = lambdaVal*ones(size(w));
reglaFunObj = @(w)penalizedL2(w,@UGM_MRF_NLL,lambda,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,@UGM_Infer_Chain);

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
w = minFunc(reglaFunObj,w, options);

% Now make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
nodePot(1,:)
edgePot(:,:,1)
fprintf('(paused)\n');


%% Do decoding and inference

decode = UGM_Decode_Chain(nodePot,edgePot,edgeStruct)

[nodeBel,edgeBel,logZ] = UGM_Infer_Chain(nodePot,edgePot,edgeStruct);
nodeBel

samples = UGM_Sample_Chain(nodePot,edgePot,edgeStruct);
figure;
imagesc(samples')
title('Samples from MRF model');
fprintf('(paused)\n');
pause

%% Do conditional decoding/inference/sampling in learned model

clamped = zeros(nNodes,1);
clamped(1:2:90) = 2;

condDecode = UGM_Decode_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Decode_Tree)
condNodeBel = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Infer_Tree)
condSamples = UGM_Sample_Conditional(nodePot,edgePot,edgeStruct,clamped,@UGM_Sample_Tree);

figure;
imagesc(condSamples')
title('Conditional samples from MRF model');
fprintf('(paused)\n');







