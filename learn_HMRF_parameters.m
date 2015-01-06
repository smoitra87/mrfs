function [] = learn_HMRF_parameters(infof, msaf, outf)
addpath(genpath(pwd))

load(infof, 'infoStruct');

%% Global params

if isfield(infoStruct, 'useMex')
   edgeStruct.useMex = infoStruct.useMex;
else
   edgeStruct.useMex = 1;
end

if isfield(infoStruct, 'options')
   options = infoStruct.options;
else
    options.LS=0;
    options.TolFun=1e-2;
    options.TolX=1e-2;
    options.Method='lbfgs';
    options.Display='on';
    options.MaxIter=400;
    options.DerivativeCheck='off';
end

if isfield(infoStruct, 'seed')
   seed = infoStruct.seed;
else
   seed = 42;
end

s = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(s);

%% Load the MSA
[names,y] = textread(msaf, '%s %s');
if(isempty(y{1}))
  y = textread(msaf, '%s');
end;
y = converttonumericmsa(y);
y = y-min(min(y));
[nInstances,nVisNodes] = size(y);

%% Process infoStruct options

% Deal with hidden variables

if isfield(infoStruct, 'hasHidden')
   hasHidden = infoStruct.hasHidden;
else
   hasHidden = 0;
end

if hasHidden
    nHidNodes = nVisNodes;
else
    nHidNodes = 0;
end

nNodes = nVisNodes + nHidNodes;
y = y+1;
y = [y zeros(nInstances, nHidNodes)];
y = int32(y);

% Get reg params
if isfield(infoStruct, 'lambdaNode')
   lambdaNode = infoStruct.lambdaNode;
else
   lambdaNode = 1;
end

if isfield(infoStruct, 'lambdaEdge')
   lambdaEdge = infoStruct.lambdaEdge;
else
   lambdaEdge = 1;
end

if isfield(infoStruct, 'inferFunc')
   inferFunc = parseInferFunc(infoStruct.inferFunc);
else
   inferFunc = @UGM_Infer_LBP;
end


if isfield(infoStruct, 'condInferFunc')
   condInferFunc = parseInferFunc(infoStruct.condInferFunc);
else
   condInferFunc = @UGM_Infer_LBP;
end

%% Make edgeStruct

if isfield(infoStruct, 'nHidStates')
    nHidStates = infoStruct.nHidStates;
else
    nHidStates = 2;
end
nStates = [ones(1, nVisNodes)*21 ones(1, nHidNodes)*nHidStates ];

if isfield(infoStruct, 'adj')
    adj = infoStruct.adj;
else
    adj = zeros(nNodes);
end
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

%% Make Node and Edge Maps

[nodeMap,edgeMap] = UGM_make_myMRFmaps(edgeStruct);

% Initialize weights
nParams = max([nodeMap(:);edgeMap(:)]);
nNodeParams = max(nodeMap(:));
nEdgeParams = max(edgeMap(:)) - nNodeParams;
w = randn(nParams,1);

%% Learn the params

% reg params
lambda = [lambdaNode*ones(nNodeParams,1) ; lambdaEdge*ones(nEdgeParams,1)];

% Evaluate NLL
if hasHidden
    nll = UGM_MRF_Hidden_NLL(w,y,nodeMap,edgeMap,edgeStruct,...
        condInferFunc, inferFunc);
    reglaFunObj = @(w)penalizedL2(w,@UGM_MRF_Hidden_NLL,lambda,y,nodeMap,...
        edgeMap,edgeStruct,condInferFunc, inferFunc);
else
    % Compute sufficient statistics
    suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);
    nll = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct, ...
        inferFunc);
    reglaFunObj = @(w)penalizedL2(w,@UGM_MRF_NLL,lambda,nInstances,...
        suffStat,nodeMap,edgeMap,edgeStruct,inferFunc);
end

% Optimize
fprintf('Init negLogLikelihood : %f\n', nll);
disp('Training...')
tic;
w = minFunc(reglaFunObj,w, options);
toc;


%% Save the params
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
save(outf,'adj','w','edgeStruct','nodeMap', 'edgeMap', 'nodePot', ...
    'edgePot', 'infoStruct');

end

