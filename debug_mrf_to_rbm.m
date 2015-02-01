function [avgLL] = debug_mrf_to_rbm(paramf, msaf, rbmparamf)

load(paramf, 'adj', 'edgeMap', 'edgePot', 'edgeStruct', 'infoStruct', ...
    'nodeMap', 'nodePot', 'w');

addpath(genpath('.'));

load(rbmparamf);

%% Load the test data file

[names,y] = textread(msaf, '%s %s');
if(isempty(y{1}))
  y = textread(msaf, '%s');
end;
y = cellstr(y{1});
y = converttonumericmsa(y);
y = y-min(min(y));

% Just take a single sequece
[nInstances,nVisNodes] = size(y);

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

%% Load inference related data

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


%% Calculate Log likelihood


suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);
nll = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,inferFunc);
%nll = UGM_MRF_Hidden_NLL(w,y,nodeMap,edgeMap,edgeStruct,...
%    condInferFunc, inferFunc);
ene = w'*suffStat;

fprintf('Energy : %f\n', ene);

%% Expanded sequence

load('aa_mrf_to_rbm.mat');
nFeats = size(L,1);

yexp = zeros(nFeats,1);
for i = 1:length(y)
    yexp((i-1)*21+y(i)) = 1;
end
yexp = yexp(Pmap);

%mrf_weights = mrf_weights(Pmap, Pmap);
%assert(sum(sum((yexp(Pmap) - Pmapmat*yexp).^2))<1e-10);
%mrf_weights = Pmapmat*mrf_weights*Pmapmat';

ene_exp =  yexp'*mrf_weights*yexp;
assert( sum(ene - ene_exp)^2 < 1e-10);




end



