function [pll, impErr] = extract_mrf_representations(paramf, msaf, outf)

load(paramf, 'adj', 'edgeMap', 'edgePot', 'edgeStruct', 'infoStruct', ...
    'nodeMap', 'nodePot', 'w');

%% Load the test data file

if isfield(infoStruct, 'hasHidden')
   hasHidden = infoStruct.hasHidden;
else
   hasHidden = 0;
end

assert(hasHidden>0, 'Representations can only be extracted for hidden var models');

[names,y] = textread(msaf, '%s %s');
if(isempty(y{1}))
  y = textread(msaf, '%s');
end;
y = converttonumericmsa(y);
y = y-min(min(y));
[nInstances,nVisNodes] = size(y);


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


if isfield(infoStruct, 'condInferFunc')
   condInferFunc = parseInferFunc(infoStruct.condInferFunc);
else
   condInferFunc = @UGM_Infer_LBP;
end


assert(isfield(infoStruct, 'nHidStates'), 'Missing nHidStates');
nHidStates = infoStruct.nHidStates;

reps = zeros(nInstances, nHidNodes * (nHidStates-1));


%% Do conditional inference
for i = 1:nInstances
    clamped = y(i,:);
    [nodeBelC,edgeBelC,logZC] = UGM_Infer_Conditional(nodePot,edgePot,...
            edgeStruct,clamped,condInferFunc);
    temp = nodeBelC((nVisNodes+1):(nVisNodes + nHidNodes),...
        1:(nHidStates-1));
    reps(i,:) = temp(:);
end

save(outf, 'reps');
end


