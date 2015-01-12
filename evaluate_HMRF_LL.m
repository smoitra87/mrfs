function [avgLL] = evaluate_HMRF_LL(paramf, msaf, outf)

load(paramf, 'adj', 'edgeMap', 'edgePot', 'edgeStruct', 'infoStruct', ...
    'nodeMap', 'nodePot', 'w');

%% Load the test data file

[names,y] = textread(msaf, '%s %s');
if(isempty(y{1}))
  y = textread(msaf, '%s');
end;
y = converttonumericmsa(y);
y = y-min(min(y));
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

if hasHidden
    [nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct);
    
    %% Do conditional inference
    logZv = zeros(nInstances, 1);
    
    for i = 1:nInstances
        clamped = y(i,:);
        [nodeBelC,edgeBelC,logZC] = UGM_Infer_Conditional(nodePot,edgePot,...
            edgeStruct,clamped,condInferFunc);
        logZC = logZC + update_logZC(y, i, nodePot, edgePot, edgeStruct.edgeEnds);
        logZv(i) = logZC;
    end
    avgLL = mean(logZv - logZ);
else
    suffStat = UGM_MRF_computeSuffStat(y,nodeMap,edgeMap,edgeStruct);
    nll = UGM_MRF_NLL(w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,inferFunc);
    %nll = UGM_MRF_Hidden_NLL(w,y,nodeMap,edgeMap,edgeStruct,...
    %    condInferFunc, inferFunc);
    avgLL = -nll / nInstances;    
end

fprintf('Average LL : %f', avgLL);
save(outf, 'avgLL');

end


function [logZC_update] = update_logZC(Y, i, nodePot, edgePot, edgeEnds)

y = Y(i,:);
logZC_update = 0;
nNodes = size(nodePot,1);
nEdges = size(edgePot,3);

for n = 1:nNodes
    if y(n)~=0
        logZC_update = logZC_update + log(nodePot(n,y(n)));
    end
end
for e = 1:nEdges
    n1 = edgeEnds(e,1);
    n2 = edgeEnds(e,2);
    if y(n1)~=0 && y(n2)~=0
        logZC_update = logZC_update + log(edgePot(y(n1),y(n2),e));
    end
end

end



