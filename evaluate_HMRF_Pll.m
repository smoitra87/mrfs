function [pll, impErr] = evaluate_HMRF_Pll(paramf, msaf, outf)

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

if isfield(infoStruct, 'usePseudo')
    pseudo = infoStruct.usePseudo;
else
    pseudo = 0;
end


%% Do conditional inference
err = 0;
pll = zeros(nInstances, 1);
for i = 1:nInstances
    
    if pseudo > 0
        nodeBel = UGM_MRF_InferPseudo(y(i,:), nodePot,edgePot, edgeStruct);
        for j = 1:nVisNodes
            [maxMarg,maxMargIdx] = max(nodeBel(j,:));
            err = err + (maxMargIdx ~= y(i,j));
            pll(i) = pll(i) + log(nodeBel(j,y(i,j)));
        end
    else
        for j = 1:nVisNodes
            clamped = y(i,:); clamped(j) = 0;
            [nodeBelC,edgeBelC,logZC] = UGM_Infer_Conditional(nodePot,edgePot,...
                edgeStruct,clamped,condInferFunc);
            [maxMarg,maxMargIdx] = max(nodeBelC(j,:));
            err = err + (maxMargIdx ~= y(i,j));
            pll(i) = pll(i) + log(nodeBelC(j,y(i,j)));
        end
    end
end
impErr = err/(nInstances*nVisNodes);
pll = mean(pll);


fprintf('Average PLL : %f\n', pll);
fprintf('Imputation Error : %f\n', impErr);
save(outf, 'pll', 'impErr');
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



