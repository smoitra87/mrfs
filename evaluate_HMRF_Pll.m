function [pll, impErr] = evaluate_HMRF_Pll(paramf, msaf, outf, varargin)

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

useBlosum = 0;
if length(varargin) > 0
    for i = 1:length(varargin)
        if strcmp(varargin{i},'blosum90')
            useBlosum = 1;
            break
        end
    end
end

blosum90 = [[5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1, -6],
[-2, 6, -1, -3, -5, 1, -1, -3, 0, -4, -3, 2, -2, -4, -3, -1, -2, -4, -3, -3, -6],
[-2, -1, 7, 1, -4, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -5, -3, -4, -6],
[-3, -3, 1, 7, -5, -1, 1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5, -6],
[-1, -5, -4, -5, 9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -6],
[-1, 1, 0, -1, -4, 7, 2, -3, 1, -4, -3, 1, 0, -4, -2, -1, -1, -3, -3, -3, -6],
[-1, -1, -1, 1, -6, 2, 6, -3, -1, -4, -4, 0, -3, -5, -2, -1, -1, -5, -4, -3, -6],
[0, -3, -1, -2, -4, -3, -3, 6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -6],
[-2, 0, 0, -2, -5, 1, -1, -3, 8, -4, -4, -1, -3, -2, -3, -2, -2, -3, 1, -4, -6],
[-2, -4, -4, -5, -2, -4, -4, -5, -4, 5, 1, -4, 1, -1, -4, -3, -1, -4, -2, 3, -6],
[-2, -3, -4, -5, -2, -3, -4, -5, -4, 1, 5, -3, 2, 0, -4, -3, -2, -3, -2, 0, -6],
[-1, 2, 0, -1, -4, 1, 0, -2, -1, -4, -3, 6, -2, -4, -2, -1, -1, -5, -3, -3, -6],
[-2, -2, -3, -4, -2, 0, -3, -4, -3, 1, 2, -2, 7, -1, -3, -2, -1, -2, -2, 0, -6],
[-3, -4, -4, -5, -3, -4, -5, -5, -2, -1, 0, -4, -1, 7, -4, -3, -3, 0, 3, -2, -6],
[-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4, 8, -2, -2, -5, -4, -3, -6],
[1, -1, 0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2, 5, 1, -4, -3, -2, -6],
[0, -2, 0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2, 1, 6, -4, -2, -1, -6],
[-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2, 0, -5, -4, -4, 11, 2, -3, -6],
[-3, -3, -3, -4, -4, -3, -4, -5, 1, -2, -2, -3, -2, 3, -4, -3, -2, 2, 8, -3, -6],
[-1, -3, -4, -5, -2, -3, -3, -5, -4, 3, 0, -3, 0, -2, -3, -2, -1, -3, -3, 5, -6],
[-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1]];

%% Do conditional inference
pll = zeros(nInstances, 1);
imperr_raw = zeros(nInstances, 1);
for i = 1:nInstances
    if pseudo > 0
        nodeBel = UGM_MRF_InferPseudo(y(i,:), nodePot,edgePot, edgeStruct);
        for j = 1:nVisNodes
            [maxMarg,maxMargIdx] = max(nodeBel(j,:));
            if useBlosum > 0
                imperr_raw(i) = imperr_raw(i) + blosum90(maxMargIdx, y(i,j));
            else
                imperr_raw(i) = imperr_raw(i) + (maxMargIdx ~= y(i,j));
            end
            pll(i) = pll(i) + log(nodeBel(j,y(i,j)));
        end
    else
        for j = 1:nVisNodes
            clamped = y(i,:); clamped(j) = 0;
            [nodeBelC,edgeBelC,logZC] = UGM_Infer_Conditional(nodePot,edgePot,...
                edgeStruct,clamped,condInferFunc);
            [maxMarg,maxMargIdx] = max(nodeBelC(j,:));
            if useBlosum > 0
                imperr_raw(i) = imperr_raw(i) + blosum90(maxMargIdx, y(i,j));
            else
                imperr_raw(i) = imperr_raw(i) + (maxMargIdx ~= y(i,j));
            end
            pll(i) = pll(i) + log(nodeBelC(j,y(i,j)));
        end
    end
    imperr_raw(i) = imperr_raw(i) / nVisNodes;
end
impErr = mean(imperr_raw);
pll_raw = pll;
pll = mean(pll);


fprintf('Average PLL : %f\n', pll);
fprintf('Imputation Error : %f\n', impErr);
save(outf, 'pll', 'impErr', 'pll_raw', 'imperr_raw');
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



