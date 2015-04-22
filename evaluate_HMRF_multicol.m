function [pll, impErr] = evaluate_HMRF_Pll(paramf, msaf, outf, multicol_file, condInferFunc, varargin)

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

% Get the conditional inference function
condInferFunc = parseInferFunc(condInferFunc);

% Load the set of columns to run multicolumn imputation
disp(exist(multicol_file))
assert(exist(multicol_file, 'file')>0, 'Multi col file not found')
load(multicol_file, 'multicols')
[nBlocks, nCols] = size(multicols)

useBlosum = 0;
for i = 1:length(varargin)
    if strcmp(varargin{i},'blosum90')
        useBlosum = 1;
        break
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
    for multidx = 1:nBlocks
        clamped = y(i,:); clamped(multicols(multidx,:)) = 0; 
        [nodeBelC,edgeBelC,logZC] = UGM_Infer_Conditional(nodePot,edgePot,...
            edgeStruct,clamped,condInferFunc);
        for j = 1:nCols
            [maxMarg,maxMargIdx] = max(nodeBelC(multicols(multidx,j),:));
            if useBlosum > 0
                imperr_raw(i) = imperr_raw(i) + blosum90(maxMargIdx, y(i,multicols(multidx, j)));
            else
                imperr_raw(i) = imperr_raw(i) + (maxMargIdx ~= y(i,multicols(multidx, j)));
            end
            pll(i) = pll(i) + log(nodeBelC(multicols(multidx, j),y(i,multicols(multidx,j))));
        end
    end
    imperr_raw(i) = imperr_raw(i) / (nBlocks * nCols);
end
impErr = mean(imperr_raw);
pll_raw = pll;
pll = mean(pll);

fprintf('Average PLL : %f\n', pll);
fprintf('Imputation Error : %f\n', impErr);
save(outf, 'pll', 'impErr', 'pll_raw', 'imperr_raw');
end

