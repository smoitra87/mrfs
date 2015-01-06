clear all;
infoStruct.useMex = 1;

options.LS=0;
options.TolFun=1e-2;
options.TolX=1e-2;
options.Method='lbfgs';
options.Display='on';
options.MaxIter=400;
options.DerivativeCheck='off';
infoStruct.options = options;

infoStruct.seed = 42;
infoStruct.hasHidden = 1;
infoStruct.lambdaNode = 0.1;
infoStruct.lambdaEdge = 1;
infoStruct.inferFunc = 'loopy';
infoStruct.condInferFunc = 'loopy';

infoStruct.nHidStates = 4;

nVisNodes = 99;
nHidNodes = 99;
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
infoStruct.adj = adj;

save('infoStruct.mat','infoStruct')