function [] = learn_parameters(fname,msaf,lambdaNode,lambdaEdge, instanceweightf)
if(~isdeployed())
  addpath(genpath(pwd))
end;
                    c = '3';
                    trainType = 'pseudo';
                    switch str2double(c)
                        case 1
                            type = 'L1-L1';
                        case 2
                            type = 'L1-L2';
                        case 3
                            type = 'L1-LInf';
                    end
                    

display = 0;
if(ischar(lambdaNode))
  lambdaNode = str2double(lambdaNode);
  lambdaEdge = str2double(lambdaEdge);
end;
if(~exist('lambdaNode', 'var'))
  lambdaNode = 1;
  lambdaEdge = 1;
end;
if(exist('instanceweightf', 'var'))
  disp([' instance weightf found ' instanceweightf ]);
  instance_weights = load(instanceweightf);
  use_instance_weights = 1;
else
  use_instance_weights = 0;
end;
if(use_instance_weights==1)
  new_out_f = [fname(1:end-4) '_weightedfix_' pretty(num2str(lambdaNode)) '_' pretty(num2str(lambdaEdge)) '.mat'];
else
  new_out_f = [fname(1:end-4) '_fix_' pretty(num2str(lambdaNode)) '_' pretty(num2str(lambdaEdge)) '.mat'];
end;
%load(fname,'adjFinal','Xedge', 'edgeWeights', 'Xnode', 'nodeWeights','edgeStruct','infoStruct');
load(fname,'adjFinal','infoStruct');
%mat = load('-ASCII', msaf);

[names,X] = textread(msaf, '%s %s');
if(isempty(X{1}))
  X = textread(msaf, '%s');
end;
mat = converttonumericmsa(X);
mat = mat-min(min(mat));
nTrain = size(mat,1);
trainNdx = 1:nTrain;
numericType = str2double(c);
nFeatures = 0; % number of features for each node
nNodes = size(mat,2); % number of nodes

ising = 0; % set to 1 to use ising potentials
testType = 'loopy'; % set to 'loopy' to test with loopy belief propagation, 'exact' for exact inference
useMex = 1; % use mex files in UGM to speed things up

%lambdaNode = 1;
%lambdaEdge = 1;

%% Generate data
nInstances = nTrain;
tied = 0;

%%%%%% Change this to 21 whenever you consistently change the TP files
nStates = infoStruct.nStates; % number of states that each node can take

y = mat + 1;
X = zeros(0,size(mat,2));
adjTrue = ones(size(mat));
adjInit = adjFinal;
edgePenaltyType = 'L2';
subDisplay = display;
%example_UGMlearnSub
example_UGMlearnSub(adjInit, X, y, nNodes, nStates, nInstances, lambdaNode, lambdaEdge, trainNdx, ising, tied, testType, useMex, new_out_f, use_instance_weights, trainType, numericType, edgePenaltyType, display, subDisplay)
