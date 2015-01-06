function [NLL,g] = UGM_MRF_Hidden_NLL(w,Y,nodeMap, ...
    edgeMap,edgeStruct,condInferFunc, inferFunc,varargin)

[nNodes,maxState] = size(nodeMap);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;
nInstances = size(Y,1);

%% Calculate Model Probs

% Make potentials
[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

% Compute marginals and logZ
[nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});


NLL = 0;
g = zeros(size(w));

%% Calculate Data dependant stats

for i = 1:nInstances
    [nodeBelC,edgeBelC,logZC] = UGM_Infer_Conditional(nodePot,edgePot,edgeStruct,Y(i,:),condInferFunc,varargin{:});

    
    
    [logZC_update, g_update] = update_grad(g, Y,i, nodePot, edgePot, nodeMap, edgeMap, ...
        nodeBelC, edgeBelC, nodeBel, edgeBel, edgeEnds, nStates);

 %   temp.g = g;
%    [logZC_updateC] = UGM_MRF_HiddenC(g, Y, int32(i), nodePot, edgePot, nodeMap, edgeMap, ...
%        nodeBelC, edgeBelC, nodeBel, edgeBel, edgeEnds, nStates);

    % Update NLL
    NLL = NLL - logZC - logZC_update + logZ;    
    g  = g + g_update;
   
end
end


function [logZC_update, g_update] = update_grad(g, Y, i, nodePot, edgePot, nodeMap, edgeMap, ...
    nodeBelC, edgeBelC, nodeBel, edgeBel, edgeEnds, nStates)

y = Y(i,:);
g_update = zeros(size(g));
logZC_update = 0;
nNodes = size(nodePot,1);
nEdges = size(edgeEnds,1);

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


if nargout > 1
    for n = 1:nNodes
        for s = 1:nStates(n)
            if nodeMap(n,s) > 0
                if y(n) == 0
                    obs = nodeBelC(n,s);
                elseif s == y(n)
                    obs = 1;
                else
                    obs = 0;
                end
                g_update(nodeMap(n,s)) = g_update(nodeMap(n,s)) + nodeBel(n,s) - obs;
            end
        end
    end
    for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
        for s1 = 1:nStates(n1)
            for s2 = 1:nStates(n2)
                if edgeMap(s1,s2,e) > 0
                    if y(n1) == 0 && y(n2) == 0
                        obs = edgeBelC(s1,s2,e);
                    elseif y(n1) == 0 && y(n2) == s2
                        obs = edgeBelC(s1,y(n2),e);
                    elseif y(n1) == s1 && y(n2) == 0
                        obs = edgeBelC(y(n1),s2,e);
                    elseif s1 == y(n1) && s2 == y(n2)
                        obs = 1;
                    else
                        obs = 0;
                    end
                    g_update(edgeMap(s1,s2,e)) = g_update(edgeMap(s1,s2,e)) + edgeBel(s1,s2,e) - obs;
                end
            end
        end
    end
end
end

