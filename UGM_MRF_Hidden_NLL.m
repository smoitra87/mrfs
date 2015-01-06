function [NLL,g] = UGM_MRF_Hidden_NLL(w,nodeMap, ...
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
    for n = 1:nNodes
        if Y(i,n)~=0
            logZC = logZC + log(nodePot(n,Y(i,n)));
        end
    end
    for e = 1:nEdges
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
        if Y(i,n1)~=0 && Y(i,n2)~=0
            logZC = logZC + log(edgePot(Y(i,n1),Y(i,n2),e));
        end
    end
    
    % Update NLL
    NLL = NLL - logZC + logZ;
    
    if nargout > 1
        for n = 1:nNodes
            for s = 1:nStates(n)
                if nodeMap(n,s) > 0
                    if Y(i,n) == 0
                        obs = nodeBelC(n,s);
                    elseif s == Y(i,n)
                        obs = 1;
                    else
                        obs = 0;
                    end
                    g(nodeMap(n,s)) = g(nodeMap(n,s)) + nodeBel(n,s) - obs;
                end
            end
        end
        for e = 1:nEdges
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);
            for s1 = 1:nStates(n1)
                for s2 = 1:nStates(n2)
                    if edgeMap(s1,s2,e) > 0
                        if Y(i,n1) == 0 && Y(i,n2) == 0
                            obs = edgeBelC(s1,s2,e);
                        elseif Y(i,n1) == 0 && Y(i,n2) == s2
                            obs = edgeBelC(s1,Y(i,n2),e);
                        elseif Y(i,n1) == s1 && Y(i,n2) == 0
                            obs = edgeBelC(Y(i,n1),s2,e);
                        elseif s1 == Y(i,n1) && s2 == Y(i,n2)
                            obs = 1;
                        else
                            obs = 0;
                        end
                        g(edgeMap(s1,s2,e)) = g(edgeMap(s1,s2,e)) + edgeBel(s1,s2,e) - obs;
                    end
                end
            end
        end
    end
end


