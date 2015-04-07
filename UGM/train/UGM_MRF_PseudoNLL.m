function [NLL,g] = UGM_MRF_PseudoNLL(w,Y,nodeMap,edgeMap,edgeStruct)

[nNodes,maxState] = size(nodeMap);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;

nInstances = size(Y,1);
NLL = 0;
g = zeros(size(w));

[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);

for i = 1:nInstances
    % Make potentials
        for n = 1:nNodes
            % Find Neighbors
            edges = E(V(n):V(n+1)-1);
            
            % Compute Probability of Each State with Neighbors Fixed
            pot = nodePot(n,1:nStates(n));
            for e = edges(:)'
                n1 = edgeEnds(e,1);
                n2 = edgeEnds(e,2);
                
                if n == edgeEnds(e,1)
                    ep = edgePot(1:nStates(n),Y(i,n2),e).';
                else
                    ep = edgePot(Y(i,n1),1:nStates(n),e);
                end
                pot = pot .* ep;
            end
            
            % Update Objective
            NLL = NLL - log(pot(Y(i,n))) + log(sum(pot));
            
            %% Update Gradient
            if nargout > 1
                nodeBel = pot/sum(pot);
                
                % Update Gradient of Node Weights
                for s = 1:nStates(n)
                    
                    if nodeMap(n,s) > 0
                        if s == Y(i,n)
                            obs = 1;
                        else
                            obs = 0;
                        end
                        g(nodeMap(n,s)) = g(nodeMap(n,s)) + (nodeBel(s) - obs);
                    end
                    
                end
                
                % Update Gradient of Edge Weights
                for e = edges(:)'
                    
                    n1 = edgeEnds(e,1);
                    n2 = edgeEnds(e,2);
                    
                    for s = 1:nStates(n)
                        if n == n1
                            s1 = s;
                            neigh = n2;
                            s2 = Y(i,neigh);
                        else
                            s2 = s;
                            neigh = n1;
                            s1 = Y(i,neigh);
                        end
                        
                        if edgeMap(s1,s2,e) > 0
                            if s == Y(i,n)
                                obs = 1;
                            else
                                obs = 0;
                            end
                            g(edgeMap(s1,s2,e)) = g(edgeMap(s1,s2,e)) + (nodeBel(s) - obs);
                        end
                        
                    end
                end
            end
        end
end
