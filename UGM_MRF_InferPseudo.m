function [nodeBel] = UGM_MRF_InferPseudo(y,nodePot,edgePot,edgeStruct)

[nNodes,maxState] = size(nodePot);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;
NLL = 0;

nodeBel = zeros(nNodes, maxState);

for n = 1:nNodes
    % Find Neighbors
    edges = E(V(n):V(n+1)-1);
    
    % Compute Probability of Each State with Neighbors Fixed
    pot = nodePot(n,1:nStates(n));
    for e = edges(:)'
        n1 = edgeEnds(e,1);
        n2 = edgeEnds(e,2);
        
        if n == edgeEnds(e,1)
            ep = edgePot(1:nStates(n),y(n2),e).';
        else
            ep = edgePot(y(n1),1:nStates(n),e);
        end
        pot = pot .* ep;
    end
    
    % Update Objective
    NLL = NLL - log(pot(y(n))) + log(sum(pot));
    nodeBel(n, 1:nStates(n)) = pot/sum(pot);
end

