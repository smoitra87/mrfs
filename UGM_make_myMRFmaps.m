function [nodeMap,edgeMap] = UGM_make_hiddenMRFmaps(edgeStruct)

nEdges = edgeStruct.nEdges;
nNodes = edgeStruct.nNodes;
nStates = edgeStruct.nStates;
maxState = max(nStates);

nodeMap = zeros(nNodes,maxState,'int32');
fNum = 0;
for n = 1:nNodes
    for s = 1:nStates(n)
        fNum = fNum+1;
        nodeMap(n,s) = fNum;
    end
end
nNodeParams = fNum;

edgeMap = zeros(maxState,maxState,nEdges,'int32');
sse = 1;
for e = 1:nEdges
    n1 = edgeStruct.edgeEnds(e,1);
    n2 = edgeStruct.edgeEnds(e,2);
    for s1 = 1:nStates(n1)
        for s2 = 1:nStates(n2)
            edgeMap(s1,s2,e) = nNodeParams+sse;
            sse = sse+1;
        end
    end
end
