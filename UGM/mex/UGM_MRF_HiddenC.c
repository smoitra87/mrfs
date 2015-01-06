#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Variables */
    
    int n, s, n1, n2, s1, s2, i,e,
            nInstances, nNodes, nEdges, maxState, 
            *edgeEnds, *nStates, *nodeMap, *edgeMap, *Y;
    double obs, *g, *nodeBel, *nodeBelC, *edgeBel, *edgeBelC, *nodePot, 
           *edgePot;
    double logZC_update = 0;

    int sizeGrad[2];

    /* Input */
    g = mxGetPr(prhs[0]);
    Y = (int*)mxGetPr(prhs[1]);
    i = (int)mxGetScalar(prhs[2])-1;
    nodePot = mxGetPr(prhs[3]);
    edgePot = mxGetPr(prhs[4]);
    nodeMap = (int*)mxGetPr(prhs[5]);
    edgeMap = (int*)mxGetPr(prhs[6]);
    nodeBelC = mxGetPr(prhs[7]);
    edgeBelC = mxGetPr(prhs[8]);
    nodeBel = mxGetPr(prhs[9]);
    edgeBel = mxGetPr(prhs[10]);
    edgeEnds = (int*)mxGetPr(prhs[11]);
    nStates = (int*)mxGetPr(prhs[12]);

    /* Calculate derivative values */
    nNodes = mxGetDimensions(prhs[3])[0];
    nEdges = mxGetDimensions(prhs[11])[0];
    maxState = getMaxState(nStates, nNodes);
    nInstances = mxGetDimensions(prhs[1])[0];

	if (!mxIsClass(prhs[1],"int32")||!mxIsClass(prhs[2],"int32")
            ||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32")
            ||!mxIsClass(prhs[11],"int32")||!mxIsClass(prhs[12],"int32"))
		mexErrMsgTxt("Y,i,nodeMap, edgeMap, edgeEnds, nStates must be int32");

     for(n = 0; n < nNodes; n++) {
        if(Y[i + nInstances*n] > 0) {
            logZC_update += log(nodePot[n + nNodes*(Y[i+nInstances*n]-1)]);
        }
    }
     
    for(e = 0; e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        
        if(Y[i + nInstances*n2] > 0 && Y[i + nInstances*n2] >0)  {
            logZC_update += log(edgePot[Y[i + nInstances*n1]-1 + 
                    maxState*(Y[i + nInstances*n2]-1  + maxState*e)]);
        }
    }

    plhs[0] = mxCreateDoubleScalar(logZC_update);

    for(n = 0; n < nNodes; n++) {
        for(s = 0; s < nStates[n]; s++) {
            if(nodeMap[n + nNodes*s] > 0) {
                if(Y[i + nInstances*n] == 0){
                    obs = nodeBelC[n + nNodes*s];
                }
                else if(s == Y[i + nInstances*n]-1) {
                    obs = 1;
                }
                else {
                    obs = 0;
                }
                g[nodeMap[n + nNodes*s]-1] += nodeBel[n + nNodes*s] - obs;
            }
        }
    }
 
    for(e = 0; e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        
        for (s1 = 0; s1 < nStates[n1]; s1++) {
            for (s2 = 0; s2 < nStates[n2]; s2++) {
                if (edgeMap[s1 + maxState*(s2 + maxState*e)] > 0) {
                    if (Y[i + nInstances*n1] == 0 && Y[i + nInstances*n2] == 0) {
                        obs = edgeBelC[s1 + maxState*(s2 + maxState*e)] ;
                    }
                    else if (Y[i + nInstances*n1]==0 && s2 == Y[i + nInstances*n2]-1) {
                        obs = edgeBelC[s1 + maxState*(Y[i + nInstances*n2]-1 + maxState*e)]; 
                    }
                    else if (s1 == Y[i + nInstances*n1]-1 && Y[i + nInstances*n2]==0) {
                        obs = edgeBelC[Y[i + nInstances*n1]-1 +
                            maxState*(s2 + maxState*e)]; 
                    }
                    else if (s1 == Y[i + nInstances*n1]-1 && s2 == Y[i + nInstances*n2]-1) {
                        obs = 1; 
                    }
                    else {
                        obs = 0;
                    }
                    g[edgeMap[s1 + maxState*(s2 + maxState*e)]-1] 
                        += edgeBel[s1 + maxState*(s2 + maxState*e)] - obs;
                }
            }
        }
    }
}
