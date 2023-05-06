#include "fem.h"




void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();


    double h = 0.015;

    int ierr;
    int idRect1 = gmshModelOccAddRectangle(0.0, 0.0, 0.0,   0.03, 0.6,   -1,0.0,&ierr); 
    int idRect2 = gmshModelOccAddRectangle(0.0, 0.57, 0.0,   0.65, 0.03,   -1,0.0,&ierr); 
    int rect1[] = {2,idRect1};
    int rect2[] = {2,idRect2};

    gmshModelOccFuse(rect1, 2, rect2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
      
        // get the local infos
        for (i = 0; i < nLocal; i++) {
            map[i] = theMesh->elem[nLocal * iElem + i];
            x[i] = theGeometry->theNodes->X[map[i]];
            y[i] = theGeometry->theNodes->Y[map[i]];
        }

        // we iterate on the integration points
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            // calcul des fonctions de formes et leurs dérivées
            femDiscretePhi2(theSpace, xsi, eta, phi); 
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0;      double dydxsi = 0;
            double dxdeta = 0;      double dydeta = 0;

            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];       dydxsi += y[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];       dydeta += y[i] * dphideta[i];
            }

            // calcul de la jacobienne
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi); 

            // calcul de dphidx et dphidy
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dydeta * dphidxsi[i] - dydxsi * dphideta[i]) / jac;
                dphidy[i] = (dxdxsi * dphideta[i] - dxdeta * dphidxsi[i]) / jac;
            }

            // Assemblage
            for (i = 0; i < theProblem->space->n;  i++) {
                double wtimesj = weight * jac;

                B[map[i] * 2    ] += 0;
                B[map[i] * 2 + 1] += - wtimesj * (phi[i] * rho * g);

                for (j = 0; j < theProblem->space->n; j++) {

                    A[map[i] * 2    ][map[j] * 2    ] += wtimesj * (a * dphidx[i] * dphidx[j] + c * dphidy[i] * dphidy[j]);;
                    A[map[i] * 2 + 1][map[j] * 2    ] += wtimesj * (b * dphidy[i] * dphidx[j] + c * dphidx[i] * dphidy[j]);;
                    A[map[i] * 2    ][map[j] * 2 + 1] += wtimesj * (b * dphidx[i] * dphidy[j] + c * dphidy[i] * dphidx[j]);;
                    A[map[i] * 2 + 1][map[j] * 2 + 1] += wtimesj * (c * dphidx[i] * dphidx[j] + a * dphidy[i] * dphidy[j]);;
                }
            }
        }
    }            
                
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
