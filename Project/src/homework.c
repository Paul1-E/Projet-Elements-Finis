#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)


/*double *theGlobalCoord;

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i, *inverse;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nNode; i++) 
                theMesh->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*theMesh->nNode);
            for (i = 0; i < theMesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->X;
            qsort(inverse, theMesh->nNode, sizeof(int), compare);
            for (i = 0; i < theMesh->nNode; i++)
                theMesh->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*theMesh->nNode);
            for (i = 0; i < theMesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->Y;
            qsort(inverse, theMesh->nNode, sizeof(int), compare);
            for (i = 0; i < theMesh->nNode; i++)
                theMesh->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}


int femMeshComputeBand(femMesh *theMesh)
{
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theMesh->number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}


double  *femBandSystemEliminate(femFullSystem *myBand, int band)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;

    /* Incomplete Cholesky factorization 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution 

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}
*/


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    int iCase = theProblem->planarStrainStress;
    int nLocal = theMesh->nLocalNode;

    double x[nLocal],y[nLocal],phi[nLocal],dphidxsi[nLocal],dphideta[nLocal],dphidx[nLocal],dphidy[nLocal];
    int iElem,iInteg,iEdge,i,j,d,map[nLocal],mapX[nLocal],mapY[nLocal];

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    //femMeshRenumber(theMesh, FEM_XNUM);
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            double xLoc = 0.0;
            
            for (i = 0; i < theSpace->n; i++) {
                xLoc   += x[i] * phi[i];    // utilisé pour l'axisymétrique
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }

            //// Assemblage ////
            // For axisym problems, we consider x --> r and y --> z
            // Therefore, x > 0, and gravity is G = -g e_z
            // and r = xLoc
            if (iCase == AXISYM) { 
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * xLoc + 
                                            dphidy[i] * c * dphidy[j] * xLoc +
                                            phi[i] * (b * dphidx[j] + a * phi[j]/xLoc) +
                                            dphidx[i] * b * phi[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * xLoc + 
                                            dphidy[i] * c * dphidx[j] * xLoc + 
                                            phi[i] * b * dphidy[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * xLoc + 
                                            dphidx[i] * c * dphidy[j] * xLoc +
                                            dphidy[i] * b * phi[j]) * jac * weight;                                                                                             
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * xLoc + 
                                            dphidx[i] * c * dphidx[j] * xLoc) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * xLoc * g * rho * jac * weight; }
            }

            else {         
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}   
        }        
    }

    // Conditions frontières
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {

        femBoundaryCondition* cnd = theProblem->conditions[i];
        femMesh * bndMesh = cnd->domain->mesh;
        femBoundaryType bndType = cnd->type;
        double * X = bndMesh->nodes->X;
        double * Y = bndMesh->nodes->Y;
        int * bndElem = cnd->domain->elem;
        int nElem = cnd->domain->nElem;
        int node0, node1;

        // Calcul des normales et tangentes
        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T ||
                cnd->type == NEUMANN_N || cnd->type == NEUMANN_T) {

            if (cnd->domain->n_t_determined == 0) { // si les n et t n'ont pas encore été calculées
                cnd->domain->n_t_determined = 1;
                // nbre de noeuds = nbre de segments + 1
                // pour chaque noeud, normale et tangent = vecteur de taille 2
                cnd->domain->normales = calloc(sizeof(double), 2 * (bndMesh->nElem + 1) ); 
                double *normales = cnd->domain->normales;
                cnd->domain->tangentes = calloc(sizeof(double),  2 * (bndMesh->nElem + 1) ); 
                double *tangentes = cnd->domain->normales;

                // Calcul des normales et tangentes
                for (int j = 0; j < nElem; j++){
                    node0 = bndMesh->elem[2 * bndElem[j]];
                    node1 = bndMesh->elem[2 * bndElem[j] + 1];

                    // contributions au noeud gauche
                    tangentes[2 * j] += X[node1] - X[node0];
                    tangentes[2 * j + 1] += Y[node1] - Y[node0];
                    normales[2 * j] += - Y[node1] + Y[node0];
                    normales[2 * j + 1] += X[node1] - X[node0];

                    // contributions au noeud droite
                    tangentes[2 * (j+1)] += X[node1] - X[node0];
                    tangentes[2 * (j+1) + 1] += Y[node1] - Y[node0];
                    normales[2 * (j+1)] += - Y[node1] + Y[node0];
                    normales[2 * (j+1) + 1] += X[node1] - X[node0];

                }

                // Normalisation des normales et tangentes :
                for (int j = 0; j < nElem + 1; j++){
                    double norm = sqrt( tangentes[2*j]*tangentes[2*j] + tangentes[2*j+1]*tangentes[2*j+1]);
                    tangentes[2*j] /= norm;
                    tangentes[2*j+1] /= norm;
                    norm = sqrt( normales[2*j]*normales[2*j] + normales[2*j+1]*normales[2*j+1]);
                    normales[2*j] /= norm;
                    normales[2*j+1] /= norm;
                }

                // petit check
                // for (int k = 0; k < nElem + 1; k++ ) {
                //     printf("normale du noeud %d = (%f ; %f)\n", k, normales[2 * k], normales[2 * k + 1]);   
                //     printf("tangente du noeud %d = (%f ; %f)\n", k, tangentes[2 * k], tangentes[2 * k + 1]);}   

                // Itération sur tous les noeuds du domaine
                for (int j = 0; j < nElem + 1; j++){

                    // technique pour récupérer tous les noeuds sur base des éléments
                    if (j == nElem) {
                        node0 = bndMesh->elem[2 * bndElem[j-1] + 1];
                    }
                    else {node0 = bndMesh->elem[2 * bndElem[j]];}
                    

                    // Combinaison linéaire des lignes de la matrice pour changer (U, V) en (N, T) (voir rapport)
                    double A_U, A_V, B_U, B_V, nx, ny, tx, ty;
                    nx = normales[2 * j];
                    ny = normales[2 * j+1];
                    tx = tangentes[2 * j];
                    ty = tangentes[2 * j+1];

                    B_U = B[2*node0];
                    B_V = B[2*node0+1];
                    B[2*node0] = nx * B_U + ny * B_V;
                    B[2*node0 + 1] = tx * B_U + ty * B_V;
                    for (int k = 0; k < theSystem->size; k++) {
                        A_U = A[2*node0][k];
                        A_V = A[2*node0+1][k];
                        A[2*node0][k]   = nx *A_U + ny * A_V;
                        A[2*node0+1][k] = tx *A_U + ty * A_V;
                    }
                }
            }       
        }  


        // Application des différentes contraintes

        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) {
            for (int j = 0; j < nElem + 1; j++){
                if (j == nElem) {
                    node0 = bndMesh->elem[2 * bndElem[j-1] + 1];
                }
                else {node0 = bndMesh->elem[2 * bndElem[j]];}       
                int shift = (cnd->type == DIRICHLET_N) ? 0 : 1;
                femFullSystemConstrain(theSystem, 2 * node0 + shift,  cnd->value);
            }
        }

        if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y) {
            for (int j = 0; j < nElem; j++){
                node0 = bndMesh->elem[2 * bndElem[j]];
                node1 = bndMesh->elem[2 * bndElem[j] + 1];                
                double jac = 0.5 * sqrt( (X[node0] - X[node1]) *  (X[node0] - X[node1]) +
                                    (Y[node0] - Y[node1]) * (Y[node0] - Y[node1]));

                if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y) {
                    int shift = (cnd->type == NEUMANN_X ) ? 0 : 1;
                    B[2 * node0 + shift] += jac * cnd->value;
                    B[2 * node1 + shift] += jac * cnd->value;
                }
                else { // NEUMANN_N ou NEUMANN_T
                    double *n_or_t = (cnd->type == NEUMANN_N) ? cnd->domain->normales : cnd->domain->tangentes;
                    // l'index du noeud gauche dans n_or_t = index de l'élement
                    // l'index du noeud gauche dans n_or_t = index de l'élement + 1
                    B[2 * node0] += jac * cnd->value *      n_or_t[2 * j];
                    B[2 * node0 + 1] += jac * cnd->value *  n_or_t[2 * j + 1];
                    B[2 * node1] += jac * cnd->value *      n_or_t[2 * (j+1)];
                    B[2 * node1 + 1] += jac * cnd->value *  n_or_t[2*(j+1) + 1] ;
  
                }

            }
        }
    }  

    int *theConstrainedNodes = theProblem->constrainedNodes; // NB : ici nodes ne correspond pas à un noeud, mais à la composante x ou y
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femBoundaryType constraintype = theProblem->conditions[theConstrainedNodes[i]]->type;
            if (constraintype == DIRICHLET_X || constraintype == DIRICHLET_Y) {
            femFullSystemConstrain(theSystem,i,value); }
            }}

    /*int band = femMeshComputeBand(theMesh);
    return femBandSystemEliminate(theSystem, band);*/

    B = femFullSystemEliminate(theSystem);  
    // CONDITIONS DIRICHLET N-T : On convertit les équations (N,T) en (U,V)
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition* cnd = theProblem->conditions[i] ;
        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) {
            int nElem = cnd->domain->nElem;
            femMesh * bndMesh = cnd->domain->mesh;
            int * bndElem = cnd->domain->elem;
            int node0, nx, ny, tx, ty;
            double *normales = cnd->domain->normales;
            double *tangentes = cnd->domain->tangentes;

            for (int j = 0; j < nElem + 1; j++){
                if (j == nElem) {
                    node0 = bndMesh->elem[2 * bndElem[j-1] + 1];
                }
                else {node0 = bndMesh->elem[2 * bndElem[j]];}

                nx = normales[2 * j];
                ny = normales[2 * j+1];          
                tx = tangentes[2 * j];
                ty = tangentes[2 * j+1];


                double B_U = B[2*node0];
                double B_V = B[2*node0+1];
                B[2*node0] = nx * B_U + ny * B_V;
                B[2*node0 + 1] = tx * B_U + ty * B_V;
            }
        }
    }

    return B;
}


double *femFindStress(femProblem *theProblem, double *displacements) {

    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C; 
    int nLocal = theMesh->nLocalNode;

    double *sigma = malloc(sizeof(double)*theNodes->nNodes*4);
    double *U = &displacements[0];
    double *V = &displacements[1];

    double xsi[nLocal], eta[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    femDiscreteXsi2(theSpace, xsi, eta);

    for (int i = 0; i < theNodes->nNodes; i++) {  // For each node, we compute the stress

        int *nextElem = malloc(sizeof(int)*10);  // Array containing index of the first node of the elements which have a common node with current node i
        double *epsilon = malloc(sizeof(double)*4);
        int idxNext = 0;
        int c;

        for (int j = 0; j < theMesh->nElem*nLocal; j++) {  // Searching for elements next to current node
            if (theMesh->elem[j] == i) {
                c = j%nLocal;
                nextElem[idxNext] = j - c;
                idxNext++;
            }
        }
        /*if (i == 0) {  TEST OK
            for (int j = 0; j < idxNext; j++) {
                printf("\n%d", nextElem[j]);
                printf("\n%d\n", theMesh->elem[nextElem[j]]);
            }
        }*/
        
        for (int j = 0; j < idxNext; j++) {  // For each element j next to the current node

            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;

            for (int k = 0; k < nLocal; k++) {  
                femDiscreteDphi2(theSpace, xsi[c], eta[c], dphidxsi, dphideta);
                dxdxsi += theNodes->X[theMesh->elem[nextElem[j]+k]]*dphidxsi[k];       
                dxdeta += theNodes->X[theMesh->elem[nextElem[j]+k]]*dphideta[k];   
                dydxsi += theNodes->Y[theMesh->elem[nextElem[j]+k]]*dphidxsi[k];   
                dydeta += theNodes->Y[theMesh->elem[nextElem[j]+k]]*dphideta[k];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (int k = 0; k < nLocal; k++) {
                dphidx[k] = (dphidxsi[k] * dydeta - dphideta[k] * dydxsi) / jac;       
                dphidy[k] = (dphideta[k] * dxdxsi - dphidxsi[k] * dxdeta) / jac; 
            }

            for (int k = 0; k < nLocal; k++) {

                double uNode = U[2*theMesh->elem[nextElem[j]+k]];
                double vNode = V[2*theMesh->elem[nextElem[j]+k]];

                epsilon[0] += dphidx[k] * uNode;
                epsilon[1] += dphidy[k] * vNode;
                epsilon[2] += dphidx[k] * vNode;
                epsilon[3] += dphidy[k] * uNode;
            }
        }

        sigma[i*4] = a*epsilon[0] + b*epsilon[1];
        sigma[i*4 + 1] = a*epsilon[1] + b*epsilon[0];
        sigma[i*4 + 2] = 2*c*epsilon[2];
        sigma[i*4 + 3] = 2*c*epsilon[3];
        
        free(epsilon);
        free(nextElem);
    }
    return sigma;
}


void femPrintStress(double *stress, int nNodes) {
    printf("\n ----------------- Stresses -----------------\n\n");
    for (int i = 0; i < nNodes*4; i+=4) {
        printf("%d : xx %14.7e, xy %14.7e, yx %14.7e, yy %14.7e\n", i/4, stress[i], stress[i + 2], stress[i + 3], stress[i + 1]);
    }
}