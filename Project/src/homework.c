#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)


double* VAL;

int compare(const void *a, const void *b) {
    if (VAL[*(int*)a] < VAL[*(int*)b]) return 1;
    if (VAL[*(int*)a] > VAL[*(int*)b]) return -1;
    return 0;
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    int *tab = malloc(sizeof(int) * theMesh->nodes->nNodes);

    for (i = 0; i < theMesh->nodes->nNodes; i++) 
        tab[i] = i;

    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM : 
            VAL = theMesh->nodes->X;
            qsort(tab, theMesh->nodes->nNodes, sizeof(int), compare);
        case FEM_YNUM : 
            VAL = theMesh->nodes->Y; 
            qsort(tab, theMesh->nodes->nNodes, sizeof(int), compare);
            break;            
        default : Error("Unexpected renumbering option"); }
    
    for (i = 0; i < theMesh->nodes->nNodes; i++) {
        theMesh->number[tab[i]] = i;
    }
    free(tab);
}

int femMeshComputeBand(femMesh *theMesh)
{   
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theMesh->number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; 
    }         
    return((++myBand)*2); // Pour chaque noeud, on a une composante U et une composante V
}

void conjugateGradient(double** A, double *b, double *x, int size) {
    double r[size], d[size], s[size];
    double alpha, beta, delta, deltaNew;

    // Initialisation
    for (int i = 0; i < size; i++) {
        r[i] = b[i];
        d[i] = r[i];
        x[i] = 0.0;
    }

    delta = 0.0;
    for (int i = 0; i < size; i++)
        delta += r[i] * r[i];

    int k = 0;
    while (delta > 1e-10) {  // Critère de convergence
        // Calcul de s
        for (int i = 0; i < size; i++) {
            s[i] = 0.0;
            for (int j = 0; j < size; j++)
                s[i] += A[i][j] * d[j];
        }

        // Calcul de alpha
        double sd = 0.0;
        for (int i = 0; i < size; i++)
            sd += d[i] * s[i];

        alpha = delta / sd;

        // Mise à jour de x et r
        for (int i = 0; i < size; i++) {
            x[i] += alpha * d[i];
            r[i] -= alpha * s[i];
        }

        deltaNew = 0.0;
        for (int i = 0; i < size; i++)
            deltaNew += r[i] * r[i];

        // Calcul de beta
        beta = deltaNew / delta;

        // Mise à jour de d
        for (int i = 0; i < size; i++)
            d[i] = r[i] + beta * d[i];

        delta = deltaNew;
        k++;
    }
    printf("Convergence after %d iterations.\n", k);
}

double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;

    //Incomplete Cholesky factorization 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band/2 + 1,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    // Back-substitution 

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band/2 + 1,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}

void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}

void femApplyBoundaryConditions(femProblem *theProblem) {

    femFullSystem *theSystem = theProblem->system;
    femMesh *theMesh         = theProblem->geometry->theElements;
    double **A  = theSystem->A;
    double *B   = theSystem->B;
    int iCase   = theProblem->planarStrainStress;

    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {

        femBoundaryCondition* cnd = theProblem->conditions[i];
        femMesh * bndMesh = cnd->domain->mesh;
        femBoundaryType bndType = cnd->type;
        double * X = bndMesh->nodes->X;
        double * Y = bndMesh->nodes->Y;
        int * bndElem = cnd->domain->elem;
        int nElem = cnd->domain->nElem;
        int node0, node1;

        // Pour des conditons en normales-tangents
        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T ||
                cnd->type == NEUMANN_N || cnd->type == NEUMANN_T) {

            if (cnd->domain->n_t_malloced == 0) { // si les n et t n'ont pas encore été calculées
                cnd->domain->n_t_malloced = 1;

                // nbre de noeuds = nbre de segments + 1, et pour chaque noeud, normale et tangent = vecteur de taille 2
                cnd->domain->normales = calloc(sizeof(double), 2 * (bndMesh->nElem + 1) ); 
                cnd->domain->tangentes = calloc(sizeof(double),  2 * (bndMesh->nElem + 1) ); 
                double *normales = cnd->domain->normales;
                double *tangentes = cnd->domain->tangentes;

                // Calcul des normales et tangentes
                for (int j = 0; j < nElem; j++){
                    node0 = bndMesh->elem[2 * bndElem[j]];
                    node1 = bndMesh->elem[2 * bndElem[j] + 1];
                    double dx = X[node1] - X[node0];
                    double dy = Y[node1] - Y[node0];

                    // contributions au noeud gauche
                    tangentes[2 * j]        += dx;
                    tangentes[2 * j + 1]    += dy;
                    normales[2 * j]         += - dy;
                    normales[2 * j + 1]     += dx;

                    // contributions au noeud droite
                    tangentes[2 * (j+1)]        += dx;
                    tangentes[2 * (j+1) + 1]    += dy;
                    normales[2 * (j+1)]         += - dy;
                    normales[2 * (j+1) + 1]     += dx;
                }

                // Normalisation des normales et tangentes :
                for (int j = 0; j < nElem + 1; j++){
                    double norm = sqrt( tangentes[2*j]*tangentes[2*j] + tangentes[2*j+1]*tangentes[2*j+1]);
                    tangentes[2*j] /= norm; tangentes[2*j+1] /= norm;
                    norm = sqrt( normales[2*j]*normales[2*j] + normales[2*j+1]*normales[2*j+1]);
                    normales[2*j] /= norm; normales[2*j+1] /= norm;
                }
            }
            
            if ((cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) &&
                    cnd->domain->n_t_matrix == 0) { // si A et B n'ont pas été adaptés en n - t
                cnd->domain->n_t_matrix = 1;
                // Itération sur tous les noeuds du domaine
                for (int j = 0; j < nElem + 1; j++){
                    // technique pour récupérer tous les noeuds sur base des éléments
                    if (j == nElem) {
                        node0 = bndMesh->elem[2 * bndElem[j-1] + 1];
                    }
                    else {node0 = bndMesh->elem[2 * bndElem[j]];}

                    // Combinaison linéaire des lignes et colonnes de la matrice pour changer (U, V) en (N, T) (voir rapport)
                    double A_U, A_V, B_U, B_V, nx, ny, tx, ty;
                    double *normales    = cnd->domain->normales;
                    double *tangentes   = cnd->domain->tangentes;
                    nx = normales[2 * j];   ny = normales[2 * j+1];
                    tx = tangentes[2 * j];  ty = tangentes[2 * j+1];

                    // Modification de B
                    B_U = B[2*node0];   B_V = B[2*node0+1];
                    B[2*node0]      = nx * B_U + ny * B_V;
                    B[2*node0 + 1]  = tx * B_U + ty * B_V;
                    // Modification des lignes de A
                    for (int k = 0; k < theSystem->size; k++) {
                        A_U = A[2*node0][k];    A_V = A[2*node0+1][k];
                        A[2*node0][k]   = nx *A_U + ny * A_V;
                        A[2*node0+1][k] = tx *A_U + ty * A_V;
                    }
                    // Modification des colonnes de A
                    for (int k = 0; k < theSystem->size; k++) {
                        A_U = A[k][2*node0];    A_V = A[k][2*node0+1];
                        A[k][2*node0]   = nx *A_U + ny * A_V;
                        A[k][2*node0+1] = tx *A_U + ty * A_V;
                    }
                }
            }
                   
        }

        // Application des différentes conditions
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

        if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y ||
            cnd->type == NEUMANN_N || cnd->type == NEUMANN_T) {
            for (int j = 0; j < nElem; j++) {
                node0 = bndMesh->elem[2 * bndElem[j]];
                node1 = bndMesh->elem[2 * bndElem[j] + 1];                
                double jac = 0.5 * sqrt( (X[node0] - X[node1]) *  (X[node0] - X[node1]) +
                                    (Y[node0] - Y[node1]) * (Y[node0] - Y[node1]));

                if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y) {
                    if (iCase == AXISYM) { // Prise en compte du cas axisymétrique
                        
                        // ce sont les deux valeurs prises par les fcts de formes, mais pas les fcts de formes
                        double phi[2] = { (1+sqrt(3))/2, (1-sqrt(3))/2}; 
                        double x[2] = {X[node0], X[node1]};
                        double xLoc[2] = {x[0] * phi[1] + x[1] * phi[0], x[0] * phi[0] + x[1] * phi[1]};

                        int shift = (cnd->type == NEUMANN_X ) ? 0 : 1;           
                        B[2 * node0 + shift] += jac * cnd->value * phi[1] * xLoc[0] ; // premier point d'intégration
                        B[2 * node0 + shift] += jac * cnd->value * phi[0] * xLoc[1] ; // deuxième point d'intégration

                        B[2 * node1 + shift] += jac * cnd->value * phi[0] * xLoc[0] ; // premier point d'intégration
                        B[2 * node1 + shift] += jac * cnd->value * phi[1] * xLoc[1] ; // deuxième point d'intégration                                                                             
                    }

                    else {
                        int shift = (cnd->type == NEUMANN_X ) ? 0 : 1;
                        B[2 * node0 + shift] += jac * cnd->value;
                        B[2 * node1 + shift] += jac * cnd->value;                      
                    }

                }
                else if (cnd->type == NEUMANN_N || cnd->type == NEUMANN_T){ // NEUMANN_N ou NEUMANN_T
                    double *n_or_t = (cnd->type == NEUMANN_N) ? cnd->domain->normales : cnd->domain->tangentes;
                    // l'index du noeud gauche dans n_or_t = index de l'élement
                    // l'index du noeud gauche dans n_or_t = index de l'élement + 1
                    B[2 * node0]        += jac * cnd->value *     n_or_t[2 * j];
                    B[2 * node0 + 1]    += jac * cnd->value *     n_or_t[2 * j + 1];
                    B[2 * node1]        += jac * cnd->value *     n_or_t[2 * (j+1)];
                    B[2 * node1 + 1]    += jac * cnd->value *     n_or_t[2*(j+1) + 1] ;
                }

            }
        }
    }  

    // Conditions DIRICHLET X-Y
    int *theConstrainedNodes = theProblem->constrainedNodes; // NB : ici nodes ne correspond pas à un noeud, mais à la composante x ou y
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femBoundaryType cndType = theProblem->conditions[theConstrainedNodes[i]]->type;
            if (cndType == DIRICHLET_X || cndType == DIRICHLET_Y) {
                int shift = (cndType ==DIRICHLET_X) ? 0 : 1;
                int node = (i - shift) / 2;
                int index = theMesh->number[node] * 2 + shift;
                femFullSystemConstrain(theSystem, index,value); }
            }
        }

}

void femEquationsN_Tto_U_V(femProblem *theProblem) {
    femFullSystem *theSystem = theProblem->system;
    femMesh *theMesh         = theProblem->geometry->theElements;
    double **A  = theSystem->A;
    double *B   = theSystem->B;

    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
    femBoundaryCondition* cnd = theProblem->conditions[i] ;
        if ( (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) &&
                (cnd->domain->n_t_matrix == 1) ) { // si on a pas encore converti de N-T en U-V
            cnd->domain->n_t_malloced = 0;
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

                nx = normales[2 * j];   ny = normales[2 * j+1];          
                tx = tangentes[2 * j];  ty = tangentes[2 * j+1];

                double B_U = B[2*node0];    double B_V = B[2*node0+1];
                B[2*node0]      = nx * B_U + ny * B_V;
                B[2*node0 + 1]  = tx * B_U + ty * B_V;
            }
        }
    }
}

void femFreeN_T(femProblem *theProblem) {
    // free tan - norm
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition* cnd = theProblem->conditions[i] ;
            if (cnd->domain->n_t_malloced == 1) {
                free(cnd->domain->normales);
                free(cnd->domain->tangentes);
            }
        }
}

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
    femRenumType renumType = theProblem->renumType;

    femMeshRenumber(theMesh, renumType);
    int myBand = femMeshComputeBand(theMesh);
    printf("\nSystem size : %d, band size : %d\n", theSystem->size, myBand);
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {

        for (j=0; j < nLocal; j++){
            map[j]  = theMesh->elem[iElem*nLocal+j];
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
            map[j]  = theMesh->number[map[j]];  // Pour la renumérotation
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
        }

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
                xLoc   += x[i] * phi[i];    // utilisé pour l'axisymétrique, r = xLoc
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }

            //// Assemblage ////
            if (iCase != AXISYM) {         
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

            else if (iCase == AXISYM) {  // For axisym problems, x --> r and y --> z
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
                B[mapY[i]] -= phi[i] * xLoc * g * rho * jac * weight; }} 
        }        
    }  
    femApplyBoundaryConditions(theProblem);

    // Résolution du système

    double *sol = malloc(sizeof(double) * theSystem->size); // solution avec la renumérotation
    if (theProblem->solver == SOLVEUR_PLEIN) {
        B = femFullSystemEliminate(theSystem);
        sol = memcpy(sol, B, sizeof(double)*theSystem->size);
    }
    if (theProblem->solver == SOLVEUR_BANDE) {
        femBandSystem *theBandSystem = femBandSystemCreate(theSystem->size, myBand);
        for (i = 0; i < theBandSystem->size; i++){
            int jmin = fmax(0, i - myBand/2);
            int jmax = fmin(theBandSystem->size, i + myBand/2 +1);
            for (j = jmin; j < jmax; j++) {
                theBandSystem->A[i][j] = A[i][j];
            }
            theBandSystem->B[i] = theSystem->B[i];
        }
        femBandSystemEliminate(theBandSystem);
        sol = B;
    }
    else if (theProblem->solver == GRADIENTS_CONJUGUES) {
        conjugateGradient(A, B, sol, theSystem->size);  
    }

    femEquationsN_Tto_U_V(theProblem);
    femFreeN_T(theProblem);

    // On remet les noeuds dans le bon ordre, suite à la renumérotation
    for (int i = 0; i < theMesh->nodes->nNodes; i++) {
        B[2*i] = sol[2* theMesh->number[i]];
        B[2*i + 1] = sol[2* theMesh->number[i] + 1];
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

    int l = 3;
    if (theProblem->planarStrainStress == AXISYM) { l = 4; }

    double *sigma = malloc(sizeof(double)*theNodes->nNodes*l);
    double *U = &displacements[0];
    double *V = &displacements[1];

    double xsi[nLocal], eta[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    femDiscreteXsi2(theSpace, xsi, eta);

    for (int i = 0; i < theNodes->nNodes; i++) {  // For each node, we compute the stress

        int *nextElem = malloc(sizeof(int)*10);  // Array containing index of the first node of the elements which have a common node with current node i
        double *epsilon = malloc(sizeof(double)*l);
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
                epsilon[2] += (dphidx[k] * vNode + dphidy[k] * uNode) / 2;
            }
        }
        sigma[i*l] = a*epsilon[0] + b*epsilon[1];
        sigma[i*l + 1] = a*epsilon[1] + b*epsilon[0];
        sigma[i*l + 2] = 2*c*epsilon[2];

        if (theProblem->planarStrainStress == AXISYM) {
            epsilon[3] += U[2*theMesh->elem[i]] / theNodes->X[i];
            sigma[i*l] += b*epsilon[3];
            sigma[i*l + 1] += b*epsilon[3];
            sigma[i*l + 3] = b*(epsilon[0] + epsilon[1]) + a*epsilon[3];
        }
        
        free(epsilon);
        free(nextElem);
    }
    return sigma;
}

double *femPlastic(femProblem *theProblem, double *sigma) {
    // Compute the von Mises equivalent constraint
    int l = 3;
    if (theProblem->planarStrainStress == AXISYM) { l = 4; }
    int n = theProblem->geometry->theNodes->nNodes;

    double *sigEq = malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++) {
        double sigtt = 0;
        if (l==4) {
            sigtt = sigma[l*i+3]; }
        // Von Mises
        sigEq[i] = sqrt((pow(sigma[l*i]-sigma[l*i+1], 2) + pow(sigma[l*i+1]-sigtt,2) + pow(sigtt-sigma[i*l],2))/2 + 3*pow(sigma[i*l+2],2));
    }
    double sigMax = femMax(sigEq, n);
    printf("\nLa limite de plasticité est de %14.7e, la contrainte équivalente maximale : %14.7e -> %14.7f %% \n", theProblem->sigmaY, sigMax, sigMax/theProblem->sigmaY); // %% interprété comme % par printf
    return sigEq;
}

void femPrintStress(double *stress, int nNodes, int l) {
    printf("\n ----------------- Stress -----------------\n\n");
    if (l == 3) {
        for (int i = 0; i < nNodes*l; i+=l) {
            printf("%d : xx %14.7e, xy %14.7e, yy %14.7e\n", i/l, stress[i], stress[i + 2], stress[i + 1]);
        }
    }
    else {
        for (int i = 0; i < nNodes*l; i+=l) {
            printf("%d : rr %14.7e, zz %14.7e, rz %14.7e, qq %14.7e\n", i/l, stress[i], stress[i + 1], stress[i + 2], stress[i+4]);
        }
    }
}
