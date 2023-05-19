 
void largeur_de_bande(){ 
    // check de la largeur de bande après assemblage
    int largeurBande = 0;
    // Parcours de la matrice pour trouver la largeur de la bande
    int size = theSystem->size;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (A[i][j] != 0) {
                int distanceDiagonale = i - j;
                int distanceAbsolue = distanceDiagonale >= 0 ? distanceDiagonale : -distanceDiagonale;
                if (distanceAbsolue > largeurBande) {
                    largeurBande = distanceAbsolue;
                }
            }
        }
    }
    largeurBande++;
    printf("La largeur de la bande de la matrice est : %d\n", largeurBande);    
}

void print_normales_tangentes() {
        // petit check
    // for (int k = 0; k < nElem + 1; k++ ) {
    //     printf("normale du noeud %d = (%f ; %f)\n", k, normales[2 * k], normales[2 * k + 1]);   
    //     printf("tangente du noeud %d = (%f ; %f)\n", k, tangentes[2 * k], tangentes[2 * k + 1]);}  
}