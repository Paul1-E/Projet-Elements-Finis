#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade 
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


double geoSize(double x, double y) {

    femGeo* theGeometry = geoGetGeometry();
    double f = 2*(exp(-pow(x-theGeometry->LxPlate/2, 2)/(theGeometry->LxPlate*1.5)) + exp(-pow(y-theGeometry->LyPlate, 2)/(theGeometry->LyPlate*1.5)));
    f += 2*(exp(-pow(x+theGeometry->LxPlate/2, 2)/(theGeometry->LxPlate*1.5)) + exp(-pow(y, 2)/(theGeometry->LyPlate*1.5)));
    return theGeometry->h/f;
}



void geoMeshGenerate_standard() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    double dh = h/12;
    double dw = w/20;
    
    int ierr;

    int idPt0, idPt1, idPt2, idPt3, idSeg0, idSeg1, idSeg2, idSeg3, idWire;

    // Surface 1
    idPt0= gmshModelOccAddPoint(0.0, 0.0, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.8*dw, 0.0, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-1.2*dw, 11*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-2*dw, 11*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    int seg[] = {idSeg0, idSeg1, idSeg2, idSeg3};
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface1 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface1[] = {2, idSurface1};

    // Surface 2
    idPt0= gmshModelOccAddPoint(-1.5*dw, 9.3*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(8.5*dw, 11.5*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8*dw, 12*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-1.5*dw, 10*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface2 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface2[] = {2, idSurface2};

    // Surface 3
    idPt0= gmshModelOccAddPoint(0.0, 0.2*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.8*dw, 0.0, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(9*dw, 9*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(9*dw, 10*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface3 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface3[] = {2, idSurface3};


    // Surface 4
    idPt0= gmshModelOccAddPoint(10*dw, 3*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(10.7*dw, 3*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8.7*dw, 12*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(8*dw, 12*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface4 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface4[] = {2, idSurface4};

    // Surface 5
    idPt0= gmshModelOccAddPoint(-8*dw, 2.2*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(-1.2*dw, 10*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-1.5*dw, 10.5*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 3*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface5 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface5[] = {2, idSurface5};

    // Surface 6
    idPt0= gmshModelOccAddPoint(-8*dw, 2*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.0, 0.0, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(0.0, 0.8*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 2.5*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface6 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface6[] = {2, idSurface6};

    // Surface 7
    idPt0= gmshModelOccAddPoint(-8*dw, 1.7*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(-7.5*dw, 1.7*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-7.5*dw, 3*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 3*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface7 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface7[] = {2, idSurface7};

    // Fuse
    gmshModelOccFuse(surface1, 2, surface2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface3, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface4, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface5, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface6, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface7, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);


    // int idRect1 = gmshModelOccAddRectangle(0.0, h - w/10,0.0,   h,w/10,     -1,0.0,&ierr);  
    // int rect0[] = {2, idRect0}; int rect1[] = {2, idRect1};
    // gmshModelOccFuse(rect0, 2, rect1, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    // const int idPt0= gmshModelOccAddPoint(0.0, 0.0, 0.0,      -1,-1,&ierr);
    // const int idPt1= gmshModelOccAddPoint(w/10, 0.0, 0.0,     -1,-1,&ierr);
    // const int idPt2= gmshModelOccAddPoint(h, h-w/10, 0.0,     -1,-1,&ierr);
    // const int idPt3= gmshModelOccAddPoint(h, h, 0.0,          -1,-1,&ierr);
    // int idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    // int idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    // int idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    // int idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    // int seg[] = {idSeg0, idSeg1, idSeg2, idSeg3};
    // int idwire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);
    // int idRect2 = gmshModelOccAddPlaneSurface(&idwire, 1, -1, &ierr);
    // int rect2[] = {2, idRect2};

    // gmshModelOccFuse(rect0, 2, rect2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    
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



void geoMeshGenerate_standard_small(){

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    double dh = h/12;
    double dw = w/20;
    
    int ierr;

    int idPt0, idPt1, idPt2, idPt3, idSeg0, idSeg1, idSeg2, idSeg3, idWire;

    // Surface 1
    idPt0= gmshModelOccAddPoint(0.0, 0.0, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.8*dw, 0.0, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-1.2*dw, 11*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-2*dw, 11*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    int seg[] = {idSeg0, idSeg1, idSeg2, idSeg3};
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface1 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface1[] = {2, idSurface1};

    // Surface 2
    idPt0= gmshModelOccAddPoint(-1.5*dw, 9.3*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(8.5*dw, 11.5*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8*dw, 12*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-1.5*dw, 10*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface2 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface2[] = {2, idSurface2};

    // Surface 3
    idPt0= gmshModelOccAddPoint(0.0, 0.2*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.8*dw, 0.0, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(9*dw, 9*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(9*dw, 10*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface3 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface3[] = {2, idSurface3};


    // Surface 4
    idPt0= gmshModelOccAddPoint(10*dw, 3*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(10.7*dw, 3*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8.7*dw, 12*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(8*dw, 12*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface4 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface4[] = {2, idSurface4};

    // Surface 5
    idPt0= gmshModelOccAddPoint(-8*dw, 2.4*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(-1.4 * dw, 10*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-1.5*dw, 10.5*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 3*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface5 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface5[] = {2, idSurface5};

    // Surface 6
    idPt0= gmshModelOccAddPoint(-8*dw, 2.1*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.0, 0.0, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(0.0, 0.5*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 2.5*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface6 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface6[] = {2, idSurface6};

    // Surface 7
    idPt0= gmshModelOccAddPoint(-8*dw, 1.7*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(-7.5*dw, 1.7*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-7.5*dw, 3*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 3*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface7 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface7[] = {2, idSurface7};

    // Fuse
    gmshModelOccFuse(surface1, 2, surface2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface3, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface4, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface5, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface6, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface7, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);


    // int idRect1 = gmshModelOccAddRectangle(0.0, h - w/10,0.0,   h,w/10,     -1,0.0,&ierr);  
    // int rect0[] = {2, idRect0}; int rect1[] = {2, idRect1};
    // gmshModelOccFuse(rect0, 2, rect1, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    // const int idPt0= gmshModelOccAddPoint(0.0, 0.0, 0.0,      -1,-1,&ierr);
    // const int idPt1= gmshModelOccAddPoint(w/10, 0.0, 0.0,     -1,-1,&ierr);
    // const int idPt2= gmshModelOccAddPoint(h, h-w/10, 0.0,     -1,-1,&ierr);
    // const int idPt3= gmshModelOccAddPoint(h, h, 0.0,          -1,-1,&ierr);
    // int idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    // int idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    // int idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    // int idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    // int seg[] = {idSeg0, idSeg1, idSeg2, idSeg3};
    // int idwire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);
    // int idRect2 = gmshModelOccAddPlaneSurface(&idwire, 1, -1, &ierr);
    // int rect2[] = {2, idRect2};

    // gmshModelOccFuse(rect0, 2, rect2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    
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


void geoMeshGenerate_BMX() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    double dh = h/12;
    double dw = w/20;
    
    int ierr;

    int idPt0, idPt1, idPt2, idPt3, idSeg0, idSeg1, idSeg2, idSeg3, idWire;

    // Surface 1
    idPt0= gmshModelOccAddPoint(0.0, 4 * dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.8*dw, 4 * dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-1.2*dw, 11*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-2*dw, 11*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    int seg[] = {idSeg0, idSeg1, idSeg2, idSeg3};
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface1 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface1[] = {2, idSurface1};

    // Surface 2
    idPt0= gmshModelOccAddPoint(-1.5*dw, 9.4*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(8.5*dw, 11.6*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8*dw, 12*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-1.5*dw, 10*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface2 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface2[] = {2, idSurface2};

    // Surface 3
    idPt0= gmshModelOccAddPoint(0.0, 4.6*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.5*dw, 4.2 * dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8.5*dw, 10.5*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(8.5*dw, 11.5*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface3 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface3[] = {2, idSurface3};


    // Surface 4
    idPt0= gmshModelOccAddPoint(11*dw, 2.5*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(11.7*dw, 2.5*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(8.7*dw, 12*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(8*dw, 12*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface4 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface4[] = {2, idSurface4};

    // Surface 5
    idPt0= gmshModelOccAddPoint(-8*dw, 2.3*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(-1.3*dw, 9.5*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-1.5*dw, 10*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 3*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface5 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface5[] = {2, idSurface5};

    // Surface 6
    idPt0= gmshModelOccAddPoint(-8*dw, 2.1*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(0.1 * dw, 4.3*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(0.1*dw, 4.8*dh , 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 2.5*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface6 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface6[] = {2, idSurface6};

    // Surface 7
    idPt0= gmshModelOccAddPoint(-8*dw, 2*dh, 0.0,      -1,-1,&ierr);
    idPt1= gmshModelOccAddPoint(-7.5*dw, 2*dh, 0.0,     -1,-1,&ierr);
    idPt2= gmshModelOccAddPoint(-7.5*dw, 3*dh, 0.0,     -1,-1,&ierr);
    idPt3= gmshModelOccAddPoint(-8*dw, 3*dh, 0.0,          -1,-1,&ierr);

    idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    seg[0] = idSeg0, seg[1] = idSeg1; seg[2] = idSeg2; seg[3] = idSeg3;
    idWire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);

    int idSurface7 = gmshModelOccAddPlaneSurface(&idWire, 1, -1, &ierr);
    int surface7[] = {2, idSurface7};

    // Fuse
    gmshModelOccFuse(surface1, 2, surface2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface3, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface4, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface5, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface6, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccFuse(surface1, 2, surface7, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);


    // int idRect1 = gmshModelOccAddRectangle(0.0, h - w/10,0.0,   h,w/10,     -1,0.0,&ierr);  
    // int rect0[] = {2, idRect0}; int rect1[] = {2, idRect1};
    // gmshModelOccFuse(rect0, 2, rect1, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    // const int idPt0= gmshModelOccAddPoint(0.0, 0.0, 0.0,      -1,-1,&ierr);
    // const int idPt1= gmshModelOccAddPoint(w/10, 0.0, 0.0,     -1,-1,&ierr);
    // const int idPt2= gmshModelOccAddPoint(h, h-w/10, 0.0,     -1,-1,&ierr);
    // const int idPt3= gmshModelOccAddPoint(h, h, 0.0,          -1,-1,&ierr);
    // int idSeg0 = gmshModelOccAddLine(idPt0, idPt1, -1,&ierr);
    // int idSeg1 = gmshModelOccAddLine(idPt1, idPt2, -1,&ierr);
    // int idSeg2 = gmshModelOccAddLine(idPt2, idPt3, -1,&ierr);
    // int idSeg3 = gmshModelOccAddLine(idPt3, idPt0, -1,&ierr);
    // int seg[] = {idSeg0, idSeg1, idSeg2, idSeg3};
    // int idwire = gmshModelOccAddWire(seg, 4, -1, 0, &ierr);
    // int idRect2 = gmshModelOccAddPlaneSurface(&idwire, 1, -1, &ierr);
    // int rect2[] = {2, idRect2};

    // gmshModelOccFuse(rect0, 2, rect2, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    
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

void geoMeshGenerate_tandem(){}



void geoMeshGenerateGeo() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(-w/2, -h/2, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( w/2, -h/2, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( w/2,  h/2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-w/2,  h/2, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(-w/2,    r, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.,      r, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.,     -r, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(-w/2,   -r, 0., lc, 8, &ierr);
    int p9 = gmshModelGeoAddPoint(0.,     0., 0., lc, 9, &ierr); // center of circle


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


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

 //   gmshFltkRun(&ierr);
}
