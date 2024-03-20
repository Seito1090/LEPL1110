#include "fem.h"


double InterpolationHermite(double d, double *dLoc, double *fLoc, double *dfLoc){
    double t = (d - dLoc[0])/(dLoc[1] - dLoc[0]);
    double h00 = 2 * t*t*t - 3*t*t + 1;
    double h10 = t*t*t - 2*t*t + t;
    double h01 = -2*t*t*t + 3*t*t; 
    double h11 = t*t*t - t*t;
    return h00*fLoc[0] + h10*(dLoc[1] - dLoc[0])*dfLoc[0] + h01*fLoc[1] + h11*(dLoc[1] - dLoc[0])*dfLoc[1];
}
double min(double a, double b){
    if(a>b){
        return b;
    }
    return a;
}

double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


//
//     A modifier !
//     
// Your contribution starts here ....
//
    
    double interpolation;
    double distnotch = sqrt(pow(x - x0, 2) + pow(y - y0, 2)) - r0;
    double disthole = sqrt(pow(x - x1, 2) + pow(y - y1, 2)) - r1;
    
    double dfLoc[2] = {0.0, 0.0};  // valeurs de f' aux CL

    if ( disthole < d1 && distnotch < d0){ // intersection des deux champs
      double dLoc0[2] = {0.0, d0}; // points conditions limites
        double fLoc0[2] = {h0, h}; // valeurs de f aux CL
        double dLoc1[2] = {0.0, d1}; // points conditions limites
        double fLoc1[2] = {h1, h}; // valeurs de f aux CL
        double HInt0 = InterpolationHermite(distnotch, dLoc0, fLoc0, dfLoc);
        double HInt1 = InterpolationHermite(disthole, dLoc1, fLoc1, dfLoc);
        interpolation = min(min(HInt0, HInt1),h);
    }
    else if ( distnotch < d0 ){
      double dLoc0[2] = {0.0, d0}; // points conditions limites
        double fLoc0[2] = {h0, h}; // valeurs de f aux CL
        double HInt = InterpolationHermite(distnotch, dLoc0, fLoc0, dfLoc);
        interpolation = min(HInt, h);
    }
    else if ( disthole < d1 ){
        double dLoc1[2] = {0.0, d1}; // points conditions limites
        double fLoc1[2] = {h1, h}; // valeurs de f aux CL
        double hermInt = InterpolationHermite(disthole, dLoc1, fLoc1, dfLoc);
        interpolation = min(hermInt, h);
    }
    else{interpolation = h;}
    return interpolation;
    /*else{interpolation = h;}
    // Check if the point is inside the notch or hole
    if (dist_to_notch < r0) {
        return h0; // Size inside the notch
    } else if (dist_to_hole < r1) {
        return h1; // Size inside the hole
    } else {
        return h; // Default size outside the notch and hole
    }
    */
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(x0, y0, 0, w, h, -1, 0,&ierr);   
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0, r0, r0, -1, NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1, y1, 0, r1, r1, -1, NULL,0,NULL,0,&ierr);    
    ErrorGmsh(ierr);
    
    int plate[] = {2, idPlate};
    int notch[] = {2, idNotch};
    int hole[]  = {2, idHole};

    gmshModelOccCut(plate, sizeof(plate)/sizeof(plate[0]), notch, sizeof(notch)/sizeof(notch[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, sizeof(plate)/sizeof(plate[0]), hole, sizeof(hole)/sizeof(hole[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
    
    
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  //chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  //chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  //chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//    gmshFltkInitialize(&ierr);
//    gmshFltkRun(&ierr);  //chk(ierr);
//
    
}
