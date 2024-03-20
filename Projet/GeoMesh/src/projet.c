#include "fem.h"
#include <math.h>

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
 
    int ierr;
    double r_outer = 1.0; // Outer radius of the donut
    double r_inner = 0.8; // Inner radius of the donut
    double wing_length = 5.0; // Length of each wing
    double wing_width = 0.2; // Width of each wing

    // Add the outer circle
    int idFuselage = gmshModelOccAddDisk(x0, y0, 0, r_outer, r_outer, -1, NULL, 0, NULL, 0, &ierr);
    ErrorGmsh(ierr);

    // Add the inner circle (the hole)
    int idInner = gmshModelOccAddDisk(x0, y0, 0, r_inner, r_inner, -1, NULL, 0, NULL, 0, &ierr);
    ErrorGmsh(ierr);

    // Calculate coordinates for wing attachment points
    double x_left_wing = x0 - r_outer - wing_length + 0.2;
    double x_right_wing = x0 + r_outer - 0.2;
    double y_wings = y0 - 0.4;

    int idLeftWing = gmshModelOccAddRectangle(x_left_wing, y_wings - wing_width / 2, 0, wing_length, wing_width, -1, 0, &ierr);
    ErrorGmsh(ierr);

    int idRightWing = gmshModelOccAddRectangle(x_right_wing, y_wings - wing_width / 2, 0, wing_length, wing_width, -1, 0, &ierr);
    ErrorGmsh(ierr);

    int fuselage[] = {2, idFuselage};
    int inner[] = {2, idInner};
    int leftWing[] = {2, idLeftWing};
    int rightWing[] = {2, idRightWing};
    
    // Cut the outer circle with the inner circle to create the hole
    gmshModelOccCut(fuselage, sizeof(fuselage)/sizeof(fuselage[0]), inner, sizeof(inner)/sizeof(inner[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    // Attach wings to the fuselage
    gmshModelOccFuse(fuselage, sizeof(fuselage)/sizeof(fuselage[0]), leftWing, sizeof(leftWing)/sizeof(leftWing[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccFuse(fuselage, sizeof(fuselage)/sizeof(fuselage[0]), rightWing, sizeof(rightWing)/sizeof(rightWing[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    // Calculate coordinates for wingtip attachment points
    double x_left_wingtip = x_left_wing - 0.2 - r_outer;
    double x_right_wingtip = x_right_wing + wing_length;

    // Define pivot point for rotation
    double pivot_x = x_right_wing + wing_length; // Pivot point x-coordinate
    double pivot_y = y_wings; // Pivot point y-coordinate

    // Define rotation angle in radians (positive for counterclockwise rotation)
    double rotation_angle = 65 * 3.14 / 180.0; // Rotate by 40 degrees

    // Add wingtip to the right wing
    int idRightWingtip = gmshModelOccAddRectangle(x_right_wingtip, y_wings - wing_width / 2, 0, 1.2, wing_width+0.02, -1, 0, &ierr);
    ErrorGmsh(ierr);

    // Add wingtip to the left wing
    int idLeftWingtip = gmshModelOccAddRectangle(x_left_wingtip, y_wings - wing_width / 2, 0, 1.2, wing_width+0.02, -1, 0, &ierr);
    ErrorGmsh(ierr);

    int rightWingtip[] = {2, idRightWingtip};
    int leftWingtip[] = {2, idLeftWingtip};

    // Rotate the wingtip
    gmshModelOccRotate(rightWingtip, sizeof(rightWingtip)/sizeof(rightWingtip[0]), pivot_x, pivot_y, 0, cos(rotation_angle), sin(rotation_angle), 0, rotation_angle, &ierr);
    ErrorGmsh(ierr);

    // Attach rotated wingtip to the right wing
    gmshModelOccFuse(rightWing, sizeof(rightWing)/sizeof(rightWing[0]), rightWingtip, sizeof(rightWingtip)/sizeof(rightWingtip[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    

    // Define pivot point for rotation of left wingtip
    pivot_x = x_left_wingtip + 0.2 + r_outer; // Update pivot point for left wingtip

    // Rotate the left wingtip
    gmshModelOccRotate(leftWingtip, sizeof(leftWingtip)/sizeof(leftWingtip[0]), pivot_x, pivot_y, 0, cos(rotation_angle), sin(-rotation_angle), 0, -rotation_angle, &ierr);
    ErrorGmsh(ierr);

    // Attach rotated wingtip to the left wing
    gmshModelOccFuse(leftWing, sizeof(leftWing)/sizeof(leftWing[0]), leftWingtip, sizeof(leftWingtip)/sizeof(leftWingtip[0]), NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
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
