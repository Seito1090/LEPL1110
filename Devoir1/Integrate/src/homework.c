#include <stdio.h>
#include <math.h>
#include "glfem.h"

double changeOfVariable(double v[3], double xi, double eta);

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    //plop the information given in the pdf
    double weights[3] = { 1/6.0, 1/6.0, 1/6.0 };
    double eta[3] = { 1/6.0, 2/3.0, 1/6.0 };
    double xi[3] = { 1/6.0, 1/6.0, 2/3.0 };

    //Now we calculate the local coordinates
    for (int a = 0; a < 3; a++){
        double xiLocal = xi[a];
        double etaLocal = eta[a];
        xLoc[a] = changeOfVariable(x, xiLocal, etaLocal);
        yLoc[a] = changeOfVariable(y, xiLocal, etaLocal);
        I += weights[a] * f(xLoc[a], yLoc[a]);
    }

    //Calculate the Jacobian cuz the teach said so during class
    double J = fabs((x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]));

//
// ... A modifier :-)
//
//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
//

  glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
  glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
  glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    return I * J;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    //Faut trouver un moyen de split le triangle en 4 mais la flemme de faire tout en hard coding, y doit y avoir un moyen d'automatiser ca ...
    //Cas de base 
    if (n == 0){
        return integrate(x, y, f);
    }

    double splitsX[4][3] = {{},{},{},{}}; //we can store the 4 triangles in there in the future
    double splitsY[4][3] = {{},{},{},{}}; //same here
    //What about the weights and the xi and eta ? :(
    //Idea : well since the weights xi and eta are taken care of in the integrate function, we don't really have to worry about them ?
    //First triangle
    splitsX[0][0] = x[0];
    splitsX[0][1] = (x[0]+x[1])/2;
    splitsX[0][2] = (x[0]+x[2])/2;
    splitsY[0][0] = y[0];
    splitsY[0][1] = (y[0]+y[1])/2;
    splitsY[0][2] = (y[0]+y[2])/2;

    //Second triangle
    splitsX[1][0] = x[1];
    splitsX[1][1] = (x[0]+x[1])/2;
    splitsX[1][2] = (x[1]+x[2])/2;
    splitsY[1][0] = y[1];
    splitsY[1][1] = (y[0]+y[1])/2;
    splitsY[1][2] = (y[1]+y[2])/2;

    //Third triangle
    splitsX[2][0] = x[2];
    splitsX[2][1] = (x[0]+x[2])/2;
    splitsX[2][2] = (x[1]+x[2])/2;
    splitsY[2][0] = y[2];
    splitsY[2][1] = (y[0]+y[2])/2;
    splitsY[2][2] = (y[1]+y[2])/2;

    //Fourth triangle
    splitsX[3][0] = (x[0]+x[1])/2;
    splitsX[3][1] = (x[0]+x[2])/2;
    splitsX[3][2] = (x[1]+x[2])/2;
    splitsY[3][0] = (y[0]+y[1])/2;
    splitsY[3][1] = (y[0]+y[2])/2;
    splitsY[3][2] = (y[1]+y[2])/2;


    //Now we calculate the integral of each triangle and sum them up
    double I = 0;

    for (int i = 0; i < 4; i++){
        I += integrateRecursive(splitsX[i], splitsY[i], f, n-1);
    }
//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    
//
//
//    
     
    return I;
}

double changeOfVariable(double v[3], double xi, double eta){ //It's here cuz this thing will be used a lot, easier on the eyes when reading it afterwards 
    return v[0] * xi + v[1] * eta + v[2] * (1 - xi - eta) ;
}