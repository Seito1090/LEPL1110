/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"

int main(void) {

  //
  //  -1- Construction de la geometrie
  //
  geoInitialize();
  femGeo *theGeometry = geoGetGeometry();

  //Lecture d'un fichier .msh
  geoMeshGenerateMshFile("../../data/avion.msh");


  geoMeshImport();
  geoSetDomainName(0, "underRW");
  geoSetDomainName(1, "RWL");
  geoSetDomainName(2, "upperRW");
  geoSetDomainName(3, "underLW");
  geoSetDomainName(4, "LWL");
  geoSetDomainName(5, "upperLW");
  geoSetDomainName(6, "VStabL");
  geoSetDomainName(7, "VStabR");
  geoSetDomainName(8, "HStabR");
  geoSetDomainName(9, "HStabL");
  geoSetDomainName(10, "Fuselage0");
  geoSetDomainName(11, "Wheels");
  geoSetDomainName(12, "Fuselage1");
  geoSetDomainName(13, "Fuselage2");
  geoSetDomainName(14, "Fuselage3");
  geoSetDomainName(15, "Fuselage4");
  geoSetDomainName(16, "Fuselage5");
  geoSetDomainName(17, "Cab");
  geoMeshWrite("../../data/mesh.txt");

  //
  //  -2- Definition du probleme
  //

  double E = 294e9;
  double nu = 0.3;
  double rho = 7.85e3;
  double gx = 0;
  double gy = 9.81;
  double lift = -3.5e6;
  double mass = 398e3;
  double weight = mass * gy;

  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, PLANAR_STRAIN);
  femElasticityAddBoundaryCondition(theProblem, "underRW", NEUMANN_Y, lift, NAN); 
  femElasticityAddBoundaryCondition(theProblem, "underRW", DIRICHLET_X, 0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "upperRW", NEUMANN_Y, weight, NAN); 
  femElasticityAddBoundaryCondition(theProblem, "underLW", NEUMANN_Y, lift, NAN);
  femElasticityAddBoundaryCondition(theProblem, "underLW", DIRICHLET_X, 0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "upperLW", NEUMANN_Y, weight, NAN);
  femElasticityAddBoundaryCondition(theProblem, "Wheels", DIRICHLET_Y, -0.014, NAN);
  femElasticityAddBoundaryCondition(theProblem, "HStabR", NEUMANN_Y, -lift/33.5, NAN);
  femElasticityAddBoundaryCondition(theProblem, "HStabL", NEUMANN_Y, -lift/33.5, NAN);
  femElasticityAddBoundaryCondition(theProblem, "VStabR", DIRICHLET_X, 0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "VStabL", DIRICHLET_X, 0, NAN);

  femElasticityPrint(theProblem);
  femElasticityWrite(theProblem, "../../data/problem.txt");

  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
