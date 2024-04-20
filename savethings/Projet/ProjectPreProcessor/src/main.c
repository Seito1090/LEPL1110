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
  geoMeshGenerateMshFile("../data/avion.msh");


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
  geoMeshWrite("../data/mesh.txt");

  //
  //  -2- Definition du probleme
  //

  double E = 211.e9;
  double nu = 0.3;
  double rho = 7.85e3;
  double gx = 0;
  double gy = -9.81;
  double lift = 300e4;
  double mass = 200e3;
  double weight = mass * gy;

  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, PLANAR_STRAIN);
  femElasticityAddBoundaryCondition(theProblem, "underRW", NEUMANN_Y, lift/2, NAN); 
  femElasticityAddBoundaryCondition(theProblem, "underRW", DIRICHLET_X, 0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "upperRW", NEUMANN_Y, weight/3, NAN); 
  femElasticityAddBoundaryCondition(theProblem, "underLW", NEUMANN_Y, lift/2, NAN);
  femElasticityAddBoundaryCondition(theProblem, "underLW", DIRICHLET_X, 0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "upperLW", NEUMANN_Y, weight/3, NAN);
  femElasticityAddBoundaryCondition(theProblem, "Wheels", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Wheels", NEUMANN_Y, weight/3, NAN);
  femElasticityAddBoundaryCondition(theProblem, "Cab", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Fuselage0", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Fuselage1", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Fuselage2", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Fuselage3", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Fuselage4", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Fuselage5", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "HStabL", DIRICHLET_NT, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "HStabR", DIRICHLET_NT, 0, 0);
  femElasticityPrint(theProblem);
  femElasticityWrite(theProblem, "../data/problem.txt");

  free(meshSizeField);
  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
