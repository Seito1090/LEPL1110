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

  // OPTION 1 : Utilisation de GMSH avec OpenCascade
  // theGeometry->h = 0.05;
  // geoMeshGenerate();

  // OPTION 2 : Utilisation de GMSH directement
  // theGeometry->h = 0.05;
  // geoMeshGenerateGeo();

  // OPTION 3 : Lecture d'un fichier .geo
  //theGeometry->h = 0.05;
  //geoMeshGenerateGeoFile("../data/mesh.geo");

  // OPTION 4 : Lecture d'un fichier .msh
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
  //femElasticityAddBoundaryCondition(theProblem, "Wheels", DIRICHLET_XY, 0, 0);
  femElasticityAddBoundaryCondition(theProblem, "Wheels", DIRICHLET_XY, 0, 0);
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
  //Bloque dirichlet pour le bas de la cellule 
  femElasticityPrint(theProblem);
  femElasticityWrite(theProblem, "../data/problem.txt");

  //
  //  -3- Champ de la taille de référence du maillage (uniquement pour la visualisation)
  //

  double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  femNodes *theNodes = theGeometry->theNodes;
  for (int i = 0; i < theNodes->nNodes; ++i)
    meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
  double hMin = femMin(meshSizeField, theNodes->nNodes);
  double hMax = femMax(meshSizeField, theNodes->nNodes);
  printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
  printf(" ==== Minimum h          : %14.7e \n", hMin);
  printf(" ==== Maximum h          : %14.7e \n", hMax);

  //
  //  -4- Visualisation
  //
/*
  int mode = 1;
  int domain = 0;
  int freezingButton = FALSE;
  double t, told = 0;
  char theMessage[MAXNAME];

  GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
  glfwMakeContextCurrent(window);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  do {
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    glfemReshapeWindows(window, theGeometry->theNodes, w, h);

    t = glfwGetTime();
    if (glfwGetKey(window, 'D') == GLFW_PRESS) {
      mode = 0;
    }
    if (glfwGetKey(window, 'V') == GLFW_PRESS) {
      mode = 1;
    }
    if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) {
      domain++;
      freezingButton = TRUE;
      told = t;
    }

    if (t - told > 0.5) {
      freezingButton = FALSE;
    }
    if (mode == 1) {
      glfemPlotField(theGeometry->theElements, meshSizeField);
      glfemPlotMesh(theGeometry->theElements);
      sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
      glColor3f(1.0, 0.0, 0.0);
      glfemMessage(theMessage);
    }
    if (mode == 0) {
      domain = domain % theGeometry->nDomains;
      glfemPlotDomain(theGeometry->theDomains[domain]);
      sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
      glColor3f(1.0, 0.0, 0.0);
      glfemMessage(theMessage);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
  } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
  */
  // Check if the ESC key was pressed or the window was closed

  free(meshSizeField);
  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
