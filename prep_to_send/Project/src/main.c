/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

int main(void) {
  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../src/data/mesh.txt");
  femProblem *theProblem = femElasticityRead(theGeometry, "../src/data/problem.txt");
  femElasticityPrint(theProblem);
  double *theSoluce = CGfemElasticitySolve(theProblem);
  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, "../src/data/UV.txt");
  femElasticityFree(theProblem);
  geoFree();
  return 0;
}
