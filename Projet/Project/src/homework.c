#include "fem.h"
#include <time.h>


/* Partie juste pour le solveur bande */
# define MAX(a,b) (( a > b ) ?( a ) :( b ) )
# define MIN(a,b) (( a < b ) ?( a ) :( b ) ) 
double *theGlobalArrayPos;

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory) //DONE -> assiette 
//  (2) Ajouter les conditions de Neumann !   (mandatory) //DONE Neumann -> force Dirichlet -> déplacement
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised) //DONE
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory) //DONE -> bande

/* Fonctions donées de base dans le tempate du projet */


void femElasticityAssembleElements(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];  
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight;
      }
    }
  }
}

void femElasticityAssembleNeumann(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T){
      continue;
    }

    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
      }

      double tx = x[1] - x[0];
      double ty = y[1] - y[0];
      double length = hypot(tx, ty);
      double jac = length / 2.0;
      
      double f_x = 0.0;
      double f_y = 0.0;
      double f_n = 0.0;
      double f_t = 0.0;
      if (type == NEUMANN_X) {
        f_x = value;
      }
      if (type == NEUMANN_Y) {
        f_y = value;
      }
      
      double nx =  ty / length;
      double ny = -tx / length;
      double txn = -ny; // Tangent en x 
      double tyn = nx;  // Tangent en y 

      if (type == NEUMANN_N) {
        f_n = value;
      }
      if (type == NEUMANN_T) {
        f_t = value;
      }

      for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double weight = theRule->weight[iInteg];
        femDiscretePhi(theSpace, xsi, phi);
        for (i = 0; i < theSpace->n; i++) {
          B[2 * map[i] + 0] += jac * weight * phi[i] * (f_x + f_n * nx + f_t * txn);
          B[2 * map[i] + 1] += jac * weight * phi[i] * (f_y + f_n * ny + f_t * tyn);
        }
      }
    }
  }
}

void femElasticityApplyDirichlet(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double norm = sqrt(nx * nx + ny * ny);
      nx /= norm;
      ny /= norm;
      femFullSystemConstrain(theSystem, 2 * node + 0, value * nx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value * ny);
    }

    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double norm = sqrt(nx * nx + ny * ny);
      double tx = -ny / norm; 
      double ty = nx / norm;  
      femFullSystemConstrain(theSystem, 2 * node + 0, value * tx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value * ty);
    }

    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double norm = sqrt(nx * nx + ny * ny);
      nx /= norm;
      ny /= norm;
      double tx = -ny; 
      double ty = nx;  
      femFullSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * tx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);
    }

  }
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femFullSystemEliminate(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}

/* Solveur axisymetrique, pas utilisé dans notre cas étant donné la nature de notre problème */


void axisymetriqueAssembly(femProblem *theProblem){
  //L'idée ici, on calcule juste une section du problème et ensuite on le fait tourner :), le problème donné est toujours un problème 2D ! (P.S. : notre problème avec l'avion n'est pas du tout adapté pour 
  //l'axisymétrique, c'est pour ça que cette fonction n'est pas utilisée autrepart)

  //Dans le cas axisymétrique, on doit adapter les matrices A et B 
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight * x[i];
        B[mapY[i]] += phi[i] * gy * rho * jac * weight * x[i];
      }

      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * x[i] + dphidy[i] * c * dphidy[j] * x[i] + dphidx[i] * b * phi[j] + phi[i] * (b * dphidx[j] + a * (phi[j]/rho))) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * x[i] + dphidy[i] * c * dphidx[j] * x[i] + phi[i] * b * dphidy[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * x[i] + dphidx[i] * c * dphidy[j] * x[i] + dphidy[i] * b * phi[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * x[i] + dphidx[i] * c * dphidx[j] * x[i]) * jac * weight;
        }
      } 
      
    }
  }
}

double *femElasticitySolveAxisymetrique(femProblem *theProblem) {
  axisymetriqueAssembly(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = femFullSystemEliminate(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}

/* Solveur bande et ses autres fonctions, ceci est basé sur un des devoirs faits durant ce quadri */


// Modify compareNodePos to use theGlobalArrayPos, nécesssaire pour pouvoir renuméroter les noeuds
int compareNodePos(const void *nodeOne, const void *nodeTwo) {
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobalArrayPos[*iOne] - theGlobalArrayPos[*iTwo];
    return (diff < 0) ? -1 : (diff > 0) ? 1 : 0; // Return -1 / 0 / 1
}

//Fonction pour renuméroter les noeuds
void bandFemMeshRenumber(femMesh *theMesh, femRenumType renumType) {
    /* Idée ici c'est de simplement renuméroter les noeuds dans le sens inverse */
    int i;
    int *inverse = (int *)malloc(sizeof(int) * theMesh->nodes->nNodes);
    for (i = 0; i < theMesh->nodes->nNodes; i++)
        inverse[i] = i;

    switch (renumType) {
        case FEM_NO:
            break;
        case FEM_XNUM:
            theGlobalArrayPos = theMesh->nodes->X;
            qsort(inverse, theMesh->nodes->nNodes, sizeof(int), compareNodePos);
            break;
        case FEM_YNUM:
            theGlobalArrayPos = theMesh->nodes->Y;
            qsort(inverse, theMesh->nodes->nNodes, sizeof(int), compareNodePos);
            break;
        default:
            Error("Unexpected renumbering option");
    }

    for (i = 0; i < theMesh->nodes->nNodes; i++)
        theMesh->nodes->number[inverse[i]] = i;
    
    free(inverse);
}


//Fonctions pour assembler les éléments
void bandFemSystemAssemble(femFullSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc){
    for ( int i = 0; i < nLoc ; i ++) {
        int myRow = map [ i ];
        for ( int j = 0; j < nLoc ; j ++) {
            int myCol = map [ j ];

            // Ici , on ne regarde que les elements superieurs de la diagonale
            if ( myCol >= myRow ) {
                myBandSystem -> A [ myRow ][ myCol ] += Aloc [ i * nLoc + j ];
            }
        }
        myBandSystem -> B [ myRow ] += Bloc [ i ];
    }
}
void bandFemElasticityAssembleElements(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];  
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
          A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
          A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
          A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        B[mapX[i]] += phi[i] * gx * rho * jac * weight;
        B[mapY[i]] += phi[i] * gy * rho * jac * weight;
      }
    }
    bandFemSystemAssemble(theSystem, A[iElem], B, map, nLocal);
  }
}


//Fonction pour résoudre le système par élimination
double  *bandFemSystemEliminate(femBandSystem *myBand){
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    /* Gaussian elimination */
    for (k =0; k < size ; k ++) {
        if (fabs(A[k][k]) <= 1e-4) { Error("Cannot eliminate with such a pivot"); }

        // Limite de la ligne rouge
        jend = MIN ( k + band , size ) ;
        for (i = k+1 ; i < jend ; i ++) {

            // On recupere l"element symetrique
            factor = A [k][i] / A [k][k];
            for (j = i ; j < jend ; j ++){
                A [i][j] = A [i][j] - A [k][j] * factor ;
            }
            B[i] = B [i] - B [k] * factor ;
        }
    }
    /* Back - substitution */
    for (i = (size -1) ; i >= 0 ; i--){
        factor = 0;

        // Limite de la ligne rouge
        jend = MIN (i + band , size) ;
        for ( j = i +1 ; j < jend ; j ++){
            factor += A [i][j] * B [j];
        }
        B[i] = (B [i] - factor) / A [i][i];
    }
    return (myBand -> B) ;
}

//Fonction solve appelée par le main
double *bandFemElasticitySolve(femProblem *theProblem) {
  bandFemMeshRenumber(theProblem->geometry->theElements, FEM_XNUM);
  bandFemElasticityAssembleElements(theProblem); 
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = bandFemSystemEliminate((femBandSystem *)theProblem->system); 
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}

/* CG + CSR */
# define TOL 10e-6
# define MAX_ITER 1000
typedef struct sparseMatrix {
  int size;
  int nnz;
  int *col;
  int *rptr;
  double *val;
} sparseMatrix;
void sparseMatrixFree(sparseMatrix *sp) {
  free(sp->col);
  free(sp->rptr);
  free(sp->val);
  free(sp);
}

sparseMatrix* to_sparse(double **A, int size) {
  int nnz = 0;
  for (int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      nnz += (A[i][j] != 0);
    }
  }
  int* col = malloc(nnz * sizeof(int));
  int* rptr = malloc((size + 1) * sizeof(int));
  double* val = malloc(nnz * sizeof(double));
  
  nnz = 0;
  for(int i = 0; i < size; i++) {
    rptr[i] = nnz;
    for(int j = 0; j < size; j++) {
      if(A[i][j] != 0) {
        col[nnz] = j;
        val[nnz] = A[i][j];
        nnz++;
      }
    }
  }
  rptr[size] = nnz;

  sparseMatrix* sp = malloc(sizeof(sparseMatrix));
  sp->size = size;
  sp->nnz = nnz;
  sp->col = col;
  sp->rptr = rptr;
  sp->val = val;
  return sp;
}

static inline void spmv(const sparseMatrix* sp, const double* x, double* y) {
  for(int i = 0; i < sp->size; i++) {
    double s = 0;
    for(int j = sp->rptr[i]; j < sp->rptr[i+1]; j++) {
      s += sp->val[j] * x[sp->col[j]];
    }
    y[i] = s;
  }
}

static inline void residual(const sparseMatrix* sp, const double* x, const double* b, double* r) {
  spmv(sp, x, r);
  for(int i = 0; i < sp->size; i++) {
    r[i] = b[i] - r[i];
  }
}

static inline double dot(const double* x, const double* y, int size) {
  double s = 0;
  for(int i = 0; i < size; i++) {
    s += x[i] * y[i];
  }
  return s;
}

static inline void axpy(double* x, const double* y, double a, int size) {
  for(int i = 0; i < size; i++) {
    x[i] += a * y[i];
  }
}

double *solve_cg(femFullSystem *mySystem) {
  int size = mySystem->size;
  double **A = mySystem->A;
  double *B = mySystem->B;
  //clock_t start, stop;

  //start = clock();
  sparseMatrix* sp = to_sparse(A, size);
  //stop = clock();
  //printf("Time to create sparse mat: %f ms\n", 1000 * (double)(stop - start) / CLOCKS_PER_SEC);

  int niter = 0;
  double *x = malloc(size * sizeof(double));
  double *r = malloc(size * sizeof(double));
  double *p = malloc(size * sizeof(double));
  double *Ap = malloc(size * sizeof(double));

  // Initialize x, r, and p
  for (int i = 0; i < size; i++) {
    x[i] = 0.0;
    r[i] = B[i];
    p[i] = r[i];
  }

  double alpha, beta, rr, rrNew;
  rr = dot(r, r, size);

  while (niter < size) {
    spmv(sp, p, Ap);
    alpha = rr / dot(p, Ap, size);
    axpy(x, p, alpha, size);
    axpy(r, Ap, -alpha, size);
    rrNew = dot(r, r, size);
    beta = rrNew / rr;
    for (int i = 0; i < size; i++) {
      p[i] = r[i] + beta * p[i];
    }
    rr = rrNew;
    niter++;
  }

  free(r);
  free(p);
  free(Ap);

  return x;
}

double *CGfemElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  double *soluce = solve_cg(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}



