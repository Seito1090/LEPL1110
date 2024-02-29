#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = 0;
    
    //  A completer :-) done
    femDomain *theBoundary = malloc(sizeof(femDomain)); // Définition de theBoundary
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains, theGeometry->nDomains * sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains - 1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary * sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name, "Boundary");
    // Parcourir les bords pour compter les noeuds de frontière
    for (int i = 0; i < theEdges->nElem; ++i) {
        for (int j = 0; j < theEdges->nLocalNode; ++j) {
            int node = theEdges->elem[i*theEdges->nLocalNode + j];
            if (node == -1) continue;  // Noeud fictif pour les bords ouverts
            // Vérifier si le noeud est déjà présent dans le domaine de frontière
            int found = 0;
            for (int k = 0; k < nBoundary; ++k) {
                if (theBoundary->elem[k] == node) {
                    found = 1;
                    break;
                }
            }
            // Ajouter le noeud s'il n'est pas déjà présent
            if (!found) {
                theBoundary->elem[nBoundary] = node;
                nBoundary++;
            }
        }
    }

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");
 
    // A completer :-) done
    // Copier les noeuds de frontière dans theBoundary->elem
    int counter = 0;
    for (int i = 0; i < theEdges->nElem; ++i) {
        for (int j = 0; j < theEdges->nLocalNode; ++j) {
            int node = theEdges->elem[i * theEdges->nLocalNode + j];
            if (node == -1) continue;
            
            int found = 0;
            for (int k = 0; k < counter; ++k) {
                if (theBoundary->elem[k] == node) {
                    found = 1;
                    break;
                }
            }
            
            if (!found) {
                theBoundary->elem[counter++] = node;
            }
        }
    }

}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{

    // A completer :-) done
    // Libérer la mémoire allouée pour theProblem->geo
    geoMeshFree(theProblem->geo);

    // Libérer la mémoire allouée pour theProblem->space
    femDiscreteFree(theProblem->space);

    // Libérer la mémoire allouée pour theProblem->rule
    femIntegrationFree(theProblem->rule);

    // Libérer la mémoire allouée pour theProblem->system
    femFullSystemFree(theProblem->system);

    // Libérer la mémoire allouée pour theProblem lui-même
    free(theProblem);
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    
    //  A completer :-) done
    int nLocalNode = theMesh->nLocalNode;
    int *elemNodes = &(theMesh->elem[iElem*nLocalNode]);

    // Copier les indices de variables
    for (int i = 0; i < nLocalNode; ++i) {
        map[i] = elemNodes[i]; // Supposant que les indices de noeuds correspondent aux indices de variables
    }

    // Copier les coordonnées x et y
    for (int i = 0; i < nLocalNode; ++i) {
        int node = elemNodes[i];
        x[i] = theMesh->X[node];
        y[i] = theMesh->Y[node];
    }

}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
    int nLocal = theMesh->nLocalNode;

    // A completer :-) done
    // Boucle sur tous les éléments du maillage
    for (iElem = 0; iElem < theMesh->nElem; ++iElem) {
        // Récupérer les informations locales pour l'élément courant
        femPoissonLocal(theProblem, iElem, map, x, y);

        // Boucle sur tous les points d'intégration
        for (iInteg = 0; iInteg < theRule->n; ++iInteg) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            
            // Calculer les fonctions de forme et leurs dérivées
            femDiscretePhi(theSpace, xsi, eta, phi);
            femDiscreteDphi(theSpace, xsi, eta, dphidxsi, dphideta);

            // Calculer les dérivées par rapport à x et y
            for (i = 0; i < nLocal; ++i) {
                dphidx[i] = dphidxsi[i] * (1.0 / theMesh->jac[iElem]) + dphideta[i] * (1.0 / theMesh->jac[iElem+nLocal]);
                dphidy[i] = dphideta[i] * (1.0 / theMesh->jac[iElem]) + dphidxsi[i] * (1.0 / theMesh->jac[iElem+nLocal]);
            }

            // Assembler le système local
            double ke = 0.0; // Matrice de rigidité élémentaire
            double fe = 0.0; // Vecteur de force élémentaire

            for (i = 0; i < nLocal; ++i) {
                for (j = 0; j < nLocal; ++j) {
                    ke += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * theMesh->jac[iElem] * theRule->weight[iInteg];
                }
            }

            // Ajouter ke à la matrice de système global
            for (i = 0; i < nLocal; ++i) {
                for (j = 0; j < nLocal; ++j) {
                    int row = map[i];
                    int col = map[j];
                    femFullSystemAdd(theSystem, row, col, ke);
                }
            }

            // Ajouter fe au vecteur de force global
            // Dans ce problème, fe est nul, car la source est nulle
        }
    }

    // Résoudre le système linéaire global
    double *solution = femFullSystemEliminate(theSystem);
    // Stocker la solution dans theProblem->system->B
    for (i = 0; i < theSystem->size; ++i) {
        theSystem->B[i] = solution[i];
    }
    free(solution);
}

# endif



