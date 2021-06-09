#include "Graph.h"
#include "AuxillaryFunctions.h"
#include "MT19937.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void simEpidemic(graph* g, int burnin, int timesteps ){

  
  // Select initial infected population
  g->infected = 0;
  g->susceptible = g->size;
  
  for ( int i = 0 ; i < (int)g->size ; i++) {
    g->nodevec[i] = 1;
    g->infected += 1;
    g->susceptible -= 1; 
  }
 
  setInfVec(g);

  // Cumalitve total infected over the course of the simulation
  int cuminf = 0;
  
  // Run time steps
  for ( int time = 0 ; time < timesteps ; time++ ){
   
    // Set vector tracking infected negihbours
    setInfVec( g );

    // newnodevec stores the changes to the current state.
    // g->nodevec stores the current state (i..e. what we are using to compute newnodevec.
    int* newnodevec = malloc(sizeof(int)*g->size);
    
    // Evaluate each node
    for ( int node = 0 ; node < g->size ; node++ ){
      
      newnodevec[node] = g->nodevec[node];

      // Infection
      if ( g->nodevec[node] == 0 ){
	
	for ( int i = 0 ; i < g->infvec[node] ; i++ ){
	  if ( compareDouble( g->beta, genrand()) ){
	    //printf("Node Infected: %d\tProb: %f\n\n", node, nodephi);
	    newnodevec[node] = 1;
	    g->infected += 1;
	    g->susceptible -=1;
	    break;
	  }
	}
      }
      // Recovery
      else if ( g->nodevec[node] == 1 ){
	if ( compareDouble( g->delta, genrand() ) ){
	  //printf("Node Recovered: %d\tProb: %f\n\n", node, prob);
	  newnodevec[node] = 0;
	  g->susceptible += 1;
	  g->infected -= 1;
	}
      }
    }

    // Replace current state with updated state.
    for ( int i = 0; i < g->size ; i++ ){
      g->nodevec[i] = newnodevec[i];
    }
    free(newnodevec);

    // Add to cumalitive infected if burn in is finished
    
    if( time > burnin ) cuminf += g->infected;
  }

  // Set average infected
  g->avginfected = (double) cuminf / (double) (timesteps - burnin) / (double) g->size;
}



void createGraph(graph* g, int size, double p_m0, double beta, double delta){ 

  g->size = size;
  g->beta = beta;
  g->delta = delta;
  g->infected = 0;
  g->susceptible = size;
  g->avginfected = 0;

  g->adjmat = malloc(sizeof(int*)*size);
  for ( int i = 0 ; i < size ; i++ ){
    g->adjmat[i] = malloc(sizeof(int)*size);
  }

  // Adjacency matrix generation; either uncomment 1. or 2. for the requisite graph type
  // 1. Generate adjacency matrix for Erdos-Renyi graph
  genERAdjMat( g->adjmat, g->size, p_m0 );

  // 2. Generate adjacency matrix for Barabasi-Albert power law graph
  // genBAAdjMat( g->adjmat, g->size, (int) p_m0 )
  
  g->nodevec = malloc(sizeof(int)*size);
  g->infvec = malloc(sizeof(int)*size);
  g->lambda = powerMethod(g->adjmat, 1000, g->size);
  g->tau = 1/g->lambda;
}

void addEdge( graph* g, int nodeA, int nodeB){

  g->adjmat[nodeA][nodeB] = 1;
  g->adjmat[nodeB][nodeA] = 1;

  return;
}

void removeEdge( graph* g, int nodeA, int nodeB){
  g->adjmat[nodeA][nodeB] = 0;
  g->adjmat[nodeB][nodeA] = 0;
  return;
}

void setInfVec( graph* g ){
  int size = g->size;
  
  for ( int i = 0; i < size; i++ ){
    int sum = 0;
    for ( int j = 0; j < size; j++){
      sum += g->adjmat[i][j] * ( g->nodevec[j] % 2);
    }
    g->infvec[i] = sum;
  }
}

void setStatus( graph* g, int node, int status ){
  g->nodevec[node] = status;
}

void removeNode( graph* g, int node ){
  int size = g->size;

  for ( int i = 0 ; i < size ; i++ ){
    g->adjmat[node][i] = 0;
    g->adjmat[i][node] = 0;
  }
}

void deleteGraph(graph* g){
  free(g->nodevec);
  free(g->infvec);
  for( int i = 0 ; i < g->size; i++ ){
    free(g->adjmat[i]);
  }
  free(g->adjmat);
  free(g);
}

void printGraph( graph* g){

  int size = g->size;
  printf("Status code:\n");
  printf("\t0 = Susceptible\n\t1 = Infected\n\t2 = immune\n\n");
  printf("======================\n");
  printf("|Node|Status|Infected Contacts|Contacts\n");
  printf("========================================\n");
  
  for ( int i = 0 ; i < size; i++){

    printf("|%4d|%6d|%17d|", i, g->nodevec[i], g->infvec[i]);
    for ( int j = 0; j < size; j++ ){
      if ( g->adjmat[i][j] == 1) printf(" %2d ", j);
    }
    printf("\n");
    printf("---------------------------------------\n");
  }
}
