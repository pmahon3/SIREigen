#include "Graph.h"
#include "MT19937.h"
#include "AuxillaryFunctions.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define NTHREADS 4

void *mainfun(void *x){

  // Get thread id
  int tid = *((int *) x);

  // Create a results file for the thread

  char buf[0x100];
  snprintf(buf, sizeof(buf), "Thread%d.csv", tid); 

  
  FILE* results  = fopen(buf, "w");

  sgenrand(time(0));

  int tsteps = 900;
  int burnin = 100;
  int size = 1000;

  int nreplicates = 3;

  // Graph initialization parameters; uncomment either 1. or 2. for the requisite graph and...
  // select the corresponding adjacency matrix generating function in the createGraph function in Graph.c.
  // 1. Edge probability for an Erdos-Renyi graph
  double p = 0.01;
  // 2. Mean degree for Barabasi-Albert graph
  //double m_0 = 2;

  // Graph initialization
  graph* g = malloc(sizeof(struct graph));
  createGraph(g, size, p, 0, 0);

  // Run simulations

  double beta = 0.05;
  double delta = 0.9;
  double score = beta / delta * g->lambda;

  while( delta >= 0.01){

    g->beta = beta;
    g->delta = delta;

    double totalavginfected = 0;
    for ( int i = 0; i < nreplicates ; i++ ){
      simEpidemic(g, burnin, tsteps);

      totalavginfected += g->avginfected;
    }

    printf("Thread %d\nLambda: %f\nBeta: %f\nDelta: %f\nScore: %f\nRatio infected: %f\n\n", tid, g->lambda, g->beta, g->delta, score, totalavginfected / (double) nreplicates );
    fprintf(results, "%f, %f,\n", score, totalavginfected / (double) nreplicates );
    delta -= 0.01;
    score = beta / delta * g->lambda;
  }

  fclose(results);

  printf("Thread %d complete.\n", tid);

  deleteGraph(g);
  return NULL;
}

int main(){

  // Parallelization code sourced from gribblelab.org/CBootCamp/A2_Parallel_Programming_in_C.html  
  pthread_t threads[NTHREADS];
  int thread_args[NTHREADS];
  int rc, i;

  /* spawn the threads */
  for (i=0; i<NTHREADS; ++i)
    {
      thread_args[i] = i;
      printf("spawning thread %d\n", i);
      rc = pthread_create(&threads[i], NULL, mainfun, (void *) &thread_args[i]);
    }

  /* wait for threads to finish */
  for (i=0; i<NTHREADS; ++i) {
    rc = pthread_join(threads[i], NULL);
  }

  return 0;
}
