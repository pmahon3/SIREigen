#include "AuxillaryFunctions.h"
#include "MT19937.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// powerMethod computes the principle eigenvalue of a given matrix using the power method or power iteration. 
double powerMethod( int** mat, int iters, int size){
  double* x_k = malloc(sizeof(double) * size );
  double* y_k = malloc(sizeof(double) * size );
  
  for ( int i = 0 ; i < size ; i++ ){
    y_k[i] = 1 + genrand();
  }

  double m_k = 0.0;

  for ( int j = 0 ; j < iters; j++){
    m_k = 0;
    for ( int i = 0 ; i < size ; i++ ){
      double sum = 0;
      for ( int j = 0 ; j < size ; j++ ){
	sum += mat[i][j] * y_k[j];
      }
      x_k[i] = sum;
      if ( compareDouble(sum, m_k) ) m_k = sum;
    }

    for ( int i = 0 ; i < size ; i++ ) y_k[i] = x_k[i] / m_k;
  }

  free(y_k);
  free(x_k);
  return m_k;
}

// compareDouble returns 1 if a is bigger than b and zero otherwise.
int compareDouble( double a, double b){
  double result = a - b;
  int sgn = copysign(1.0, result);

  if ( (int) sgn > 0 ){
    return 1;
  }
  else return 0;
}

// genERAdjMat generates an Erdos-Renyi graph where each edge has a probability, p, of existing.

void genERAdjMat(int** adjmat, int size, double p ){
  for ( int i = 0 ; i < size; i++ ){
    for ( int j = i+1 ; j < size ; j++) {
	if ( compareDouble( p, genrand() ) ){
	  adjmat[i][j] = 1;
	  adjmat[j][i] = 1;
	}
    }
  }
}

// genBAAdjMat generates an adjacency matrix with mean node degree m0, with node degrees distributed according to the probability P(X=x) = x^(-3)
// For algorithm details see:
//      Barabasi and Albert. (1999) Emergence of scaling in random networs. Science 286, 509-512
//      Albert and Barabasi. (2002) Statistical mechanics of complex networks. Reviews of Modern Physics 74, 48-85
//      Prettejohn et al. (2011) Methods for generating complex networks with selected structural properties for simulations: a review and tutorial for neuroscientists. Fronties in Computational Neuroscience

void genBAAdjMat(int** adjmat, int size, int m0 ){

  // Initialize the matrix
  int i,j;

  // Intialize vector to store node degrees
  int* degmat = malloc(sizeof(int)*size);
  
  for ( i = 0; i < size; i++ ){
    degmat[i] = 0;
  }

  //Set initial connected nodes
  int totaldeg = 0;
  for ( i = 0 ; i < m0 ; i++ ){
    for ( j = i + 1; j < m0; j++ ){
      adjmat[i][j] = 1;
      adjmat[j][i] = 1;
      degmat[i]++;
      degmat[j]++;
      totaldeg++;
    }
  }
  
  // Add edges to remaining nodes with preference given to nodes with higher degree
  // For all yet to be connected nodes...
  for ( i = m0 ; i < size ; i++ ){

    // while the current nodes degree is less than the seed...
    while( degmat[i] < m0 ) {

      // select a random node j from the set of nodes already connected, not equal to i and not adjacent to i,...
      j = (int)(i * genrand());
            
      while ( adjmat[i][j] == 1 || j == i ){
	j = i * genrand();
      }

      // compute the probability of adding an between i and j 
      double prob = (double)degmat[j]/(double)totaldeg;
 
      // add the edge if a random number on [0,1] is less than the probability.
      if ( genrand() < prob ){
	adjmat[i][j] = 1;
	adjmat[j][i] = 1;
	degmat[i]++;
	degmat[j]++;
	totaldeg++;
      }
    }
  }
  free(degmat);
  }
