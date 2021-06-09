typedef struct graph{
  int size;

  double beta;
  double delta;
  double tau;
  double lambda;
  
  int infected;
  int susceptible;
  double avginfected; 

  int **adjmat;
  int *nodevec;
  int *infvec;
}graph;

void createGraph(graph* g, int size, double p_m0, double beta, double delta);
void setInfVec( graph* g );
void setStatus( graph* g, int node, int status);
void addEdge( graph* g, int nodeA, int nodeB);
void removeEdge( graph* g, int nodeA, int nodeB);
void removeNode( graph* g, int node );
void printGraph( graph* g);
void deleteGraph( graph* g);
void simEpidemic( graph* g, int burnin, int timesteps );
