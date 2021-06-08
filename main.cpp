#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include "defs.h"
#include <iostream>
#include <limits>
using namespace std;

char* inFilename;
char* outFilename;
uint32_t rootNumberToValidate;
int nIters;
// const   double DBL_MAX = std::numeric_limits<double>::max();

#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

void usage(int argc, char **argv)
{
    printf("Usage:\n");
    printf("    %s -in <input> [options]\n", argv[0]);
    printf("Options:\n");
    printf("    -in <input> -- input graph filename\n");
    printf("    -out <output> -- output filename. Default output is '<input>.v'\n");
    exit(1);
}

void init (int argc, char** argv, graph_t* G)
{
    inFilename = outFilename = NULL;
    nIters = -1;
    rootNumberToValidate = 0;
    for (int i = 1; i < argc; ++i) {
   		if (!strcmp(argv[i], "-in")) {
            inFilename = argv[++i];
        }
   		if (!strcmp(argv[i], "-out")) {
            outFilename = argv[++i];
        }
		if (!strcmp(argv[i], "-root")) {
			rootNumberToValidate = (int) atoi(argv[++i]);
        }
		if (!strcmp(argv[i], "-nIters")) {
			nIters = (int) atoi(argv[++i]);
        }
    }
    if (!inFilename) usage(argc, argv);
    if (!outFilename) {
        outFilename = (char *)malloc((strlen(inFilename) + 3) * sizeof(char));
        sprintf(outFilename, "%s.v", inFilename);
    }
}

void readGraph(graph_t *G, char *filename)
{
    edge_id_t arity;
    uint8_t align;
    FILE *F = fopen(filename, "rb");
    if (!F) error(EXIT_FAILURE, 0, "Error in opening file %s", filename);

    assert(fread(&G->n, sizeof(vertex_id_t), 1, F) == 1);
    G->scale = log(G->n) / log (2);   
    assert(fread(&arity, sizeof(edge_id_t), 1, F) == 1); 
    assert(fread(&G->directed, sizeof(bool), 1, F) == 1);
    assert(fread(&align, sizeof(uint8_t), 1, F) == 1);

    G->m = G->n * arity;

	G->rowsIndices = (edge_id_t *)malloc((G->n+1) * sizeof(edge_id_t));
    assert(G->rowsIndices);

	assert(fread(G->rowsIndices, sizeof(edge_id_t), G->n+1, F) == (G->n+1));
    
    G->endV = (vertex_id_t *)malloc(G->rowsIndices[G->n] * sizeof(vertex_id_t));
    assert(G->endV);

    assert(fread(G->endV, sizeof(vertex_id_t), G->rowsIndices[G->n], F) == G->rowsIndices[G->n]);

    assert(fread(&G->nRoots, sizeof(uint32_t), 1, F) == 1);
    
    G->roots = (vertex_id_t *)malloc(G->nRoots * sizeof(vertex_id_t));
    assert(G->roots);
    G->numTraversedEdges = (edge_id_t *)malloc(G->nRoots * sizeof(edge_id_t));
    assert(G->numTraversedEdges);

    assert(fread(G->roots, sizeof(vertex_id_t), G->nRoots, F) == G->nRoots);
    assert(fread(G->numTraversedEdges, sizeof(edge_id_t), G->nRoots, F) == G->nRoots);

    if (nIters == -1) nIters = G->nRoots;

    G->weights = (weight_t *)malloc(G->m * sizeof(weight_t));
    assert(G->weights);

    assert(fread(G->weights, sizeof(weight_t), G->m, F) == G->m);

    fclose(F);
}

/* write distances from root vertex to each others to output file. -1 = infinity */
void writeDistance(char* filename, weight_t *dist, vertex_id_t n)
{
    FILE *F = fopen(filename, "wb");
    assert(fwrite(dist, sizeof(weight_t), n, F) == n);
    fclose(F);
}
/* write number of vertices at each level */
void writeLevels (char* filename, vertex_id_t *validateNLevels, int validateNLevelsLength)
{
    FILE *F = fopen(filename, "w"); 
    fprintf(F, "%d\n", validateNLevelsLength);
    for (int i = 0; i < validateNLevelsLength; ++i) 
        fprintf(F, "%d\n", validateNLevels[i]);
}   




void freeGraph(graph_t *G)
{
    free(G->rowsIndices);
    free(G->endV);
    free(G->weights);
    free(G->roots);
    free(G->numTraversedEdges);
}

int main(int argc, char **argv) 
{
    weight_t* dist = NULL;
    graph_t g;
    struct timespec start_ts, finish_ts;
    double *perf;

    /* initializing and reading the graph */
    init(argc, argv, &g); 
    readGraph(&g, inFilename);
    
    init_sssp(&g);
    dist = (weight_t *)malloc(g.n * sizeof(weight_t));
    
    perf = (double *)malloc(g.nRoots * sizeof(double));
 
    nIters = (nIters < g.nRoots) ? nIters : g.nRoots; 
    // printf("start algorithm iterations...\n");
    rootNumberToValidate = 0;
    for (uint32_t i = 0; i < nIters; ++i) {
        /* initializing, DBL_MAX == infinity */

        for (vertex_id_t i = 0; i < g.n; i++) {
            dist[i] = DBL_MAX;
            // cout << i << ' ' << dist[i] << endl;
        }
       
        uint64_t nedges;
        clock_gettime(CLOCK, &start_ts);
        sssp(g.roots[i], &g, dist, &nedges);
        // cout << "endsssp" << endl;
        clock_gettime(CLOCK, &finish_ts);
        double time = (finish_ts.tv_nsec - (double)start_ts.tv_nsec) * 1.0e-9 + (finish_ts.tv_sec - (double)start_ts.tv_sec);
        
        perf[i] = time;
        printf("time = %.4f\n", time);

        if (rootNumberToValidate == i) {
            /* writing for validation */
            writeDistance(outFilename, dist, g.n);
        } 
    }
    // printf("algorithm iterations finished.\n");

    /* final print */
    double min_perf, max_perf, avg_perf;
    int nNonZeroRoots = 0;
    max_perf = avg_perf = 0;
    min_perf = DBL_MAX;
    for (uint32_t i = 0; i < nIters; ++i) {
        /* we consider search only with large number of traversed edges. Random root can point to connected component with very small number of edges */
        {
            ++nNonZeroRoots;
    	    avg_perf += perf[i];
            if (perf[i] < min_perf) min_perf = perf[i];
            if (perf[i] > max_perf) max_perf = perf[i];
        }
    }
    if (nNonZeroRoots != 0) avg_perf /= nNonZeroRoots;
    else error(EXIT_FAILURE, 0, "Number of roots with large traversed edges number is zero");

    // printf("%s: nIters = %d SSSP performance min = %.4f avg = %.4f max = %.4f \n", inFilename, nNonZeroRoots, min_perf, avg_perf, max_perf);

    freeGraph(&g);
    finalize_sssp();

    return 0;
}
