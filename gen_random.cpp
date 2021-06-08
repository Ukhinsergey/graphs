#include "defs.h"
#include <stdio.h>
#include <cstdlib>
#include <assert.h>
#include <string.h>
#include <error.h>
#include <vector>
#include <random>
#include <iostream>

using namespace std;

char outFilename[FNAME_LEN];

/* helper */
void usage(int argc, char **argv)
{
    printf("Random graph generator\n");
    printf("Usage:\n");
    printf("%s -s <scale> [other options]\n", argv[0]);
    printf("Options:\n");
    printf("   -s <scale>, number of vertices is 2^<scale>\n");
    printf("   -k <half the average vertex degree>, default value is 16\n");
    printf("   -out <output filename>, file for the graph storage\n");
    exit(1);
}

/* initialization */
void init(int argc, char **argv, graph_t *G)
{
	bool no_out_filename = true;
    G->scale = -1;
    G->permute_vertices = true;
    G->min_weight = 0;
    G->max_weight = 1;
    /* default value */
    G->avg_vertex_degree = DEFAULT_ARITY;
    if (argc == 1) {
        usage(argc, argv);
    }
    
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-s")) {
            G->scale = (int)atoi(argv[++i]);
        }
        
        if (!strcmp(argv[i], "-k")) {
            G->avg_vertex_degree = (int)atoi(argv[++i]);
        }
        
		if (!strcmp(argv[i], "-out")) {
            no_out_filename = false;
			sprintf(outFilename, argv[++i]);
        }
    }
    
	if (no_out_filename) {
    	sprintf(outFilename, "random-%d", G->scale);
	}
	
    if (G->scale == -1) {
        usage(argc, argv);
    }
    
    G->n = (vertex_id_t)1 << G->scale;
    G->m = G->n * G->avg_vertex_degree;
}

/* random graph generator */
void gen_random_graph(graph_t *G)
{
    /* init */
    vertex_id_t n;
    edge_id_t m;
    edge_id_t offset;
    bool permute_vertices;
    vertex_id_t *permV, tmpVal;
    vertex_id_t u, v;
    vertex_id_t *src;
    vertex_id_t *dest;
    unsigned *degree;
    permute_vertices = G->permute_vertices;
    double *dbl_weight;
    double min_weight, max_weight;
    n = G->n;
    m = G->m;
    src = new vertex_id_t[m];
    assert(src != NULL);
    dest = new vertex_id_t[m];
    assert(dest != NULL);
    degree = new unsigned[n];
    assert(degree != NULL);
    memset(degree, 0, sizeof(unsigned) * n);
    
    dbl_weight = (double *) malloc(m * sizeof(double));
    assert(dbl_weight != NULL);

    srand48(2387);
    
    /* generate edges */
    for (edge_id_t i = 0; i < m; i++) {
        vertex_id_t u = rand() % n;
        vertex_id_t v = rand() % n;
        src[i] = u;
        dest[i] = v;
    }
    
    /* reshuffle */
    if (permute_vertices) {
        srand48(4791);
        permV = new vertex_id_t[n];
        assert(permV != NULL);
        
        for (vertex_id_t i = 0; i < n; i++) {
            permV[i] = i;
        }
        
        for (vertex_id_t i = 0; i < n; i++) {
            vertex_id_t j = n * drand48();
            tmpVal = permV[i];
            permV[i] = permV[j];
            permV[j] = tmpVal;
        }
        
        for (edge_id_t i = 0; i < m; i++) {
            src[i] = permV[src[i]];
            dest[i] = permV[dest[i]];
        }
        
        delete[] permV;
    }
    
    min_weight = G->min_weight;
    max_weight = G->max_weight;

    /* Generate edge weights */
    for (edge_id_t i=0; i<m; i++) {
        dbl_weight[i]  = min_weight + (max_weight-min_weight)*drand48();
    }

    /* update graph data structure */
    for (edge_id_t i = 0; i < m; i++) {
        degree[src[i]]++;
        degree[dest[i]]++;
    }
    
    G->endV = new vertex_id_t[2 * m];
    assert(G->endV != NULL);
    
    G->rowsIndices = new edge_id_t[n + 1];
    assert(G->rowsIndices != NULL);
    
    G->n = n;
    /* undirected graph, each edge is stored twice; if edge is (u, v), then it's
     * stored at the vertex u and at the vertex v */
    G->m = 2 * m;
    
    G->weights = (double *) malloc(G->m * sizeof(double));       
    assert(G->weights != NULL);

    G->rowsIndices[0] = 0;
    for (vertex_id_t i = 1; i <= G->n; i++) {
        G->rowsIndices[i] = G->rowsIndices[i - 1] + degree[i - 1];
    }
    
    for (edge_id_t i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        offset = degree[u]--;
        G->endV[G->rowsIndices[u] + offset - 1] = v;
        G->weights[G->rowsIndices[u]+offset-1] = dbl_weight[i];
        offset = degree[v]--;
        G->endV[G->rowsIndices[v] + offset - 1] = u;
        G->weights[G->rowsIndices[v]+offset-1] = dbl_weight[i];

    }
    G->nRoots = 10;
    G->roots = (vertex_id_t *)malloc(G->nRoots * sizeof(vertex_id_t));
    for (int i = 0; i < G->nRoots; ++i) {
        G->roots[i] = static_cast<vertex_id_t>( random() % n);
    }
    delete[] src;
    delete[] dest;
    delete[] degree;
}

/* write graph to file */
void writeGraph(graph_t *G, char *filename)
{
    FILE *F = fopen(filename, "wb");
    if (!F) error(EXIT_FAILURE, 0, "Error in opening file %s", filename);
	
    assert(fwrite(&G->n, sizeof(vertex_id_t), 1, F) == 1);
    
    edge_id_t arity = G->m / G->n;
    assert(fwrite(&arity, sizeof(edge_id_t), 1, F) == 1);
    assert(fwrite(&G->directed, sizeof(bool), 1, F) == 1);
    uint8_t align = 0;
    assert(fwrite(&align, sizeof(uint8_t), 1, F) == 1);

    assert(fwrite(G->rowsIndices, sizeof(edge_id_t), G->n+1, F) == G->n+1);
    assert(fwrite(G->endV, sizeof(vertex_id_t), G->rowsIndices[G->n], F) == G->rowsIndices[G->n]);
    
    assert(fwrite(&G->nRoots, sizeof(uint32_t), 1, F) == 1);
    assert(fwrite(G->roots, sizeof(vertex_id_t), G->nRoots, F) == G->nRoots);
    assert(fwrite(&G->numTraversedEdges, sizeof(edge_id_t), G->nRoots, F) == G->nRoots);

    assert(fwrite(G->weights, sizeof(weight_t), G->m, F) == G->m);
    fclose(F);
}

/* print graph */
void printGraph(graph_t *G)
{
	int i,j;
	for (i = 0; i < (int)G->n; ++i) {
		printf("%d:", i);
		for (j=G->rowsIndices[i]; j < (int)G->rowsIndices[i+1]; ++j)
			printf("%d (%f), ", G->endV[j], G->weights[j]);
		printf("\n");
	}
}

int main(int argc, char **argv) {
    graph_t g;
    init(argc, argv, &g);
    gen_random_graph(&g);
    // printGraph(&g);
    writeGraph(&g, outFilename);
    
    return 0;
}
