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
#include <set>
#include <map>
#include <vector>
#include <unordered_map>
#include "shmem.h"


using namespace std;

char* inFilename;
char* outFilename;
uint32_t rootNumberToValidate;
int nIters;
int lgsize;
const double delta = 0.01;

const int INF_BUCKET = 10000000;
set<int> curbuckets;
map<int, set<vertex_id_t> > vertex_in_bucket;
int bufsize = 1000000;


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

void readGraph(graph_t *G, char *filename, vector<int> &offsets)
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
    int n = G->n;
    G->m = G->n * arity;
    int m = G->m;
    edge_id_t* rowsIndices = new edge_id_t[n + 1];
    assert(fread(rowsIndices, sizeof(edge_id_t), n + 1, F) == n + 1);
    vertex_id_t* endV = new vertex_id_t[rowsIndices[n]];
    assert(fread(endV, sizeof(vertex_id_t), rowsIndices[n], F) == rowsIndices[n]);

    int rank = shmem_my_pe();
    int size = shmem_n_pes();
    G->nproc = size;
    G->rank = rank;
    for (lgsize = 0; lgsize < size; ++lgsize) {
        if ((1 << lgsize) == size) {
            break;
        }
    }
    G->local_n = n / size;
    G->rowsIndices = new edge_id_t [G->local_n + 1];
    G->rowsIndices[0] = 0;
    vertex_id_t idx = 1;
    for (vertex_id_t i = 0; i < n; ++i) {
        if (VERTEX_OWNER2(i) == rank) {
            G->rowsIndices[idx] = G->rowsIndices[idx - 1] + rowsIndices[i + 1] - rowsIndices[i];
            ++idx;
        }
        offsets[VERTEX_OWNER2(i)]++;
    }
    G->endV = new vertex_id_t [G->rowsIndices[G->local_n]];

    assert(fread(&G->nRoots, sizeof(uint32_t), 1, F) == 1);
    
    G->roots = (vertex_id_t *)malloc(G->nRoots * sizeof(vertex_id_t));
    assert(G->roots);
    G->numTraversedEdges = (edge_id_t *)malloc(G->nRoots * sizeof(edge_id_t));
    assert(G->numTraversedEdges);

    assert(fread(G->roots, sizeof(vertex_id_t), G->nRoots, F) == G->nRoots);
    assert(fread(G->numTraversedEdges, sizeof(edge_id_t), G->nRoots, F) == G->nRoots);

    if (nIters == -1) nIters = G->nRoots;

    weight_t* weights = new weight_t [m];
    assert(fread(weights, sizeof(weight_t), m, F) == G->m);
    G->weights = new weight_t [G->rowsIndices[G->local_n]];
    idx = 0;
    for (vertex_id_t i = 0; i < n; ++i) {
        if (VERTEX_OWNER2(i) == rank) {
            for (edge_id_t j = rowsIndices[i]; j < rowsIndices[i + 1]; ++j) {
                G->endV[idx] = endV[j];
                G->weights[idx] = weights[j];
                ++idx;
            }
        }
    }

    delete[] weights;
    delete[] endV;
    delete[] rowsIndices;
    fclose(F);
}




/* write distances from root vertex to each others to output file. -1 = infinity */
void writeDistance(char* filename, weight_t *dist, vertex_id_t n, vector<int> & offsets)
{
    if (!my_pe()) {
        FILE *F = fopen(filename, "wb");
        vector<double> alldists(n);
        for(int i = 0; i < shmem_n_pes(); ++i) {
            if (i != 0) {
                double tempdist[offsets[i]];
                shmem_double_get(tempdist, (const double *) dist, offsets[i], i);
                for(int j = 0; j < offsets[i]; ++j) {
                    alldists[VERTEX_TO_GLOBAL(i, j)] = tempdist[j];
                }
            } else {
                for(int j = 0; j < offsets[i]; ++j) {
                    alldists[VERTEX_TO_GLOBAL(i, j)] = dist[j];
                }
            }
        }
        assert(fwrite(alldists.data(), sizeof(weight_t), n, F) == n);
        fclose(F);
    } 
    shmem_barrier_all();
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

void move_to_new_bucket(int bucket, vertex_id_t v, double dist) {
    vertex_in_bucket[bucket].erase(v);
    vertex_in_bucket[(int) (dist / delta)].insert(v);
    if (bucket != INF_BUCKET && vertex_in_bucket[bucket].size() == 0) {
        vertex_in_bucket.erase(bucket);
    }
    curbuckets.insert((int) (dist / delta));
}

void shmem_sssp(vertex_id_t root, graph_t *G, weight_t *dist) {
    struct timespec start_ts, finish_ts;
    int *numdist[shmem_n_pes()];
    double *changeddist[shmem_n_pes()];
    int *changedindex[shmem_n_pes()];
    unsigned long long numsends = 0;
    for(int i = 0 ; i < shmem_n_pes(); ++i) {
        numdist[i] = (int *) shmalloc(sizeof(int));
        numdist[i][0] = 0;
        changeddist[i] = (double *) shmalloc (sizeof(double) * bufsize);
        changedindex[i] = (int *) shmalloc(sizeof(int) * bufsize);
    }
    int *mymem = (int*) shmalloc(sizeof(int));
    int *othermem = (int*) shmalloc(sizeof(int));
    long *psync = (long *) shmalloc(SHMEM_REDUCE_SYNC_SIZE * sizeof(long));
    int *pwrk = (int *) shmalloc(SHMEM_REDUCE_MIN_WRKDATA_SIZE * sizeof(int));
    if (VERTEX_OWNER2(root) == my_pe()) {
        dist[VERTEX_LOCAL2(root)] = 0;
        curbuckets.insert(0);
        vertex_in_bucket[INF_BUCKET].erase(VERTEX_LOCAL2(root));
        vertex_in_bucket[0].insert(VERTEX_LOCAL2(root));
    }
    *mymem = *(curbuckets.begin());
    *othermem = -1;
    for(int i= 0 ; i < SHMEM_REDUCE_SYNC_SIZE; ++i) {
        psync[i] = _SHMEM_SYNC_VALUE;
    }
    shmem_barrier_all();
    shmem_int_min_to_all(othermem, mymem, 1, 0, 0, shmem_n_pes(), pwrk, psync);
    shmem_flush();
    shmem_barrier_all();
    shmem_free(pwrk);
    shmem_free(psync);
    int numbucket = *othermem;
    while(numbucket != INF_BUCKET) {
        vector<vertex_id_t> my_upd;
        vector<int> not_my_upd1[shmem_n_pes()];
        vector<double> not_my_upd2[shmem_n_pes()];
        vector<int> numbuckets;
        if (vertex_in_bucket.find(numbucket) != vertex_in_bucket.end()) {
            for(std::set<vertex_id_t>::iterator from = vertex_in_bucket[numbucket].begin(); from != vertex_in_bucket[numbucket].end(); ++from) {
                for(edge_id_t idx = G->rowsIndices[*from]; idx < G->rowsIndices[*from + 1]; ++idx) {
                    vertex_id_t to = G->endV[idx];
                    const double weight = G->weights[idx];
                    if (VERTEX_OWNER2(to) == my_pe()) {
                        to = VERTEX_LOCAL2(to);
                        if (dist[to] > dist[*from] + weight) {
                            if (abs(dist[to] - DBL_MAX) < 0.00001) {
                                numbuckets.push_back(INF_BUCKET);
                            } else {
                                numbuckets.push_back((int) (dist[to] / delta));
                            }
                            dist[to] = dist[*from] + weight;
                            my_upd.push_back(to);
                        }
                    // } else {
                    //     int owner = VERTEX_OWNER2(to);
                    //     double distto;
                    //     shmem_double_get(&distto, (const double *) &(dist[VERTEX_LOCAL2(to)]), 1, owner);
                    //     shmem_flush();
                    //     numsends++;
                    //     if (distto > dist[*from] + weight) {
                    //         distto = dist[*from] + weight;
                    //         not_my_upd1[owner].push_back(VERTEX_LOCAL2(to));
                    //         not_my_upd2[owner].push_back(distto);
                    //     }
                    // }
                    } else {
                        int owner = VERTEX_OWNER2(to);
                        not_my_upd1[owner].push_back(VERTEX_LOCAL2(to));
                        not_my_upd2[owner].push_back(dist[*from] + weight);
                    }
                }
            }
        }
        
        for(int i = 0 ; i < shmem_n_pes(); ++i) {
            int num = not_my_upd1[i].size();
            if (i != my_pe() && num != 0) {
                shmem_int_put(numdist[my_pe()], &num, 1, i);
                shmem_flush();
                shmem_int_put(changedindex[my_pe()], not_my_upd1[i].data(), num, i);
                shmem_flush();
                shmem_double_put(changeddist[my_pe()], not_my_upd2[i].data(), num, i);
                shmem_flush();
            }
        }
        shmem_barrier_all();
        unordered_map<int, double> all_changes;
        for(int i = 0 ; i < shmem_n_pes(); ++i) {
            if (i != my_pe()) {
                for(int j = 0 ; j < numdist[i][0]; ++j) {
                    if(all_changes.find(changedindex[i][j]) == all_changes.end()) {
                        all_changes[changedindex[i][j]] = changeddist[i][j]; 
                    } else {
                        all_changes[changedindex[i][j]] = min(all_changes[changedindex[i][j]], changeddist[i][j]);
                    }
                }
            }
        }
        for(unordered_map<int, double>::iterator i = all_changes.begin(); i != all_changes.end(); ++i ) {
            if (dist[i->first] > i->second) {
                if (abs(dist[i->first] - DBL_MAX) < 0.00001) {
                    numbuckets.push_back(INF_BUCKET);
                } else {
                    numbuckets.push_back((int) (dist[i->first] / delta));
                }
                dist[i->first] = i->second;
                my_upd.push_back(i->first);
            }
        }
        curbuckets.erase(numbucket);
        int k = 0;
        for(vector<vertex_id_t>::iterator i = my_upd.begin(); i != my_upd.end(); ++i) {
            move_to_new_bucket(numbuckets[k], *i, dist[*i]);
            ++k;
        }

        *mymem = *(curbuckets.begin());
        *othermem = INF_BUCKET;
        psync = (long *) shmalloc(SHMEM_REDUCE_SYNC_SIZE * sizeof(long));
        pwrk = (int *) shmalloc(SHMEM_REDUCE_MIN_WRKDATA_SIZE * sizeof(int));
        for(int i= 0 ; i < SHMEM_REDUCE_SYNC_SIZE; ++i) {
            psync[i] = _SHMEM_SYNC_VALUE;
        }
        shmem_barrier_all();
        shmem_int_min_to_all(othermem, mymem, 1, 0, 0, shmem_n_pes(), pwrk, psync);
        shmem_flush();
        shmem_free(pwrk);
        shmem_free(psync);
        shmem_barrier_all();
        numbucket = *othermem;
        // if (my_pe() == 0) {
        //     cout << "numbucket = " << numbucket << endl;
        // }
    }
    shmem_free(mymem);
    shmem_free(othermem);
    for(int i = 0 ; i < shmem_n_pes(); ++i) {
        shmem_free(numdist[i]);
        shmem_free(changeddist[i]);
        shmem_free(changedindex[i]);
    }

    if (my_pe() == 0) {
        cout << "numsends = " << numsends << endl;
    }

}


int main(int argc, char **argv) 
{
    shmem_init();
    weight_t* dist = NULL;
    graph_t g;
    struct timespec start_ts, finish_ts;
    int rootNumberToValidate = 0;

    /* initializing and reading the graph */
    init(argc, argv, &g); 
    vector<int> offsets(shmem_n_pes());
    for(int i = 0 ; i < offsets.size(); ++i) {
        offsets[i] = 0;
    }
    readGraph(&g, inFilename, offsets);
    
    // init_sssp(&g);
    dist = (weight_t *) shmalloc(g.local_n * sizeof(weight_t));
    
 
    double alltime = 0;
    for (uint32_t i = 0; i < 1; ++i) {
        /* initializing, DBL_MAX == infinity */
        curbuckets.clear();
        curbuckets.insert(INF_BUCKET);
        vertex_in_bucket.clear();

        for (vertex_id_t j = 0; j < g.local_n; j++) {
            dist[j] = DBL_MAX;
            vertex_in_bucket[INF_BUCKET].insert(j);
        }
       
        clock_gettime(CLOCK, &start_ts);
        
        shmem_sssp(g.roots[i], &g, dist);
        shmem_barrier_all();
        clock_gettime(CLOCK, &finish_ts);
        double time = (finish_ts.tv_nsec - (double)start_ts.tv_nsec) * 1.0e-9 + (finish_ts.tv_sec - (double)start_ts.tv_sec);
        cout << my_pe() << " time = " << time << endl;
        if (rootNumberToValidate == i) {
            /* writing for validation */
            writeDistance(outFilename, dist, g.n, offsets);
        } 
        alltime += time;
    }
    if (my_pe() == 0) {
        cout << "alltime = " << alltime << endl;
    }

    freeGraph(&g);
    finalize_sssp();
    shmem_free(dist);
    shmem_finalize();
    return 0;
}
