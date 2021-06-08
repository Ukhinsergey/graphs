#include <iostream>
#include <vector>
#include <algorithm>
#include "defs.h"

using namespace std;

extern "C" void init_sssp(graph_t *G)
{
}

extern "C" void finalize_sssp()
{
}

/* Reference SSSP implementation --- Dijkstra algorithm 
 * it is needed to carefully count the number of traversed edges
 * dist must be initialized as -1
*/
extern "C" void sssp(vertex_id_t root, graph_t *G, weight_t *dist, uint64_t *traversed_edges)
{
    /* (distance,node) */
    // cout << __FILE__ <<endl;
    vector< pair<weight_t, vertex_id_t> > pq(G->n);
    // cout << "vector" << endl;
    uint64_t nedges = 0;
    pq.at(0) = make_pair(0, root); 
    // cout << "make_pair" << endl;
    int len = 1; 
    // cout << "before dist " << root<< endl;
    dist[root] = 0;
    // cout << "after dist" << endl;
    // cout << len << endl;
    while ( len ) {
        pair<weight_t, vertex_id_t> s = pq[0];
        pop_heap(pq.begin(), pq.begin()+len, greater<pair<weight_t, vertex_id_t> >()); 
        len--;
        if ( (dist[s.second] >= s.first) || (dist[s.second] == -1) ) {
            for ( unsigned int i = G->rowsIndices[s.second]; i < G->rowsIndices[s.second+1]; ++i ) {
                if ( (dist[G->endV[i]] > dist[s.second] + G->weights[i]) || (dist[G->endV[i]] == -1) ) {
                    dist[G->endV[i]] = dist[s.second] + G->weights[i];
                    pq.at(len) = make_pair(dist[G->endV[i]], G->endV[i]); 
                    len++;
                    if ((unsigned)len+1 > pq.size()) { pq.resize(2*len); }
                    push_heap(pq.begin(), pq.begin()+len, greater<pair<weight_t, vertex_id_t> >());
                }
	            /* we count every traversed edge */
	            ++nedges;
            }
        }
    }
    *traversed_edges = nedges;
}

