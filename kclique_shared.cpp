/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.

To compile:
"gcc DDegColNodeParallel.c -O9 -o DDegColNodeParallel -fopenmp".

To execute:
"./DDegColNodeParallel p k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
k is the size of the k-cliques
p is the number of threads
Will print the number of k-cliques.
*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <unordered_map>

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed


#define KAR_DEBUG
#undef KAR_DEBUG

unsigned sharedVar;

typedef struct edge{
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
				   //unsigned *map;//oldID newID correspondance NOT USED IN THIS VERSION
    unsigned *sgExist; //does the edge exist in a particular subgraph
} edgelist;

typedef struct {
	unsigned n;
    unsigned e;
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors
    std::pair<unsigned, unsigned> *adjEid;
    //std::vector<std::pair<unsigned, unsigned>> adjEid;
    unsigned *eid;
	unsigned core;//core value of the graph
} graph;

typedef struct {
	unsigned *n;//n[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *adj;//truncated list of neighbors
	unsigned char *lab;//lab[i] label of node i
	unsigned **nodes;//sub[l]: nodes in G_l
	unsigned core;
	unsigned *color;
} subgraph;

typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;

typedef struct {
	unsigned id;
	unsigned color;
} idcolor;

int cmp(const void* a, const void* b)
{
	iddegree *x = (iddegree*)a, *y = (iddegree*)b;
	return y->degree - x->degree;
}

int cmpadj(const void* a, const void* b)
{
	idcolor *x = (idcolor*)a, *y = (idcolor*)b;
	return y->color - x->color;
}

void free_edgelist(edgelist *el) {
	free(el->edges);
	free(el->rank);
	free(el);
}

void free_graph(graph *g) {
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph *sg, unsigned char k) {
	unsigned char i;
	free(sg->n);
	for (i = 2; i < k; i++) {
		free(sg->d[i]);
		free(sg->nodes[i]);
	}
	free(sg->d);
	free(sg->nodes);
	free(sg->lab);
	free(sg->adj);
	free(sg);
}


//Compute the maximum of three unsigned integers.
unsigned int max3(unsigned int a, unsigned int b, unsigned int c);
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

edgelist* readedgelist(char* input) {
	unsigned e1 = NLINKS;
	edgelist *el = (edgelist*) malloc(sizeof(edgelist));
	FILE *file;

	el->n = 0;
	el->e = 0;
    unsigned s,t;
	file = fopen(input, "r");
	el->edges = (edge*) malloc(e1 * sizeof(edge));
	while (fscanf(file, "%u %u", &s, &t) == 2) {//Add one edge
        if (s==t) continue;
        el->edges[el->e].s = s; el->edges[el->e].t = t;
		el->n = max3(el->n, el->edges[el->e].s, el->edges[el->e].t);
		el->e++;
		if (el->e == e1) {
			e1 += NLINKS;
			el->edges = (edge*) realloc(el->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	el->n++;

	el->edges = (edge*) realloc(el->edges, el->e * sizeof(edge));

	return el;
}

void relabel(edgelist *el) {
	unsigned i, source, target, tmp;

	for (i = 0; i < el->e; i++) {
		source = el->rank[el->edges[i].s];
		target = el->rank[el->edges[i].t];
		if (source < target) {
			tmp = source;
			source = target;
			target = tmp;
		}
		el->edges[i].s = source;
		el->edges[i].t = target;
	}

}

///// CORE ordering /////////////////////
typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;

bheap *construct(unsigned n_max) {
	unsigned i;
	bheap *heap = (bheap*) malloc(sizeof(bheap));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned*) malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = (keyvalue*) malloc(n_max * sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap, unsigned i, unsigned j) {
	keyvalue kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

void bubble_up(bheap *heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n) {
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap, keyvalue kv) {
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_up(heap, heap->n - 1);
}

void update(bheap *heap, unsigned key) {
	unsigned i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value)--;
		bubble_up(heap, i);
	}
}

keyvalue popmin(bheap *heap) {
	keyvalue min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n, unsigned *v) {
	unsigned i;
	keyvalue kv;
	bheap* heap = construct(n);
	for (i = 0; i < n; i++) {
		kv.key = i;
		kv.value = v[i];
		insert(heap, kv);
	}
	return heap;
}

void freeheap(bheap *heap) {
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(edgelist* el) {
	unsigned i, j, r = 0, n = el->n, e = el->e;
	keyvalue kv;
	bheap *heap;

	unsigned *d0 = (unsigned*) calloc(el->n, sizeof(unsigned));
	unsigned *cd0 = (unsigned*) malloc((el->n + 1) * sizeof(unsigned));
	unsigned *adj0 = (unsigned*) malloc(2 * el->e * sizeof(unsigned));
	for (i = 0; i < e; i++) {
		d0[el->edges[i].s]++;
		d0[el->edges[i].t]++;
	}
	cd0[0] = 0;
	for (i = 1; i < n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];
		d0[i - 1] = 0;
	}
	for (i = 0; i < e; i++) {
		adj0[cd0[el->edges[i].s] + d0[el->edges[i].s]++] = el->edges[i].t;
		adj0[cd0[el->edges[i].t] + d0[el->edges[i].t]++] = el->edges[i].s;
	}

	heap = mkheap(n, d0);

	el->rank = (unsigned*) malloc(n * sizeof(unsigned));
	for (i = 0; i < n; i++) {
		kv = popmin(heap);
		el->rank[kv.key] = n - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
			update(heap, adj0[j]);
		}
	}

	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}


//////////////////////////
//Building the special graph
graph* mkgraph(edgelist *el) {
	unsigned i, max;
	unsigned *d;
    unsigned *din;
	graph* g = (graph*) malloc(sizeof(graph));

	d = (unsigned*) calloc(el->n, sizeof(unsigned));
    

	for (i = 0; i < el->e; i++) {
		d[el->edges[i].s]++;
	}

	g->cd = (unsigned*) malloc((el->n + 1) * sizeof(unsigned));
	g->cd[0] = 0;
	max = 0;
	for (i = 1; i < el->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		d[i - 1] = 0;
	}

	//printf("core value (max truncated degree) = %u\n", max);

	g->adj = (unsigned*) malloc(el->e * sizeof(unsigned));

	for (i = 0; i < el->e; i++) {
		g->adj[g->cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
	}

	free(d);
	g->core = max;
	g->n = el->n;
	return g;
}



/*****************************
    MY CODE HERE
   
Extract smaller subgraphs (cyclic node IDs)
Filter based on trussness

*****************************/


unsigned updateGlobalQueue (unsigned locQWrPtr, unsigned locQSize, unsigned* globQWrPtr, unsigned* locBuff, unsigned* globBuff)
{
    if (locQWrPtr >= locQSize)
    {
        unsigned tempIdx = __sync_fetch_and_add(globQWrPtr, locQSize);
        for (unsigned bufIdx = 0; bufIdx < locQSize; bufIdx++)
            globBuff[tempIdx + bufIdx] = locBuff[bufIdx];
        locQWrPtr = 0;
    }
    return locQWrPtr;
}


/*******************************************************************************************************/
void trussScan(unsigned numEdges, int *EdgeSupport, unsigned level, unsigned *curr, unsigned *currTail, bool *InCurr) {
    // Size of cache line
    const unsigned BUFFER_SIZE_BYTES = 2048;
    const unsigned BUFFER_SIZE = BUFFER_SIZE_BYTES/sizeof(unsigned);

    #pragma omp single
    sharedVar = 0;
    #pragma omp barrier

    unsigned buff[BUFFER_SIZE];
    unsigned index = 0;


    #pragma omp for schedule(static) 
    for(long i = 0; i < numEdges; i++) {
	    if( EdgeSupport[i] <= level ) {
            buff[index++] = i;
            InCurr[i] = true;

            index = updateGlobalQueue(index, BUFFER_SIZE, currTail, buff, curr);
	    }
    }

    if(index > 0) {
	    unsigned tempIdx = __sync_fetch_and_add(currTail, index);

	    for(unsigned j = 0; j < index; j++) {
	        curr[tempIdx+j] = buff[j];
	    }
    }


#pragma omp barrier

}


//Process a sublevel in a level using intersection based approach
//void PKT_processSubLevel_intersection(graph *g, unsigned *curr, bool *InCurr, unsigned currTail, int *EdgeSupport, int level, unsigned *next, unsigned *nextTail, bool *processed, std::vector<edge>& edgeIdtoEdge) {
void PKT_processSubLevel_intersection(graph *g, unsigned *curr, bool *InCurr, unsigned currTail, int *EdgeSupport, int level, unsigned *next, unsigned *nextTail, bool *processed, edge* edgeIdtoEdge) {

    //Size of cache line
    const unsigned BUFFER_SIZE_BYTES = 2048;
    const unsigned BUFFER_SIZE = BUFFER_SIZE_BYTES/sizeof(unsigned);

    unsigned buff[BUFFER_SIZE];
    unsigned index = 0;

    #pragma omp single
    sharedVar=0;
    #pragma omp barrier

#pragma omp for schedule(dynamic,4)
    for (unsigned i = 0; i < currTail; i++) {

	    //process edge <u,v>
        unsigned e1 = curr[i]; 

	    edge ev = edgeIdtoEdge[e1];  

	    unsigned u = ev.s;
	    unsigned v = ev.t;

#ifdef KAR_DEBUG
        assert(e1 < g->e/2);
        assert(EdgeSupport[e1] <= level);
        assert(u < g->n); assert(v < g->n);
#endif

	    unsigned uStart = g->cd[u], uEnd = g->cd[u+1];
        unsigned vStart = g->cd[v], vEnd = g->cd[v+1];

        unsigned int numElements = (uEnd - uStart) + (vEnd - vStart);
        unsigned j_index = uStart, k_index = vStart;

	    for(unsigned int innerIdx = 0; innerIdx < numElements; innerIdx ++) {
	        if( j_index >= uEnd) {
		        break;
	        }
            else if( k_index >= vEnd ) {
		        break;
            }
            else if( g->adj[j_index] == g->adj[k_index] ) {
		        unsigned e2 = g->eid[ k_index ];  //<v,w>
                unsigned e3 = g->eid[ j_index ];  //<u,w>

#ifdef KAR_DEBUG
                assert(e2 < g->e/2); assert(e3 < g->e/2);
#endif


                //If e1, e2, e3 forms a triangle
		        if( (!processed[e2]) && (!processed[e3]) ) {

		            //Decrease support of both e2 and e3
		            if( EdgeSupport[e2] > level && EdgeSupport[e3] > level) {

			            //Process e2
			            int supE2 = __sync_fetch_and_sub( &EdgeSupport[e2], 1);
			            if( supE2 == (level+1) ) 
			                buff[index++] = e2;
			            if( supE2 <= level ) 
			                __sync_fetch_and_add(&EdgeSupport[e2],1);

                        index = updateGlobalQueue(index, BUFFER_SIZE, nextTail, buff, next);

			            //Process e3
			            int supE3 = __sync_fetch_and_sub(&EdgeSupport[e3], 1);

			            if( supE3 == (level +1) ) 
			                buff[index++] = e3;

			            if(supE3 <= level ) {
			                __sync_fetch_and_add(&EdgeSupport[e3],1);
			            }

                        index = updateGlobalQueue(index, BUFFER_SIZE, nextTail, buff, next);
		            }
		            else if(EdgeSupport[e2] > level ) {
			            //process e2 only if e1 < e3
			            if((e1 < e3) || (!InCurr[e3])) {
			                int supE2 = __sync_fetch_and_sub(&EdgeSupport[e2], 1);
                            if (supE2 == (level+1)) buff[index++] = e2;
                            if (supE2 <= level) __sync_fetch_and_add(&EdgeSupport[e2], 1);
                            index = updateGlobalQueue(index, BUFFER_SIZE, nextTail, buff, next);
			                	
			            }
		            }
		            else if(EdgeSupport[e3] > level ) {
			            //process e3 only if e1 < e2
			            if((e1 < e2) || (!InCurr[e2])) {
			                int supE3 = __sync_fetch_and_sub(&EdgeSupport[e3], 1);
                            if (supE3 == (level+1)) buff[index++] = e3;
                            if (supE3 <= level) __sync_fetch_and_add(&EdgeSupport[e3], 1);
                            index = updateGlobalQueue(index, BUFFER_SIZE, nextTail, buff, next);
			                	
			            }
		            }

		       }	

               j_index ++;
               k_index ++;
            }
            else if( g->adj[j_index] < g->adj[k_index] ) {
                j_index++;
#ifdef KAR_DEBUG
                if (j_index != uEnd)
                    assert(g->adj[j_index] > g->adj[j_index-1]);
#endif
            }
            else if( g->adj[k_index] < g->adj[j_index] ) {
                k_index++;
#ifdef KAR_DEBUG
                if (k_index != vEnd)
                    assert(g->adj[k_index] > g->adj[k_index-1]);
#endif
            }
	    }
    }

    if (index > 0) {
        long tempIdx =  __sync_fetch_and_add(nextTail, index);;
        for (long bufIdx = 0; bufIdx < index; bufIdx++)
            next [tempIdx + bufIdx] = buff[bufIdx];
    }

#pragma omp barrier

#pragma omp for schedule(static)
    for (long i = 0; i < currTail; i++) {
        unsigned e = curr[i];  
	    processed[e] = true;
	    InCurr[e] = false;
    }

#pragma omp barrier

}



void serialPrefix(unsigned* in, unsigned* out, unsigned start, unsigned end)
{
    assert(start!=0);
    out[start] = in[start-1]; 
    for (unsigned i = start+1; i <= end; i++)
        out[i] = out[i-1]+in[i-1]; 
}

void thrdPrefix(unsigned* arr, unsigned BS, unsigned NTHRD, unsigned len)
{
    arr[0] = 0;
    for (unsigned i = 0; i < NTHRD; i++)
        arr[std::min((i+1)*BS, len)] += arr[i*BS]; 
}

void applyThrdPrefix(unsigned* arr, unsigned start, unsigned end)
{
    unsigned prevSum = arr[start-1];
    for (unsigned i = start; i < end; i++)
        arr[i] += prevSum;
}

void triangleCount(graph* g, unsigned* ds, int* supp, unsigned maxOutDeg)
{ 
    long long unsigned int totCount = 0;
    #pragma omp parallel
    {
        thread_local std::unordered_map<unsigned, unsigned> neighSet (2*maxOutDeg);
        #pragma omp for reduction(+:totCount)
        for (unsigned i = 0; i < g->n; i++)
        {
            //hash neighbors
            for (unsigned j = g->cd[i]; j<g->cd[i]+ds[i]; j++)
                neighSet.insert(std::make_pair(g->adj[j], g->eid[j]));

            //join with neighbors of neighbors
            for (unsigned j = g->cd[i]; j<g->cd[i]+ds[i]; j++)
            {
                unsigned neigh = g->adj[j];
                unsigned e1 = g->eid[j];
                for (unsigned k = g->cd[neigh]; k < g->cd[neigh]+ds[neigh]; k++)
                {
                    unsigned neighOfNeigh = g->adj[k];
                    std::unordered_map<unsigned, unsigned>::const_iterator got = neighSet.find(neighOfNeigh);
                    if (got == neighSet.end())
                        continue;
                    unsigned e2 = g->eid[k];
                    unsigned e3 = got->second;
#ifdef KAR_DEBUG
                    assert(e1!=e2); assert(e2!=e3); assert(e1!=e3);
                    assert(e1 <= g->e/2); assert(e2 <= g->e/2); assert(e3 <= g->e/2);
#endif
                    __sync_fetch_and_add(&supp[e1], 1);
                    __sync_fetch_and_add(&supp[e2], 1);
                    __sync_fetch_and_add(&supp[e3], 1);
                    totCount++;
                } 
            }

            for (unsigned j = g->cd[i]; j<g->cd[i] + ds[i]; j++)
                neighSet.erase(g->adj[j]);
        }
    }
    //printf("total count = %llu\n", totCount);
}

//overloaded function, discards deleted edges
void triangleCount(graph* g, unsigned* ds, bool* deleted, int* supp, unsigned maxOutDeg)
{ 
    #pragma omp for
    for (unsigned i = 0; i < g->e; i++)
        supp[i] = 0;

    thread_local std::unordered_map<unsigned, unsigned> neighSet (2*maxOutDeg);
    #pragma omp for
    for (unsigned i = 0; i < g->n; i++)
    {
            //hash neighbors
            for (unsigned j = g->cd[i]; j<g->cd[i]+ds[i]; j++)
            {
                if (!deleted[g->eid[j]])
                    neighSet.insert(std::make_pair(g->adj[j], g->eid[j]));
            }

            //join with neighbors of neighbors
            for (unsigned j = g->cd[i]; j<g->cd[i]+ds[i]; j++)
            {
                unsigned neigh = g->adj[j];
                unsigned e1 = g->eid[j];
                if (deleted[e1]) continue;
                for (unsigned k = g->cd[neigh]; k < g->cd[neigh]+ds[neigh]; k++)
                {
                    unsigned neighOfNeigh = g->adj[k];
                    unsigned e2 = g->eid[k];
                    if (deleted[e2]) continue;
                    std::unordered_map<unsigned, unsigned>::const_iterator got = neighSet.find(neighOfNeigh);
                    if (got == neighSet.end())
                        continue;
                    unsigned e3 = got->second;
                    __sync_fetch_and_add(&supp[e1], 1);
                    __sync_fetch_and_add(&supp[e2], 1);
                    __sync_fetch_and_add(&supp[e3], 1);
                } 
            }

            for (unsigned j = g->cd[i]; j<g->cd[i] + ds[i]; j++)
                neighSet.erase(g->adj[i]);
    }
}

bool comparePair(std::pair<unsigned, unsigned> p1, std::pair<unsigned, unsigned> p2)
{
    return p1.first > p2.first;
}


graph* extractSub(graph* dag, unsigned startV, unsigned stride, unsigned thresh)
{
    assert(stride > 0);

    bool *vExist = (bool *)malloc(dag->n*sizeof(bool));
    unsigned *ds = (unsigned *)malloc(dag->n*sizeof(unsigned));
    unsigned *dp = (unsigned *)malloc(dag->n*sizeof(unsigned));
    graph* g = (graph*) malloc(sizeof(graph));
    edge* eIdToEdge;
    unsigned maxOutDeg = 0;
    unsigned NTHRD = omp_get_num_threads();
    unsigned BS = (dag->n-1)/NTHRD + 1;

    //printf("num vertices = %d, num edges = %d\n", dag->n, dag->cd[dag->n]);

    int *supp;
    unsigned *uniqE;
    unsigned *currFrontier;
    unsigned *nxtFrontier;
    bool* inCurr;
    bool* processed;
    unsigned currFrontierSize, nxtFrontierSize;


    #pragma omp parallel num_threads(NTHRD)
    {
        unsigned tid = omp_get_thread_num();
        #pragma omp for
        for (unsigned i = 0; i < dag->n; i++)
        {
            vExist[i] = false;
            ds[i] = 0;
            dp[i] = 0;
        }

        //FINDING VERTICES IN INDUCED SUBGRAPH
        #pragma omp for
        for (unsigned i = startV; i < dag->n; i+=stride)
        {
            vExist[i] = true;
            for (unsigned j = dag->cd[i]; j < dag->cd[i+1]; j++)
            {
                vExist[dag->adj[j]] = true;
#ifdef KAR_DEBUG
                assert(dag->adj[j] < i);
#endif
            }
        }

        //#pragma omp single
        //printf("found vertices\n");

        //CREATING DEGREE ARRAYS OF VERTICES
        #pragma omp for reduction (max:maxOutDeg)
        for (unsigned i = 0; i < dag->n; i++)
        {
            if (vExist[i])
            {
                for (unsigned j = dag->cd[i]; j < dag->cd[i+1]; j++)
                {
                    unsigned neigh = dag->adj[j];
                    if (vExist[neigh])
                    {
                        ds[i]++;
                        __sync_fetch_and_add(&dp[neigh], 1);
                    }
                }
                maxOutDeg = std::max(maxOutDeg, ds[i]);
                __sync_fetch_and_add(&dp[i], ds[i]);
            }
        }


        //#pragma omp single
        //printf("computed degrees\n");

        #pragma omp barrier 
        //PREFIX SCAN START//
        #pragma omp single
        {
            g->cd = (unsigned *)malloc((dag->n+1)*sizeof(unsigned)); g->cd[0]=0;
            uniqE = (unsigned *)malloc((dag->n+1)*sizeof(unsigned)); uniqE[0]=0;
            if (dag->n < 5*NTHRD)
            { 
                serialPrefix(dp, g->cd, 1, dag->n);
                serialPrefix(ds, uniqE, 1, dag->n);
            }
        }
        
        #pragma omp barrier

        if (dag->n >= 5*NTHRD)
        {
            unsigned start = tid*BS + 1; unsigned end = std::min((tid+1)*BS, dag->n); 
            serialPrefix(dp, g->cd, start, end); serialPrefix(ds, uniqE, start, end);
            #pragma omp barrier
            #pragma omp single
            {
                thrdPrefix(g->cd, BS, NTHRD, dag->n);
                thrdPrefix(uniqE, BS, NTHRD, dag->n);
            }
            #pragma omp barrier
            applyThrdPrefix(g->cd, start, end);
            applyThrdPrefix(uniqE, start, end);
        }
        //PREFIX SCAN END//
}

    //printf("computed csr offsets. Edges = %u\n", uniqE[dag->n]);

    supp = (int  *)malloc(uniqE[dag->n]*sizeof(int));
    eIdToEdge = (edge *)malloc(uniqE[dag->n]*sizeof(edge)); 
    g->adjEid = (std::pair<unsigned, unsigned>*) malloc(g->cd[dag->n]*sizeof(std::pair<unsigned, unsigned>));
    //printf("size of edge = %d, size of adjEid = %d\n", sizeof(edge), sizeof(std::pair<unsigned, unsigned>));
    g->adj = (unsigned *)malloc(g->cd[dag->n]*sizeof(unsigned));
    g->eid = (unsigned *)malloc(g->cd[dag->n]*sizeof(unsigned));

    #pragma omp parallel num_threads(NTHRD)
    {
        unsigned tid = omp_get_thread_num();

        #pragma omp for
        for(unsigned i = 0; i < dag->n; i++)
            dp[i] = 0;

    
        //construct subgraph with fwd and bkwd edges
        #pragma omp for 
        for (unsigned i = 0; i < dag->n; i++)
        {
            if (vExist[i])
            {
                unsigned deg = 0;
                for (unsigned j = dag->cd[i]; j < dag->cd[i+1]; j++)
                {
                    unsigned neigh = dag->adj[j];
#ifdef KAR_DEBUG
                    //no duplicate edges
                    if (j > dag->cd[i]) assert (neigh != dag->adj[j-1]);
                    //no self edges
                    assert(neigh < i);
#endif

                    if (vExist[neigh])
                    {
                        unsigned edgeId = uniqE[i]+deg;

                        eIdToEdge[edgeId].s = i;
                        eIdToEdge[edgeId].t = neigh;

                        g->adjEid[g->cd[i]+deg] = std::make_pair(neigh, edgeId);

                        unsigned prev = __sync_fetch_and_add(&dp[neigh], 1);
                        g->adjEid[g->cd[neigh]+ds[neigh]+prev] = std::make_pair(i, edgeId);

#ifdef KAR_DEBUG
                        assert(edgeId < uniqE[dag->n]);
                        assert(g->cd[i]+deg < g->cd[dag->n]);
                        assert(g->cd[neigh]+ds[neigh]+prev < g->cd[dag->n]);
#endif
                        deg++;
                    }
                }
            }
        }

        //#pragma omp single
        //printf("constructed csr edge array\n");

        #pragma omp for
        for (unsigned i = 0; i < dag->n; i++)
        {
            if (vExist[i])
            {
                std::sort(g->adjEid+g->cd[i], g->adjEid+g->cd[i+1]);
                for (unsigned j = g->cd[i]; j < g->cd[i+1]; j++)
                {
                    g->adj[j] = g->adjEid[j].first;
                    g->eid[j] = g->adjEid[j].second;
                }
            }
        }

        //#pragma omp single
        //printf("sorted adjacencies\n");

    } 

    free(g->adjEid);


    g->e = g->cd[dag->n]; g->n = dag->n;
    assert(g->e/2 == uniqE[dag->n]);
    g->core = dag->core;

    #pragma omp parallel for
    for (unsigned i=0; i<g->e/2; i++)
        supp[i] = 0;

    triangleCount(g, ds, supp, maxOutDeg);

    currFrontier = (unsigned *)malloc((g->e/2)*sizeof(unsigned));
    nxtFrontier = (unsigned *)malloc((g->e/2)*sizeof(unsigned));
    inCurr = (bool *)malloc((g->e/2)*sizeof(bool));
    processed = (bool *)malloc((g->e/2)*sizeof(bool));
    currFrontierSize = 0; nxtFrontierSize = 0;
    #pragma omp parallel num_threads(NTHRD)
    {
        #pragma omp for schedule (static)
        for (unsigned i = 0; i < uniqE[g->n]; i++)
        {
            processed[i] = false;
            inCurr[i] = false;
        }

        //Remove_undesired_edges();
        trussScan(g->e/2, supp, thresh-1, currFrontier, &currFrontierSize, inCurr);
        //#pragma omp single
        //printf("edges peeled = %u\n", currFrontierSize);
        while(currFrontierSize > 0)
        {
	        PKT_processSubLevel_intersection(g, currFrontier, inCurr, currFrontierSize, supp, thresh-1, nxtFrontier, &nxtFrontierSize, processed, eIdToEdge);
            //#pragma omp single
            //printf("done level. New size = %u\n", nxtFrontierSize);
            #pragma omp for
            for (unsigned i = 0; i < nxtFrontierSize; i++)
                inCurr[nxtFrontier[i]] = true;

            #pragma omp single
            {
                unsigned* tempFrontier = currFrontier;
                currFrontier = nxtFrontier;
                nxtFrontier = tempFrontier;
                 
                currFrontierSize = nxtFrontierSize;
                nxtFrontierSize = 0;
            }

            #pragma omp barrier
            
        }
    }


    free(dp);
    free(eIdToEdge);
    free(currFrontier);
    free(nxtFrontier);
    free(inCurr);
    free(processed);
    free(uniqE);
    

    unsigned *newCd;
    unsigned *newAdj;
    unsigned *tmpAdj;
    unsigned *tmpCd;
    sharedVar = 0;
    //Reconstruct Graph
    #pragma omp parallel num_threads(NTHRD)
    {
        unsigned tid = omp_get_thread_num();
        #pragma omp for reduction (+:sharedVar)
        for (unsigned i = 0; i < dag->n; i++)
        {
            if (vExist[i])
            {
                unsigned deg = 0;
                for (unsigned j = g->cd[i]; j<g->cd[i]+ds[i]; j++)
                {
                    unsigned eid = g->eid[j];
                    if (supp[eid] >= thresh)
                        deg++;
                }
                ds[i] = deg;
                sharedVar += deg;
            }
            else
                ds[i] = 0;
        }

        //#pragma omp single
        //printf("num edges remaining = %u\n", sharedVar);

        //PREFIX SCAN
        #pragma omp single
        {
            
            newCd = (unsigned *)malloc((g->n+1)*sizeof(unsigned)); newCd[0]=0;
            if (g->n < 5*NTHRD)
                serialPrefix(ds, newCd, 1, g->n);
        }

        
        #pragma omp barrier

        if (dag->n >= 5*NTHRD)
        {
            unsigned start = tid*BS + 1; unsigned end = std::min((tid+1)*BS, g->n); 
            serialPrefix(ds, newCd, start, end);
            #pragma omp barrier
            #pragma omp single
            {
                thrdPrefix(newCd, BS, NTHRD, g->n);
            }
            #pragma omp barrier
            applyThrdPrefix(newCd, start, end);
        }
    

        #pragma omp barrier

        //PREFIX SCAN END

        #pragma omp single
        newAdj = (unsigned *)malloc(newCd[g->n]*sizeof(unsigned));
        #pragma omp barrier
        
        #pragma omp for
        for (unsigned i = 0; i < g->n; i++)
        {
            if (vExist[i])
            {
                unsigned deg = 0;
                for (unsigned j = g->cd[i]; j<g->cd[i+1]; j++)
                {
                    if (g->adj[j]>i) break;
                    unsigned eid = g->eid[j];
                    if (supp[eid] >= thresh)
                        newAdj[newCd[i] + deg++] = g->adj[j];
                }
            }
        }


    }    

    g->e = newCd[g->n];


    tmpAdj = g->adj;
    g->adj = newAdj;
    newAdj = tmpAdj;

    tmpCd = g->cd;
    g->cd = newCd;
    newCd = tmpCd;

    free(newAdj);
    free(newCd);
    free(vExist);
    free(supp);
    free(ds);
    free(g->eid);

    //printf("computed filtered graph\n");

    return g;
}



subgraph* allocsub(graph *g, unsigned char k) {
	unsigned i;
	subgraph* sg = (subgraph*) malloc(sizeof(subgraph));
	sg->n = (unsigned*) calloc(k, sizeof(unsigned));
	sg->d = (unsigned**) malloc(k * sizeof(unsigned*));
	sg->nodes = (unsigned**) malloc(k * sizeof(unsigned*));
	sg->adj = (unsigned*) malloc(g->core*g->core * sizeof(unsigned));
	for (i = 2; i < k; i++) {
		sg->d[i] = (unsigned*) malloc(g->core * sizeof(unsigned));
		sg->nodes[i] = (unsigned*) malloc(g->core * sizeof(unsigned));
	}
	sg->lab = (unsigned char*) calloc(g->core, sizeof(unsigned char));

	sg->core = g->core;
	return sg;
}

void mksub(graph* g, unsigned u, subgraph* sg, unsigned char k) {
	unsigned i, j, l, v, w;

	static unsigned *old = NULL, *mynew = NULL;//to improve
#pragma omp threadprivate(mynew,old)

	if (old == NULL) {
		mynew = (unsigned*) malloc(g->n * sizeof(unsigned));
		old = (unsigned*) malloc(g->core * sizeof(unsigned));
		for (i = 0; i < g->n; i++) {
			mynew[i] = -1;
		}
	}

	for (i = 0; i < sg->n[k - 1]; i++) {
		sg->lab[i] = 0;
	}

	j = 0;
	for (i = g->cd[u]; i < g->cd[u + 1]; i++) {
		v = g->adj[i];
		mynew[v] = j;
		old[j] = v;
		sg->lab[j] = k - 1;
		sg->nodes[k - 1][j] = j;
		sg->d[k - 1][j] = 0;//new degrees
		j++;
	}

	sg->n[k - 1] = j;

	unsigned *d0 = (unsigned*) calloc(j, sizeof(unsigned));
	for (i = 0; i < sg->n[k - 1]; i++) {//reodering adjacency list and computing new degrees
		v = old[i];
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			w = g->adj[l];
			j = mynew[w];
			if (j != -1) {
				sg->adj[sg->core*i + sg->d[k - 1][i]++] = j;
				sg->adj[sg->core*j + sg->d[k - 1][j]++] = i;
				d0[i]++;
				d0[j]++;

			}
		}
	}
	unsigned *C = (unsigned*) calloc(sg->n[k - 1], sizeof(unsigned));
	int *color = (int*) malloc(sg->n[k - 1] * sizeof(int));
	unsigned *Index = (unsigned*) malloc(sg->n[k - 1] * sizeof(unsigned));
	iddegree *ig;
	ig = (iddegree*) malloc(sg->n[k - 1] * sizeof(iddegree));
	for (i = 0; i < sg->n[k - 1]; i++)
	{
		color[i] = -1;
		ig[i].id = i;
		ig[i].degree = d0[i];
	}
	qsort(ig, sg->n[k - 1], sizeof(ig[0]), cmp);

	for (i = 0; i < sg->n[k - 1]; i++)
		Index[ig[i].id] = i;


	//color ordering
	color[0] = 0;
	int colorNum = 0;

	for (int i = 1; i < sg->n[k - 1]; i++)
	{
		int tmpdegree = ig[i].degree, tmpid = ig[i].id;

		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[sg->adj[sg->core*tmpid + j]];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < ig[0].degree + 1; j++)
			if (C[j] == 0)
			{
				color[i] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}

		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[sg->adj[sg->core*tmpid + j]];
			if (color[now] != -1)
				C[color[now]] = 0;
		}

	}

	sg->color = (unsigned*) malloc(sg->n[k - 1] * sizeof(unsigned));

	for (int i = 0; i < sg->n[k - 1]; i++)
	{
		sg->d[k - 1][i] = 0;
		sg->color[i] = color[Index[i]];
	}

	//relabel
	for (i = 0; i < sg->n[k - 1]; i++) {
		v = old[i];
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			w = g->adj[l];
			j = mynew[w];
			if (j != -1) {

				if (color[Index[i]] > color[Index[j]])
				{
					sg->adj[sg->core*i + sg->d[k - 1][i]++] = j;
				}
				else
				{
					sg->adj[sg->core*j + sg->d[k - 1][j]++] = i;
				}

			}
		}
	}


	for (i = g->cd[u]; i < g->cd[u + 1]; i++) {
		mynew[g->adj[i]] = -1;
	}
}

void kclique_thread(unsigned char l, subgraph *sg, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (i = 0; i < sg->n[2]; i++) {//list all edges
			u = sg->nodes[2][i];
			(*n) += sg->d[2][u];
			/*
			end = u*sg->core + sg->d[2][u];
			for (j = u*sg->core; j < end; j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
			*/
		}
		return;
	}

	if (l > sg->n[l])
		return;
	for (i = 0; i < sg->n[l]; i++) {
		u = sg->nodes[l][i];
        /***********************
        MY CODE HERE
        Degree based filtering
        ***********************/
		if ((sg->color[u] < l - 1) || (sg->d[l][u] < l - 1))
			continue;

		sg->n[l - 1] = 0;
		end = u*sg->core + sg->d[l][u];
		for (j = u*sg->core; j < end; j++) {//relabeling nodes and forming U'.
			v = sg->adj[j];
			if (sg->lab[v] == l) {
				sg->lab[v] = l - 1;
				sg->nodes[l - 1][sg->n[l - 1]++] = v;
				sg->d[l - 1][v] = 0;//new degrees
			}
		}
		for (j = 0; j < sg->n[l - 1]; j++) {
			v = sg->nodes[l - 1][j];
			end = sg->core*v + sg->d[l][v];
			int Index = sg->core*v;
			for (k = sg->core*v; k < end; k++) {
				w = sg->adj[k];
				if (sg->lab[w] == l - 1) {
					sg->d[l - 1][v]++;
				}
				else {
					sg->adj[k--] = sg->adj[--end];
					sg->adj[end] = w;
				}
			}
		}

		kclique_thread(l - 1, sg, n);

		for (j = 0; j < sg->n[l - 1]; j++) {//restoring labels
			v = sg->nodes[l - 1][j];
			sg->lab[v] = l;
		}

	}
}

unsigned long long kclique_main(unsigned char k, graph *g) {
	int u;
	unsigned long long n = 0;
	subgraph *sg;
#pragma omp parallel private(sg,u) reduction(+:n)
	{
		sg = allocsub(g, k);
#pragma omp for schedule(dynamic, 1) nowait
		for (u = 0; u < g->n; u++) {
			mksub(g, u, sg, k);
			kclique_thread(k - 1, sg, &n);
		}

	}
	return n;
}

unsigned long long kclique_main(unsigned char k, unsigned startV, unsigned stride, graph *g) {
	int u;
	unsigned long long n = 0;
	subgraph *sg;
#pragma omp parallel private(sg,u) reduction(+:n)
	{
		sg = allocsub(g, k);
#pragma omp for schedule(dynamic, 1) nowait
		for (u = startV; u < g->n; u+=stride) {
			mksub(g, u, sg, k);
			kclique_thread(k - 1, sg, &n);
		}

	}
	return n;
}

int main(int argc, char** argv) {
	edgelist* el;
	graph* g;

    if (argc < 4)
    {
        printf("Usage: ./DDegColNodeParallel <num_threads> <k> <graph_file>\n");
        exit(1);
    }

	unsigned char k = atoi(argv[2]);
	unsigned long long n;

	omp_set_num_threads(atoi(argv[1]));

	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	//printf("Reading edgelist from file %s\n", argv[3]);

	el = readedgelist(argv[3]);
	//printf("Number of nodes = %u\n", el->n);
	//printf("Number of edges = %u\n", el->e);

	t2 = time(NULL);
	//printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	//printf("Building the graph structure\n");
	ord_core(el);
	relabel(el);
	g = mkgraph(el);

	//printf("Number of nodes (degree > 0) = %u\n", g->n);

	free_edgelist(el);

	t2 = time(NULL);
	//printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	//printf("Iterate over all cliques\n");


    unsigned stride = 3;
    graph* gFilt;
    n = 0;
    for (unsigned i = 0; i < stride; i++)
    {
        gFilt = extractSub(g, i, stride, k-2); 
        unsigned long long locCount = kclique_main(k, i, stride, gFilt);
        n += locCount;
        
        std::cout << "n count: " << locCount << "val: " << i << std::endl;

        free_graph(gFilt);
    }

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	//printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	free_graph(g);

	//printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}
