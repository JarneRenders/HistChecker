/**
 * histChecker.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE "Usage: ./histChecker [-1|-f] [-c] [-hv]"
#define HELPTEXT "Helptext: Program for checking HISTs in a graph.\n\
\n\
Graph are read from stdin in graph6 format. Graphs are sent to stdout in graph6\n\
format. For more information on the format, see \n\
http://users.cecs.anu.edu.au/~bdm/data/formats.txt.\n\
\n\
Without any parameters the program outputs graphs which are HIST-free to stdout.\n\
\n\
  -1  : outputs graphs which are HIST-critical; cannot be used with -f\n\
  -f  : count the number of HISTs in each graph and give information on the\n\
        smallest counts encountered; combine with -o# to output graphs with a\n\
        specific number of HISTs\n\
  -o# : use only with -f; outputs the graphs with # HISTs\n\
  -c  : outputs the complement of graphs which would normally be output\n\
  -h  : print this helptext\n\
  -v  : output extra information to stderr\n"

#define HISTCOUNTSTOKEEP 3

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include "readGraph/readGraph6.h"
#include "bitset.h"

//******************************************************************************
//
//                          Structs
//
//******************************************************************************

struct graph {
    int numberOfVertices;
    int numberOfEdges;
    bitset *adjacencyList;
    int *edgeIndices; // Mapping of edge to index
    int *indexToEdge; // Mapping of index to edge
};

struct tree {
    int numberOfVertices;
    int numberOfEdges;
    bitset usedVertices;
    bitset availableVertices;
    bitset *adjacencyList;
};

struct options {
    bool complementFlag;
    bool countFlag;
    bool hasModResPair;
    bool histCriticalFlag;
    bool outputFlag; // Only if countFlag
    bool verboseFlag;
};

struct counters {
    long long unsigned int skippedGraphs;
    long long unsigned int fewestNumberOfHISTs;
};

//******************************************************************************
//
//                  Macros for dealing with graphs
//
//******************************************************************************

//  Add one edge. 
#define addEdge(g,i,j) {\
 add((g)->adjacencyList[i], j); add((g)->adjacencyList[j],i);\
 (g)->numberOfEdges++;\
}

//  Remove one edge.
#define removeEdge(g,i,j) {\
 removeElement((g)->adjacencyList[i], j);\
 removeElement((g)->adjacencyList[j],i);\
 (g)->numberOfEdges--;\
}

//******************************************************************************
//
//                      Writing graphs
//
//******************************************************************************

// Readable format to stderr
void printGraph(bitset adjacencyList[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

// Readable format to stderr
void printTree(struct tree *T) {
    printGraph(T->adjacencyList, T->numberOfVertices);
}

// Graph6 format to stdout
void writeToG6(bitset adjacencyList[], int numberOfVertices) {
    char graphString[8 + numberOfVertices*(numberOfVertices - 1)/2];
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(numberOfVertices <= 62) {
        graphString[pointer++] = (char) numberOfVertices + 63;
    }
    else if(numberOfVertices <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else if(numberOfVertices <= 68719476735) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's 
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < numberOfVertices; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(contains(adjacencyList[i], j)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline and end of string character.
    graphString[pointer++] = '\n';
    graphString[pointer++] = '\0';
    printf("%s", graphString);
}

//******************************************************************************
//
//               Methods for creating and destroying graphs 
//
//******************************************************************************

// Compute the mapping of edges to indices. Every edge of the graph G gets a
// label from 0 to |E(G)| - 1.
void initEdgeIndices(struct graph *g) {
    g->edgeIndices = 
     malloc(g->numberOfVertices * g->numberOfVertices * sizeof(int));
    g->indexToEdge = 
     malloc(g->numberOfVertices * (g->numberOfVertices - 1) * sizeof(int));
    int idx = 0;
    for(int i = 0; i < g->numberOfVertices; i++) {
        forEachAfterIndex(j, g->adjacencyList[i], i) {
            g->edgeIndices[g->numberOfVertices * i + j] = idx;
            g->edgeIndices[g->numberOfVertices * j + i] = idx;
            g->indexToEdge[2*idx] = i;
            g->indexToEdge[2*idx + 1] = j;
            idx++;
        }
    }
    g->numberOfEdges = idx;
}

// Load the graph struct g with the correct data obtained from the graphString
// in graph6 format.
int readGraph(const char *graphString, struct graph *g, struct options *options,
 struct counters *counters) {

    g->numberOfVertices = getNumberOfVertices(graphString);
    if(g->numberOfVertices == -1 || g->numberOfVertices > MAXBITSETSIZE) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph! Too many vertices.\n");
        }
        counters->skippedGraphs++;
        return 1;
    }

    g->adjacencyList = malloc(sizeof(bitset)*g->numberOfVertices);
    if(loadGraph(graphString, g->numberOfVertices, g->adjacencyList) == -1) {
        if(options->verboseFlag){
            fprintf(stderr,
             "Skipping invalid graph! Some error occurred while reading.\n");
        }
        counters->skippedGraphs++;
        return 1;
    }

    initEdgeIndices(g);

    return 0;
}

void createEmptyTree(struct graph *g, struct tree *T) {
    T->numberOfVertices = g->numberOfVertices; 
    T->adjacencyList = calloc(g->numberOfVertices, sizeof(bitset));
    T->usedVertices = complement(EMPTY, g->numberOfVertices);
    T->numberOfEdges = 0;
    T->availableVertices = EMPTY;
}

void freeGraph(struct graph *g) {
    free(g->adjacencyList);
    free(g->edgeIndices);
    free(g->indexToEdge);
}

void freeTree(struct tree *T) {
    free(T->adjacencyList);
}

//******************************************************************************
//
//                      Methods for choosing next edge
//
//******************************************************************************

//  Find the vertex v in set for which d_g(v) - d_T(v) is the smallest.
int getVertexWithSmallestDegree(struct graph *g, struct tree *T, bitset set) {

    if(isEmpty(set)) {
        return -1;
    }

    int smallest = next(set, -1);
    int smallestDifference = 
     size(intersection(g->adjacencyList[smallest], T->usedVertices)) - 
     size(T->adjacencyList[smallest]);

    forEachAfterIndex(el, set, smallest) { 

        int difference = 
         size(intersection(g->adjacencyList[el], T->usedVertices)) -
         size(T->adjacencyList[el]);
        if(smallestDifference > difference) {
            smallest = el;
            smallestDifference = difference;
        }
    }

    return smallest;
}

//  The next edge e we consider will be incident to a vertex v in which
//  d_g(v) - d_T(v) is smallest and incident to the neighbour w of v for which
//  d_g(w) - d_T(w) is smallest such that e is not yet in T. 
bool getNextEdge(struct graph *g, struct tree *T, int edge[],
 bool *bothVerticesInTree) {

    edge[0] = getVertexWithSmallestDegree(g, T, T->availableVertices);
    if(edge[0] == -1) {
        return false;
    }

    bitset neighbours = 
     difference(intersection(g->adjacencyList[edge[0]], T->usedVertices),
      T->adjacencyList[edge[0]]);

    edge[1] = getVertexWithSmallestDegree(g, T, neighbours);
    if(edge[1] == -1) {
        return false;
    }

    // We will only remove e if the following is true.
    *bothVerticesInTree = size(T->adjacencyList[edge[1]]);

    return true;
}

//******************************************************************************
//
//                Methods for the recursive algorithm
//
//******************************************************************************

void recursion(struct graph *g, struct tree *T, struct options *options,
 struct counters *counters, long long unsigned int *counter, 
 long long unsigned int *intermediateGraphs, bitset *verticesDeletedHISTNotFound);

// Every time we add an edge to T or remove one from g, we need to update which
// vertices still have incident edges in g not yet in T. These are called the
// availableVertices.
void updateAvailableVertices(struct graph *g, struct tree *T, int edge[]) {

    for(int i = 0; i < 2; i++) {

        int degreeInT = size(T->adjacencyList[edge[i]]);
        int degreeInG = 
         size(intersection(g->adjacencyList[edge[i]], T->usedVertices));

        if(degreeInT > 0 && (degreeInG > degreeInT)) {
            add(T->availableVertices, edge[i]);
            continue;
        }
        removeElement(T->availableVertices, edge[i]);
    }
}

// Check if T is homeomorphically irreducible.
bool isHIT(struct tree *T) {
    forEach(i, T->usedVertices) {
        if(size(T->adjacencyList[i]) == 2) {
            return false;
        }
    }
    return true;
}

//  Add edge to tree, check if spanning, if not, either prune or continue with
//  recursive algorithm.
void addEdgeTree(struct graph *g, struct tree *T, struct options *options,
 struct counters *counters, int edge[], long long unsigned int *counter, 
 long long unsigned int *intermediateGraphs, 
 bitset *verticesDeletedHISTNotFound) {

    addEdge(T, edge[0], edge[1]);
    bitset origAvailable = T->availableVertices;
    updateAvailableVertices(g, T, edge);

    //  Check if spanning and HIST
    if(T->numberOfEdges == size(T->usedVertices) - 1) {
        if(isHIT(T)) {
            (*counter)++;
            // writeToG6(T->adjacencyList, T->numberOfVertices);
            if(options->verboseFlag) {
                printTree(T);
            }
        }
        removeEdge(T, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    //  If checking for histCriticals, we will first check if G is HIST-free. If we
    //  encounter any HIST of G-v for some v, store that we found it.
    if(options->histCriticalFlag && T->numberOfEdges == g->numberOfVertices - 2) {
        if(isHIT(T)) {

            //  Check for all vertices which has no neighbours in the tree.
            //  There must be at least one since T has n-2 edges. Might be
            //  faster if we dynamically save the not yet used vertices, but it
            //  is not really a bottleneck.
            int unusedVertex = -1;
            for(int i = 0; unusedVertex == -1; i++) {
                if(size(T->adjacencyList[i]) == 0) {
                    unusedVertex = i;
                }
            }
            if(options->verboseFlag &&
             contains(*verticesDeletedHISTNotFound, unusedVertex)) {
                fprintf(stderr, "For G - %d:\n", unusedVertex);
                printTree(T);
            }
            removeElement(*verticesDeletedHISTNotFound, unusedVertex);
        }
    }

    // Check if can prune: if vertex has degree 2 in graph and in tree it will
    // always have degree 2 in tree, so we can prune here. Only need to check
    // it for endpoints of last added edge. 
    if((size(intersection(g->adjacencyList[edge[0]], T->usedVertices)) == 2 && 
        size(T->adjacencyList[edge[0]]) == 2) ||
     (size(intersection(g->adjacencyList[edge[1]], T->usedVertices)) == 2 &&
      size(T->adjacencyList[edge[1]]) == 2)) {
        removeEdge(T, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    recursion(g, T, options, counters, counter, intermediateGraphs,
     verticesDeletedHISTNotFound);

    //  Reset
    removeEdge(T, edge[0], edge[1]);
    T->availableVertices = origAvailable;
}

//  If a g has a degree 2 vertex v (e.g. after removal of some edge), and v has
//  degree 1 in T, then we can remove other edge from graph. ONLY PERFORM IF
//  deg_T(v) != 2
void handleDegree2Vertex(struct graph *g, struct tree *T,
 struct options *options, int v, struct counters *counters, 
 long long unsigned int *counter, long long unsigned int *intermediateGraphs, 
 bitset *verticesDeletedHISTNotFound) {

    int nbr1 = next(intersection(g->adjacencyList[v], T->usedVertices), -1);
    int nbr2 = next(intersection(g->adjacencyList[v], T->usedVertices), nbr1);
    int edge1[2] = {v, nbr1};
    int edge2[2] = {v, nbr2};
    bitset origAvailable = T->availableVertices;

    //  We assume that at most one edge incident to v belongs to tree. We are
    //  sure of this since we only call this function after doing the pruning
    //  step, which will prune if v has deg 2 in both g and T.
    if(contains(T->adjacencyList[v], nbr1)) {
        removeEdge(g, edge2[0], edge2[1]);
        updateAvailableVertices(g, T, edge2);

        recursion(g, T, options, counters, counter, intermediateGraphs,
         verticesDeletedHISTNotFound);

        //  Reset
        addEdge(g, edge2[0], edge2[1]);
        T->availableVertices = origAvailable;
        return;
    }

    if(contains(T->adjacencyList[v], nbr2)) {
        removeEdge(g, edge1[0], edge1[1]);
        updateAvailableVertices(g, T, edge1);

        recursion(g, T, options, counters, counter, intermediateGraphs,
         verticesDeletedHISTNotFound);

        //  Reset
        addEdge(g, edge1[0], edge1[1]);
        T->availableVertices = origAvailable;
        return;
    }

}

//  Remove edge from graph. No changes are made to T so do not need to check if
//  it is spanning. Check if we can prune or if we get a degree 2 vertex, if
//  not, continue with recursive algorithm.
void removeEdgeTree(struct graph *g, struct tree *T, struct options *options,
 struct counters *counters, int edge[], long long unsigned int *counter, 
 long long unsigned int *intermediateGraphs, 
 bitset *verticesDeletedHISTNotFound) {

    removeEdge(g, edge[0], edge[1]);
    bitset origAvailable = T->availableVertices;
    updateAvailableVertices(g, T, edge);

    //  Check if can prune: we can if removing edge created isolated vertex. Of
    //  course, we only need to check the endpoints of the removed edge.
    int degreeVertex1 = 
     size(intersection(g->adjacencyList[edge[0]], T->usedVertices));
    int degreeVertex2 = 
     size(intersection(g->adjacencyList[edge[1]], T->usedVertices));

    if(degreeVertex1 == 0 || degreeVertex2 == 0) {
        addEdge(g, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    //  Check if can prune: we can if there is a vertex v which has degree 2 in
    //  g and in T. Only need to check endpoints of removed edge.
    if((degreeVertex1 == 2 && size(T->adjacencyList[edge[0]]) == 2) ||
     (degreeVertex2 == 2 && size(T->adjacencyList[edge[1]]) == 2)) {
        addEdge(g, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    //  If one of the endpoints of the removed edge became a degree 2 vertex,
    //  and it has degree 1 in tree, we can remove the other edge incident to
    //  the degree 2 vertex. It is important we do this AFTER the pruning step
    //  above. 
    if(degreeVertex1 == 2 && contains(T->availableVertices, edge[0])) {
        handleDegree2Vertex(g, T, options, edge[0], counters, counter,
         intermediateGraphs, verticesDeletedHISTNotFound);
        addEdge(g, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    if(degreeVertex2 == 2 && contains(T->availableVertices, edge[1])) {
        handleDegree2Vertex(g, T, options, edge[1], counters, counter,
         intermediateGraphs, verticesDeletedHISTNotFound);
        addEdge(g, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    recursion(g, T, options, counters, counter, intermediateGraphs,
     verticesDeletedHISTNotFound);

    //  Reset
    addEdge(g, edge[0], edge[1]);
    T->availableVertices = origAvailable;
}

void recursion(struct graph *g, struct tree *T, struct options *options,
 struct counters *counters, long long unsigned int *counter, 
 long long unsigned int *intermediateGraphs, 
 bitset *verticesDeletedHISTNotFound) {

    // If countFlag is false, we are checking for HIST-freeness, so return when
    // we have found 1 HIST.
    (*intermediateGraphs)++; // Counts number of recursive calls.
    if(!options->countFlag && (*counter)) { 
        return;
    }
    // When trying the find graph with the fewest number of HISTs in a list, we
    // can stop once the number of HISTs of the current graph is already more
    // than the number of HISTs of the graphs we want to store. E.g. if we want
    // to store the 3 graphs with fewest HISTs from a list and already found
    // ones with 100,120,140 HISTs, then we can stop once we found the 140th
    // HIST of the current graph.  
    // However, if outputting graphs with specific number of HISTs, we cannot
    // stop yet. Could be made more efficient if we stop at the maximum of
    // fewestNumberOfHISTs and the number of HISTs to output.
    else { 
        if((*counter) > counters->fewestNumberOfHISTs && !options->outputFlag) {
            return;
        }
    }

    int edge[2];
    bool bothVerticesInTree = false;

    //  Check if there is a next edge and store in edge[]
    if(getNextEdge(g, T, edge, &bothVerticesInTree)) {

        //  If edge makes cycle in T do not add it to T.
        if(!bothVerticesInTree) {
            addEdgeTree(g, T, options, counters, edge, counter,
             intermediateGraphs, verticesDeletedHISTNotFound);
        }

        removeEdgeTree(g, T, options, counters, edge, counter,
         intermediateGraphs, verticesDeletedHISTNotFound);
    }
}

long long unsigned int initGeneration(struct graph *g, struct options *options,
 struct counters *counters, bitset forbiddenVertices, 
 bitset *verticesDeletedHISTNotFound) {

    struct tree T;
    createEmptyTree(g, &T); 
    long long unsigned int counter = 0;
    long long unsigned int intermediateGraphs = 0;
    T.usedVertices = difference(T.usedVertices, forbiddenVertices);
    if(size(T.usedVertices) == 0) {
        freeTree(&T);
        return 0;
    }

    //  Get the first allowed vertex v and denote its neighbours by u1,...uk.
    //  First we start the recursion for all trees with vu1, then for all trees
    //  without vu1 with vu2, then for all trees without vu1,vu2 and with vu3
    //  etc. We will obtain all HISTs in this way.
    int vertex1 = next(T.usedVertices, -1);
    bitset nbrs = g->adjacencyList[vertex1]; // Make copy of original nbrs.

    // Check if starting vertex is not isolated.
    if(isEmpty(intersection(nbrs, T.usedVertices))) {
        return 0;
    }

    forEach(vertex2, intersection(nbrs, T.usedVertices)) {

        int edge[2] = {vertex1, vertex2};
        addEdge(&T, vertex1, vertex2);
        updateAvailableVertices(g,&T,edge);

        //  If there are 2 usedVertices, then we already have a spanning tree.
        if(T.numberOfEdges == size(T.usedVertices) - 1) {
            freeTree(&T);
            return 1;
        }

        recursion(g, &T, options, counters, &counter, &intermediateGraphs,
         verticesDeletedHISTNotFound);

        removeEdge(&T, vertex1, vertex2);
        updateAvailableVertices(g,&T,edge);
        
        removeEdge(g, vertex1, vertex2);
        updateAvailableVertices(g,&T,edge);

    } 

    //  Restore our original graph
    forEach(nbr, nbrs) {
        addEdge(g, vertex1, nbr);
    }

    if(options->verboseFlag && options->countFlag) {
        fprintf(stderr, "Number of HISTs found: %llu\n", counter);
        fprintf(stderr, "Times recursion step was performed: %llu\n",
         intermediateGraphs);
    }

    freeTree(&T);
    return counter;
}

//******************************************************************************
//
//            Methods for keeping track of graphs with few HISTs
//
//******************************************************************************

//  When we are done counting the HISTs of a graph, check whether or not it
//  should be in the list of graphs with fewest HISTs.
void updateHISTFrequencies(struct counters *counters, 
 long long unsigned int fewestHISTFrequencies[], 
 long long unsigned int nHISTS) {

    //  Update frequency if number of HISTs is in array or replace largest
    //  number of HISTs in array by nHISTS if nHISTS is smaller.
    int idxWithMostHISTS = 0;
    for(int i = 0; i < HISTCOUNTSTOKEEP; i++) {
        if(fewestHISTFrequencies[2*i] == nHISTS) {
            fewestHISTFrequencies[2*i+1]++;
            return;
        }
        if(fewestHISTFrequencies[2*i] > fewestHISTFrequencies[idxWithMostHISTS]) {
            idxWithMostHISTS = 2*i;
        }
    }
    if(fewestHISTFrequencies[idxWithMostHISTS] > nHISTS) {
        fewestHISTFrequencies[idxWithMostHISTS] = nHISTS;
        fewestHISTFrequencies[idxWithMostHISTS+1] = 1;
    }

    //  Find most number of HISTs in our array
    //  Not efficient, but not a bottleneck
    idxWithMostHISTS = 0;
    for(int i = 0; i < HISTCOUNTSTOKEEP; i++) {
        if(fewestHISTFrequencies[2*i] > fewestHISTFrequencies[idxWithMostHISTS]) {
            idxWithMostHISTS = 2*i;
        }
    }
    counters->fewestNumberOfHISTs = fewestHISTFrequencies[idxWithMostHISTS];
}

//******************************************************************************
//
//               Methods for starting the various algorithms
//
//******************************************************************************

//  Count number of HISTs and keep track of the HISTCOUNTSTOKEEP graphs with
//  fewest number of HISTs.
long long unsigned int countHISTs(struct graph *g, struct options *options,
 struct counters *counters, long long unsigned int fewestHISTFrequencies[]) {

    bitset empty = EMPTY;
    long long unsigned int nHISTS = 
     initGeneration(g, options, counters, EMPTY, &empty);
    updateHISTFrequencies(counters, fewestHISTFrequencies, nHISTS);

    return nHISTS;
}

//  initGeneration will return 1 if a HIST has been found and 0 otherwise in
//  this case.
bool isHISTFree(struct graph *g, struct options *options,
 struct counters *counters) {
    bitset empty = EMPTY;
    return !initGeneration(g, options, counters, EMPTY, &empty);
}

//  HIST-free but vertex-deleted subgraphs do contain a HIST.
bool isHISTCritical(struct graph *g, struct options *options,
 struct counters *counters) {

    //  Do HIST-freeness check. Remove v from leftToCheckForHIST for any HIST of
    //  G-v we accidentally encounter already.
    bitset leftToCheckForHIST = complement(EMPTY, g->numberOfVertices);
    if(initGeneration(g, options, counters, EMPTY, &leftToCheckForHIST)) {
        if(options->verboseFlag) {
            fprintf(stderr, "%s\n", "Not HIST-free.");
        }
        return false;
    }    

    //  Check whether G-v contains a HIST for those that were not yet found.
    bitset empty = EMPTY;
    if(options->verboseFlag) {
        fprintf(stderr, "Still need to check G - v for v in: ");
        forEach(i, leftToCheckForHIST) {
            fprintf(stderr, "%d ", i);
        }
        fprintf(stderr, "\n\n");
    }
    forEach(i, leftToCheckForHIST) {
        if(options->verboseFlag) {
            fprintf(stderr, "For G - %d:\n", i);
        }
        if(!initGeneration(g, options, counters, singleton(i), &empty)) {
            if(options->verboseFlag) {
                fprintf(stderr, "G - %d does not have a HIST\n",i);
            }
            return false;
        }
    }
    return true;
}

int main(int argc, char ** argv) {
    struct counters counters = {0};
    struct options options = {0};
    int opt;
    int mod = 1;
    int res = 0;
    long long unsigned int output = 0;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"HIST-critical", no_argument, NULL, '1'},
            {"complement", no_argument, NULL, 'c'},
            {"fewest-HISTs", no_argument, NULL, 'f'},
            {"help", no_argument, NULL, 'h'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "1cfho:v", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case '1':
                options.histCriticalFlag = true;
                break;
            case 'c':
                options.complementFlag = true;
                break;
            case 'f':
                options.countFlag = true;
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'o':
                options.outputFlag = true;
                output = strtoull(optarg, (char **)NULL, 10);
                if(output < 0) {
                    fprintf(stderr, "Error: output value should be >= 0.\n");
                    return 1;
                }
                break;
            case 'v':
                options.verboseFlag = true;
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "\n%s\n\n", USAGE);
                fprintf(stderr,
                 "Use ./histChecker -h for more detailed instructions.\n");
                return 1;
        }
    }

    //  Loop over non-option arguments.
    while (optind < argc) {
        if(sscanf(argv[optind], "%d/%d", &res, &mod) == 2) {
            if(options.hasModResPair) {
                fprintf(stderr,
                 "Error: You can only add one mod/res pair as an argument.\n");
                fprintf(stderr, "\n%s\n\n", USAGE);
                fprintf(stderr,
                 "Use ./histChecker -h for more detailed instructions.\n");
                return 1;
            }
            options.hasModResPair = true;
        }
        else {
            fprintf(stderr,"Error: Unknown argument: %s\n", argv[optind]);
            fprintf(stderr, "\n%s\n\n", USAGE);
            fprintf(stderr,
             "Use ./histChecker -h for more detailed instructions.\n");
            return 1;
        }
        optind++;
    }

    unsigned long long int total = 0;
    unsigned long long int counter = 0;
    unsigned long long int passedGraphs = 0;
    unsigned long long int fewestHISTFrequencies[2*HISTCOUNTSTOKEEP];
    char *graphStringsWithFewestHISTs[HISTCOUNTSTOKEEP] = {NULL};
    clock_t start = clock();

    if(options.countFlag) {
        for(int i = 0; i < HISTCOUNTSTOKEEP; i++) {
            fewestHISTFrequencies[2*i] = ULLONG_MAX;
        }
        counters.fewestNumberOfHISTs = ULLONG_MAX;
    } 

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {

        if(total++ % mod != res) {
            continue;
        }

        struct graph g;
        int errCode = readGraph(graphString, &g, &options, &counters);
        if(errCode != 0) {
            continue;
        }
        counter++;

        if(options.verboseFlag) {
            fprintf(stderr, "Looking at:\n");
            printGraph(g.adjacencyList, g.numberOfVertices);
        }

        if(options.histCriticalFlag) {
            if(isHISTCritical(&g, &options, &counters)) {
                if(!options.complementFlag) {
                    passedGraphs++;
                    printf("%s", graphString);
                }
            }
            else if(options.complementFlag) {
                passedGraphs++;
                printf("%s", graphString);
            }
            freeGraph(&g);
            continue;
        }

        if(options.countFlag) {
            long long unsigned int nHISTS;
            nHISTS = countHISTs(&g, &options, &counters, fewestHISTFrequencies);
            for(int i = 0; i < HISTCOUNTSTOKEEP; i++) {
                if(fewestHISTFrequencies[2*i] == nHISTS) {

                    // Free previously saved strings
                    if(graphStringsWithFewestHISTs[i]) {
                        free(graphStringsWithFewestHISTs[i]);
                    }

                    // Malloc copy of graphString
                    graphStringsWithFewestHISTs[i] = strdup(graphString); 
                }
            }
            if(options.outputFlag && output == nHISTS) {
                if(!options.complementFlag) {
                    passedGraphs++;
                    printf("%s", graphString);
                }
            }
            else if(options.complementFlag) {
                passedGraphs++;
                printf("%s", graphString);

            }
            freeGraph(&g);
            continue;
        }

        if(isHISTFree(&g, &options, &counters)) {
            if(!options.complementFlag) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }
        else if(options.complementFlag) {
            passedGraphs++;
            printf("%s", graphString);
        }

        freeGraph(&g);
    }
    free(graphString);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    if(options.countFlag) {
        fprintf(stderr, "#Graphs with fewest HISTs:\n");
        for(int i = 0; i < HISTCOUNTSTOKEEP; i++) {
            if(fewestHISTFrequencies[2*i] != ULLONG_MAX) {
                fprintf(stderr, "\tHaving %llu HISTS: %llu graphs.",
                 fewestHISTFrequencies[2*i], fewestHISTFrequencies[2*i+1]);
                fprintf(stderr, " E.g.: %s", graphStringsWithFewestHISTs[i]);
                free(graphStringsWithFewestHISTs[i]);
            }
        }
        fprintf(stderr, "\n");
    }

    if(options.histCriticalFlag) {
        fprintf(stderr,
         "Checked %lld graphs in %f seconds: %llu are %sHIST-critical.\n",
         counter, time_spent, passedGraphs, 
         options.complementFlag ? "not " : "");
    }
    else if(options.countFlag) {
        if(options.outputFlag) {
            fprintf(stderr,
             "Checked %lld graphs in %f seconds: %llu %s %llu HISTs.\n",
             counter, time_spent, passedGraphs,
             options.complementFlag ? "do not have" : "have",
             output);
        }            
        else {
            fprintf(stderr,
             "Checked %lld graphs in %f seconds.\n",
             counter, time_spent);
        }
    }
    else {
        fprintf(stderr,
         "Checked %lld graphs in %f seconds: %llu are %sHIST-free.\n",
         counter, time_spent, passedGraphs, 
         options.complementFlag ? "not " : "");
    }

    if(counters.skippedGraphs > 0) {
        fprintf(stderr,
         "Warning: %lld graphs were skipped. Use -v for more information.\n",
         counters.skippedGraphs);
    }

    return 0;
}