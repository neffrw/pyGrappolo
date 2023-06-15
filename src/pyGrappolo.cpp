#include <pybind11/pybind11.h>

// #include "grappolo/DefineStructure/defs.h"
#include "defs.h"
#include "input_output.h"
#include "basic_comm.h"
#include "basic_util.h"
#include "utilityClusteringFunctions.h"
#include "color_comm.h"
#include "sync_comm.h"

namespace py = pybind11;

// Create function that converts python dict to graph type
void list_to_graph_undirected(graph* G, const long NV, const long NE,
                                const py::list &pyedgeList) {
    printf("Parsing a NetworkX undirected edgelist...\n");
    printf("WARNING: Assumes that the graph is undirected -- an edge is stored ONLY ONCE!.\n");
    int nthreads = 0;
    
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    double time1, time2;

    ///////////
    time1 = omp_get_wtime();
    long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
    assert(edgeListPtr != NULL);
    edge *edgeList = (edge *) malloc(2 * NE * sizeof(edge)); //Every edge stored twice
    assert( edgeList != NULL);
    time2 = omp_get_wtime();
    printf("Time for allocating memory for storing graph = %lf\n", time2 - time1);

#pragma omp parallel for
    for (long i=0; i <= NV; i++)
        edgeListPtr[i] = 0; //For first touch purposes
    
    //Build the EdgeListPtr Array: Cumulative addition
    time1 = omp_get_wtime();
// #pragma omp parallel for
    for(long i=0; i<NE; i++) {
        __sync_fetch_and_add(&edgeListPtr[pyedgeList[i].cast<py::tuple>()[0].cast<long>()+1], 1); //Leave 0th position intact
        __sync_fetch_and_add(&edgeListPtr[pyedgeList[i].cast<py::tuple>()[1].cast<long>()+1], 1); //Leave 0th position intact
    }
    for (long i=0; i<NV; i++) {
        edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum
    }
    //The last element of Cumulative will hold the total number of characters
    time2 = omp_get_wtime();
    printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
    printf("Sanity Check: |E| = %ld, edgeListPtr[NV]= %ld\n", 2*NE, edgeListPtr[NV]);
    assert(2*NE == edgeListPtr[NV]);

    //Keep track of how many edges have been added for a vertex:
    printf("About to allocate for added vector: %ld\n", NV);
    long  *added  = (long *)  malloc(NV  * sizeof(long));
    assert( added != NULL);
#pragma omp parallel for
    for (long i = 0; i < NV; i++)
        added[i] = 0;
    
    printf("About to build edgeList...\n");
    //Build the edgeList from edgeListTmp:
// #pragma omp parallel for
    for(long i=0; i<NE; i++) {
        long head      = pyedgeList[i].cast<py::tuple>()[0].cast<long>();
        long tail      = pyedgeList[i].cast<py::tuple>()[1].cast<long>();
        const py::dict weightDict = pyedgeList[i].cast<py::tuple>()[2].cast<py::dict>();
        double weight  = weightDict["weight"].cast<double>();
        
        long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);
        edgeList[Where].head = head;
        edgeList[Where].tail = tail;
        edgeList[Where].weight = weight;
        //Add the other way:
        Where = edgeListPtr[tail] + __sync_fetch_and_add(&added[tail], 1);
        edgeList[Where].head = tail;
        edgeList[Where].tail = head;
        edgeList[Where].weight = weight;
    }
    //time2 = omp_get_wtime();
    printf("Time for building edgeList = %lf\n", time2 - time1);
    
    G->sVertices    = NV;
    G->numVertices  = NV;
    G->numEdges     = NE;
    G->edgeListPtrs = edgeListPtr;
    G->edgeList     = edgeList;
    
    //Clean up:
    free(added);
}

// Create python function to accept arguments
py::list grappolo(const py::dict &input_graph, const py::kwargs& kwargs) {

    // Parse kwargs into clustering_parameters
    clustering_parameters opts;
    opts.inFile = input_graph["inFile"].cast<std::string>().c_str();
    opts.minGraphSize = input_graph["minGraphSize"].cast<long>();
    opts.threshold = input_graph["threshold"].cast<double>();
    opts.C_thresh = input_graph["C_thresh"].cast<double>();
    opts.percentage = input_graph["percentage"].cast<int>();
    opts.numColors = input_graph["numColors"].cast<int>();
    opts.ftype = input_graph["ftype"].cast<int>();
    opts.output = input_graph["output"].cast<bool>();
    opts.strongScaling = input_graph["strongScaling"].cast<bool>();
    opts.VF = input_graph["VF"].cast<bool>();
    opts.coloring = input_graph["coloring"].cast<int>();
    opts.replaceMap = input_graph["replaceMap"].cast<bool>();
    opts.syncType = input_graph["syncType"].cast<int>();
    opts.basicOpt = input_graph["basicOpt"].cast<int>();

    int nT = 1; //Default is one thread
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
    if (nT < 1) {
        // printf("The number of threads should be greater than one.\n");
        throw std::runtime_error("The number of threads should be greater than one.");
    }
    // File Loading
    double time1, time2;
    graph* G = (graph *) malloc (sizeof(graph));
    const long NV = input_graph["numVertices"].cast<long>();
    const long NE = input_graph["numEdges"].cast<long>();
    list_to_graph_undirected(G, NV, NE, input_graph["edgeList"]);

    displayGraphCharacteristics(G);
    int threadsOpt = 0;
    if(opts.threadsOpt)
        threadsOpt =1;
    threadsOpt =1;
    
    int replaceMap = 0;
    if(opts.basicOpt == 1)
        replaceMap = 1;
    
    /* Vertex Following option */
    if( opts.VF ) {
        printf("Vertex following is enabled.\n");
        time1 = omp_get_wtime();
        long numVtxToFix = 0; //Default zero
        long *C = (long *) malloc (G->numVertices * sizeof(long)); assert(C != 0);
        numVtxToFix = vertexFollowing(G,C); //Find vertices that follow other vertices
        if( numVtxToFix > 0) {  //Need to fix things: build a new graph
            printf("Graph will be modified -- %ld vertices need to be fixed.\n", numVtxToFix);
            graph *Gnew = (graph *) malloc (sizeof(graph));
            long numClusters = renumberClustersContiguously(C, G->numVertices);
            buildNewGraphVF(G, Gnew, C, numClusters);
            //Get rid of the old graph and store the new graph
            free(G->edgeListPtrs);
            free(G->edgeList);
            free(G);
            G = Gnew;
        }
        free(C); //Free up memory
        printf("Graph after modifications:\n");
        displayGraphCharacteristics(G);
    }//End of if( VF == 1 )

    // Datastructures to store clustering information
    long *C_orig = (long *) malloc (NV * sizeof(long)); assert(C_orig != 0);
    graph* G_orig = (graph *) malloc (sizeof(graph)); //The original version of the graph
    duplicateGivenGraph(G, G_orig);
    
    //Call the clustering algorithm:
    //Call the clustering algorithm:
    if ( opts.strongScaling ) { //Strong scaling enabled
        //Retain the original copy of the graph:
        graph* G_original = (graph *) malloc (sizeof(graph)); //The original version of the graph
        time1 = omp_get_wtime();
        duplicateGivenGraph(G, G_original);
        time2 = omp_get_wtime();
        printf("Time to duplicate : %lf\n", time2-time1);
        
        //Run the algorithm in powers of two for the maximum number of threads available
        int curThread = 2; //Start with two threads
        while (curThread <= nT) {
            printf("\n\n***************************************\n");
            printf("Starting run with %d threads.\n", curThread);
            printf("***************************************\n");
            //Call the clustering algorithm:
#pragma omp parallel for
            for (long i=0; i<G->numVertices; i++) {
                C_orig[i] = -1;
            }
            if(opts.coloring != 0){
                runMultiPhaseColoring(G, C_orig, opts.coloring, opts.numColors, replaceMap, opts.minGraphSize, opts.threshold, opts.C_thresh, curThread, threadsOpt);
            }else if(opts.syncType != 0){
                runMultiPhaseSyncType(G, C_orig, opts.syncType, opts.minGraphSize, opts.threshold, opts.C_thresh, curThread, threadsOpt);
            }else{
                runMultiPhaseBasic(G, C_orig, opts.basicOpt, opts.minGraphSize, opts.threshold, opts.C_thresh, curThread,threadsOpt);
            }
            //Increment thread and revert back to original graph
            if (curThread < nT) {
                //Skip copying at the very end
                //Old graph is already destroyed in the above function
                G = (graph *) malloc (sizeof(graph)); //Allocate new space
                duplicateGivenGraph(G_original, G); //Copy the original graph to G
            }
            curThread = curThread*2; //Increment by powers of two
        }//End of while()
    } else { //No strong scaling -- run once with max threads

#pragma omp parallel for
        for (long i=0; i<NV; i++) {
            C_orig[i] = -1;
        }
        if(opts.coloring != 0){
            runMultiPhaseColoring(G, C_orig, opts.coloring, opts.numColors, replaceMap, opts.minGraphSize, opts.threshold, opts.C_thresh, nT, threadsOpt);
        }else if(opts.syncType != 0){
            runMultiPhaseSyncType(G, C_orig, opts.syncType, opts.minGraphSize, opts.threshold, opts.C_thresh, nT,threadsOpt);
        }else{
            runMultiPhaseBasic(G, C_orig, opts.basicOpt, opts.minGraphSize, opts.threshold, opts.C_thresh, nT,threadsOpt);
        }
    }

    // Return C_orig as python list
    py::list C_orig_list;
    for (long i=0; i<NV; i++) {
        C_orig_list.append(C_orig[i]);
    }


    if(C_orig != 0) free(C_orig);
    //Do not free G here -- it will be done in another routine.
    
    return C_orig_list;
}

// Create python module
PYBIND11_MODULE(grappolo, m) {
    m.doc() = "Grappolo C++ interface"; // optional module docstring

    // Create function
    m.def("grappolo", &grappolo, "A function which calls grappolo");
}
