import networkx as nx
import grappolo

# Get karate graph
Karate = nx.karate_club_graph()
grappolo_args = {
    "ftype": 7,
    "strongScaling": False,
    "output": False,
    "VF": False,
    "coloring": 0,
    "numColors": 16,
    "percentage": 100,
    "syncType": 0,
    "threadsOpt": False,
    "basicOpt": 0,
    "C_thresh": 0.01,
    "minGraphSize": 100000,
    "threshold": 0.000001,
    "compute_metrics": False,
    "numVertices" : Karate.number_of_nodes(),
    "numEdges" : Karate.size(),
    "inFile" : "",
    "replaceMap" : False,
    "edgeList" : nx.to_edgelist(Karate)
}
Grappolo_edgelist = grappolo.grappolo(grappolo_args)
print(Grappolo_edgelist)
