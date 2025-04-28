
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

// Import your SimpleGraph definition here
// using SimpleGraph = std::vector<std::unordered_map<int, Edge>>;

SimpleGraph loadGraphFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return SimpleGraph();
    }
    
    int numNodes, numEdges;
    file >> numNodes >> numEdges;
    
    SimpleGraph graph(numNodes);
    
    for (int i = 0; i < numEdges; ++i) {
        int from, to;
        float weight;
        file >> from >> to >> weight;
        
        // Add edge in both directions (for undirected graph)
        addEdge(graph, from, to, weight);
    }
    
    std::cout << "Loaded graph with " << numNodes << " nodes and " 
              << numEdges << " edges." << std::endl;
    
    return graph;
}

// Usage example in main:
// SimpleGraph graph = loadGraphFromFile("pittsburgh_graph.txt");
// or for testing:
// SimpleGraph graph = loadGraphFromFile("pittsburgh_small.txt");
    