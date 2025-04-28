#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <unordered_set>
#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <omp.h>

struct Node {
    int value;
    Node* next;

    Node(int v, Node* n = nullptr) : value(v), next(n) {}
};

struct LinkedList {
    Node* start;
    Node* end;

    LinkedList() : start(nullptr), end(nullptr) {}

    void push_back(int value) {
        Node* new_node = new Node(value);
        if (start == nullptr) {
            start = new_node;
            end = new_node;
        } else {
            end->next = new_node;
            end = new_node;
        }
    }

    void append(LinkedList* other) {
        //for (Node* current = other->start; current != nullptr; current = current->next) {
        //    push_back(current->value);
        //}

        if (other->start == nullptr) {
            return;
        }
        
        if (start == nullptr) {
            start = other->start;
            end = other->end;
        } else {
            end->next = other->start;
            end = other->end;
        }
    }

    void print() const {
        Node* current = start;
        while (current != nullptr) {
            std::cout << current->value << " ";
            current = current->next;
        }
    }
};

// A custom structure to represent a weighted edge
struct Edge {
    float weight;

    // Constructor for convenience
    Edge(float w = 0.0f) : weight(w) {}
};

// A custom structure to represent a weighted edge with a path (linked list)
struct ShortcutEdge {
    float weight;
    LinkedList* path;

    // Constructor for convenience
    ShortcutEdge(float w = 0.0f) : weight(w), path(new LinkedList()) {}
    ShortcutEdge(float w, LinkedList *p) : weight(w), path(p) {}
};

// Type 1: Simple Weighted Graph
// Each vertex is represented by an index, and its neighbors are stored in an unordered_map
// The unordered_map maps from neighbor vertex (int) to edge weight (float)
using SimpleGraph = std::vector<std::unordered_map<int, Edge>>;

// Type 2: Contracted Graph
// Contains both the graph data and a set of active nodes
struct ContractedGraph {
    std::vector<std::unordered_map<int, ShortcutEdge>> edges;
    std::unordered_set<int> activeNodes;

    // Constructor to create from a SimpleGraph
    explicit ContractedGraph(const SimpleGraph& simple) {
        // Convert simple edges to shortcut edges
        edges.resize(simple.size());
        for (int v = 0; v < simple.size(); ++v) {
            for (const auto& [neighbor, edgeData] : simple[v]) {
                LinkedList *path = new LinkedList();
                edges[v][neighbor] = ShortcutEdge(edgeData.weight, path);
            }
            // Initially, all nodes are active
            activeNodes.insert(v);
        }
    }
    
    // Constructor with specified size
    explicit ContractedGraph(int size) {
        edges.resize(size);
        for (int i = 0; i < size; ++i) {
            activeNodes.insert(i);
        }
    }
};

// Function to load a SimpleGraph from a Pittsburgh graph file
SimpleGraph loadGraphFromFile(const std::string& filename) {
    // Open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return SimpleGraph();
    }
    
    // Read the number of nodes and edges
    int numNodes, numEdges;
    file >> numNodes >> numEdges;
    
    // Create graph with appropriate size
    SimpleGraph graph(numNodes);
    
    // Read each edge
    int numSelfLoops = 0;
    for (int i = 0; i < numEdges; ++i) {
        int from, to;
        float weight;
        
        file >> from >> to >> weight;

        if (from == to) {
            numSelfLoops++;
            continue;
        }
        
        // Add the edge to the graph (both directions for undirected graph)
        graph[from][to] = Edge(weight);
        graph[to][from] = Edge(weight);
    }
    numEdges -= numSelfLoops;
    
    std::cout << "Loaded graph with " << numNodes << " nodes and " << numEdges << " edges" << std::endl;
    return graph;
}

// Helper function to create a simple graph with a specified number of vertices
SimpleGraph createSimpleGraph(int numVertices) {
    return SimpleGraph(numVertices);
}

// Helper function to add an edge to the simple graph
void addEdge(SimpleGraph& graph, int from, int to, float weight) {
    graph[from][to] = Edge(weight);
    graph[to][from] = Edge(weight);
}

// Helper function to create a grid graph with rows*cols vertices
SimpleGraph createGridGraph(int rows, int cols) {
    int numVertices = rows * cols;
    SimpleGraph graph = createSimpleGraph(numVertices);
    
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> dist(1.0f, 10.0f);
    
    // Helper function to convert (row, col) to vertex index
    auto vertexIndex = [cols](int row, int col) {
        return row * cols + col;
    };
    
    // Add edges for the grid
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int current = vertexIndex(row, col);
            
            // Add edge to right neighbor
            if (col + 1 < cols) {
                int right = vertexIndex(row, col + 1);
                float weight = dist(rng);
                addEdge(graph, current, right, weight);
            }
            
            // Add edge to bottom neighbor
            if (row + 1 < rows) {
                int bottom = vertexIndex(row + 1, col);
                float weight = dist(rng);
                addEdge(graph, current, bottom, weight);
            }
            
            // Add some diagonal edges for more interesting shortcuts
            if (row + 1 < rows && col + 1 < cols && dist(rng) > 7.0f) {
                int diagonal = vertexIndex(row + 1, col + 1);
                float weight = dist(rng);
                addEdge(graph, current, diagonal, weight);
            }
        }
    }
    
    return graph;
}


// Helper function to contract a node in the contracted graph
void contractNode(ContractedGraph& graph, int nodeToContract) {    
    for (const auto& [start, edge_in] : graph.edges[nodeToContract]) {

        for (const auto& [end, edge_out] : graph.edges[nodeToContract]) {
            
            if (start < end) {
                float new_weight = edge_in.weight + edge_out.weight;
                
                if (graph.edges[start].find(end) == graph.edges[start].end() || graph.edges[start][end].weight > new_weight) {
                    LinkedList *new_path = new LinkedList();

                    new_path->append(graph.edges[start][nodeToContract].path);
                    new_path->push_back(nodeToContract);
                    new_path->append(edge_out.path);

                    graph.edges[start][end] = ShortcutEdge(new_weight, new_path);

                    LinkedList *new_path_reverse = new LinkedList();
                    new_path_reverse->append(graph.edges[end][nodeToContract].path);
                    new_path_reverse->push_back(nodeToContract);
                    new_path_reverse->append(edge_in.path);

                    graph.edges[end][start] = ShortcutEdge(new_weight, new_path_reverse);

                }
            }
        }
    }

    for (const auto& [neighbor, edge] : graph.edges[nodeToContract]) {
        graph.edges[neighbor].erase(nodeToContract);
    }

    graph.activeNodes.erase(nodeToContract);
}

ContractedGraph makeContraction(SimpleGraph& graph) {
    ContractedGraph contractedGraph(graph);
    
    // Contract selected nodes in parallel
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < graph.size(); ++i) {
        if (i % 11 != 0) {  // Your existing selection logic
            // Each transaction handles one node contraction
            #pragma omp transaction
            {
                contractNode(contractedGraph, i);
            }
        }
    }

    return contractedGraph;
}

// Helper function to print a simple graph
void printSimpleGraph(const SimpleGraph& graph) {
    for (int v = 0; v < graph.size(); ++v) {
        std::cout << "Vertex " << v << ":\n";
        for (const auto& [neighbor, edge] : graph[v]) {
            std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ")\n";
        }
    }
}

// Helper function to print a contracted graph
void printContractedGraph(const ContractedGraph& graph) {
    std::cout << "Active nodes: ";
    for (int node : graph.activeNodes) {
        std::cout << node << " ";
    }
    std::cout << "\n";

    for (int node : graph.activeNodes) {
        std::cout << "Vertex " << node << ":\n";
        for (const auto& [neighbor, edge] : graph.edges[node]) {
            std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ", path: ";
            edge.path->print();
            std::cout << ")\n";
        }
    }
    /*for (int v = 0; v < graph.edges.size(); ++v) {
        std::cout << "Vertex " << v;
        if (graph.activeNodes.count(v) == 0) {
            std::cout << " (contracted)";
        }
        std::cout << ":\n";
        
        for (const auto& [neighbor, edge] : graph.edges[v]) {
            std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ", path: ";
            edge.path->print();
            std::cout << ")\n";
        }
    }*/
}

bool verifyGraph(SimpleGraph& graph) {
    int numNodes = graph.size();
    for (int i = 0; i < numNodes; ++i) {
        for (const auto& [neighbor, edge] : graph[i]) {
            if (neighbor == i || neighbor < 0 || neighbor >= numNodes) {
                std::cout << "Invalid neighbor: " << neighbor << std::endl;
                return false;
            }
            if (graph[neighbor].find(i) == graph[neighbor].end() || graph[neighbor][i].weight != edge.weight) {
                std::cout << "Invalid edge weight: " << graph[neighbor][i].weight << " != " << edge.weight << std::endl;
                return false;
            }
        }
    }
    return true;
}


// Example usage
int main() {
    // Test with small graph
    /*
    std::cout << "=== Testing with small graph ===\n";
    SimpleGraph simpleGraph = createSimpleGraph(5);
    
    // Add some edges
    addEdge(simpleGraph, 0, 1, 4.5);
    addEdge(simpleGraph, 0, 2, 3.1);
    addEdge(simpleGraph, 1, 2, 2.0);
    addEdge(simpleGraph, 1, 3, 1.7);
    addEdge(simpleGraph, 2, 3, 2.5);
    addEdge(simpleGraph, 3, 4, 5.0);
    
    // Print the simple graph
    std::cout << "Simple Graph representation:\n";
    printSimpleGraph(simpleGraph);
    */

    // Test with grid graph
    std::cout << "\n=== Testing with grid graph ===\n";
    SimpleGraph gridGraph = createGridGraph(10, 10);
    printSimpleGraph(gridGraph);
    assert(verifyGraph(gridGraph));

    ContractedGraph contractedGrid = makeContraction(gridGraph);
    printContractedGraph(contractedGrid);

    // Test loading Pittsburgh graph
    std::cout << "\n=== Loading Pittsburgh graph ===\n";
    std::string filepath = "../gettingGraph/pittsburgh_graph.txt";
    SimpleGraph pittsburghGraph = loadGraphFromFile(filepath);
    
    if (!pittsburghGraph.empty()) {
        std::cout << "Successfully loaded Pittsburgh graph\n";
        
        assert(verifyGraph(pittsburghGraph));
        
        // Test contraction on Pittsburgh graph
        std::cout << "\nStarting contraction of Pittsburgh graph...\n";
        auto startTime = std::chrono::high_resolution_clock::now();
        
        ContractedGraph contractedPittsburgh = makeContraction(pittsburghGraph);
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        std::cout << "Contraction complete in " << duration << " ms.\n";
    }
    
    return 0;
}
