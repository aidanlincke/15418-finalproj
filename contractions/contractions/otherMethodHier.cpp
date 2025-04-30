// This is the "Directed graph" approach mentioned in the writeup
// Parallelize by keeping accesses to a given node's map isolated
// to just 1 thread.

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
#include <cassert>
//#include <immintrin.h>
#include <tbb/concurrent_unordered_map.h>

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
        for (Node* current = other->start; current != nullptr; current = current->next) {
            push_back(current->value);
        }
        /*
        if (other->start == nullptr) {
            return;
        }
        
        if (start == nullptr) {
            start = other->start;
            end = other->end;
        } else {
            end->next = other->start;
            end = other->end;
        }*/
    }

    void print() const {
        Node* current = start;
        while (current != nullptr) {
            std::cout << current->value << " ";
            current = current->next;
        }
    }
};

struct Edge {
    float weight;
    Edge(float w = 0.0f) : weight(w) {}
};

struct ShortcutEdge {
    float weight;
    LinkedList* path;

    ShortcutEdge(float w = 0.0f) : weight(w), path(new LinkedList()) {}

    ShortcutEdge(float w, LinkedList *p) : weight(w), path(p) {}
};

using SimpleGraph = std::vector<std::unordered_map<int, Edge>>;

struct ContractedGraph {
    std::vector<tbb::concurrent_unordered_map<int, ShortcutEdge>> edges;
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

SimpleGraph loadGraphFromFile(const std::string& filename) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return SimpleGraph();
    }
    
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


void addEdge(SimpleGraph& graph, int from, int to, float weight) {
    graph[from][to] = Edge(weight);
    graph[to][from] = Edge(weight);
}

// Helper function to contract a node in the contracted graph
void contractNode(ContractedGraph& graph, int nodeToContract) {
    std::vector<std::pair<int, ShortcutEdge>> edgeList(
        graph.edges[nodeToContract].begin(),
        graph.edges[nodeToContract].end()
    );

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < edgeList.size(); ++i) {
        const auto& [start, edge_in] = edgeList[i];
        //for (const auto& [start, edge_in] : graph.edges[nodeToContract]) {

        for (const auto& [end, edge_out] : graph.edges[nodeToContract]) {
            
            if (start != end) {
                float new_weight = edge_in.weight + edge_out.weight;
                
                if (graph.edges[start].find(end) == graph.edges[start].end() || graph.edges[start][end].weight > new_weight) {
                    LinkedList *new_path = new LinkedList();

                    new_path->append(graph.edges[start][nodeToContract].path);
                    new_path->push_back(nodeToContract);
                    new_path->append(edge_out.path);

                    graph.edges[start][end] = ShortcutEdge(new_weight, new_path);

                    // LinkedList *new_path_reverse = new LinkedList();
                    // new_path_reverse->append(graph.edges[end][nodeToContract].path);
                    // new_path_reverse->push_back(nodeToContract);
                    // new_path_reverse->append(edge_in.path);

                    // graph.edges[end][start] = ShortcutEdge(new_weight, new_path_reverse);

                }
            }
        }
    }

    for (const auto& [neighbor, edge] : graph.edges[nodeToContract]) {
        graph.edges[neighbor].erase(nodeToContract);
    }

    graph.activeNodes.erase(nodeToContract);
}


ContractedGraph makeContractionBetter(SimpleGraph& graph) {
    ContractedGraph contractedGraph(graph);

    for (int i = 0; i < graph.size(); ++i) {
        float max_edge_len = 0;
        if (graph[i].size() < 4) {
            contractNode(contractGraph, i);
        }

    }
}

ContractedGraph makeContraction(SimpleGraph& graph) {
    ContractedGraph contractedGraph(graph);

    for (int i = 0; i < graph.size(); ++i) {
        //if (contractedGraph.edges[i].size() < 8) {
        if (i % 4 !=0) {
            //std::cout << "Contracting node " << i << std::endl;
            contractNode(contractedGraph, i);
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
    // for (int v = 0; v < graph.edges.size(); ++v) {
    //     std::cout << "Vertex " << v;
    //     if (graph.activeNodes.count(v) == 0) {
    //         std::cout << " (contracted)";
    //     }
    //     std::cout << ":\n";
        
    //     for (const auto& [neighbor, edge] : graph.edges[v]) {
    //         std::cout << "  -> " << neighbor << " (weight: " << edge.weight << ", path: ";
    //         edge.path->print();
    //         std::cout << ")\n";
    //     }
    // }
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

SimpleGraph simpleGraph1() {
    SimpleGraph simpleGraph = createSimpleGraph(5);
    
    // Add some edges
    addEdge(simpleGraph, 0, 1, 4.5);
    addEdge(simpleGraph, 0, 2, 3.1);
    addEdge(simpleGraph, 1, 2, 2.0);
    addEdge(simpleGraph, 1, 3, 1.7);
    addEdge(simpleGraph, 2, 3, 2.5);
    addEdge(simpleGraph, 3, 4, 5.0);

    return simpleGraph;
}


bool verifyContraction(ContractedGraph& contraction, SimpleGraph& graph) {
    int maxNode = contraction.edges.size();
    if (maxNode != graph.size()) {
        std::cout << "Graph size missmatch" << std::endl;
        return false;
    }

    for (const auto& node : contraction.activeNodes) {
        for (const auto& [neighbor, edge] : contraction.edges[node]) {
            if (neighbor == node || neighbor < 0 || neighbor >= maxNode) {
                std::cout << "Invalid neighbor" << std::endl;
                return false;
            }

            if (contraction.edges[neighbor][node].weight != edge.weight) {
                std::cout << "Weight missmatch" << std::endl;
                return false;
            }

            float weight = 0.0f;

            int currNode = node;
            int currNeighbor;
            Node *currBlock;
            
            for (currBlock = edge.path->start; currBlock != nullptr; currBlock = currBlock->next) {
                currNeighbor = currBlock->value;

                if (graph[currNode].find(currNeighbor) == graph[currNode].end()) {
                    std::cout << "Neighbor not found in " << node << " to " << neighbor << " from " << currNode << " to " << currNeighbor << std::endl;
                    return false;
                }

                weight += graph[currNode][currNeighbor].weight;

                currNode = currNeighbor;
            }

            currNeighbor = neighbor;

            if (graph[currNode].find(currNeighbor) == graph[currNode].end()) {
                    std::cout << "Neighbor not found in " << node << " to " << neighbor << " from " << currNode << " to " << currNeighbor << std::endl;
                    return false;
                }

            weight += graph[currNode][currNeighbor].weight;

            if (abs(edge.weight-weight) > 1) {
                std::cout << "Weight and path do not match: " << weight << " " << edge.weight << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Example usage
int main() {
    omp_set_num_threads(8);
    std::cout << "Num threads: " << omp_get_num_threads() << std::endl;

    // SimpleGraph gridGraph = createGridGraph(10, 10);
    
    // printSimpleGraph(gridGraph);
    
    // ContractedGraph gridContraction = makeContraction(gridGraph);

    // printContractedGraph(gridContraction);

    // assert(verifyContraction(gridContraction, gridGraph));

    

    // Test loading Pittsburgh graph
    std::cout << "\n=== Loading Pittsburgh graph ===\n";
    std::string filepath = "../gettingGraph/pittsburgh_graph.txt";
    SimpleGraph pittsburghGraph = loadGraphFromFile(filepath);
    
    //printSimpleGraph(pittsburghGraph);

    // Test contraction on Pittsburgh graph
    std::cout << "\nStarting contraction of Pittsburgh graph...\n";
    auto startTime = std::chrono::high_resolution_clock::now();
    
    ContractedGraph contractedPittsburgh = makeContraction(pittsburghGraph);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    
    std::cout << "Contraction complete in " << duration << " ms.\n";

    std::cout << "Num remaining nodes: " << contractedPittsburgh.activeNodes.size() << std::endl;
    //printContractedGraph(contractedPittsburgh);

    assert(verifyContraction(contractedPittsburgh, pittsburghGraph));
    
    return 0;
}