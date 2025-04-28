#include <vector>
#include <utility>
#include <random>
#include <algorithm>
//#include <omp.h>
#include <iostream>
#include <unordered_map>

// Represents a graph using an adjacency list
// Each vertex is represented by an index, and its neighbors are stored in a vector
// The pair<int, bool> represents (neighbor_vertex, flip)
using Graph = std::vector<std::unordered_map<int, bool>>;

// Function to process a graph represented as an adjacency list
// Parameters:
//   graph: The input graph represented as an adjacency list
//   num_vertices: The number of vertices in the graph (defaults to graph.size())
// Returns:
//   A processed version of the graph (specific return type to be determined based on use case)
Graph processGraph(const Graph& graph, int num_vertices = -1) {
    if (num_vertices == -1) {
        num_vertices = graph.size();
    }
    // Implementation to be added
    return graph;
}

// Sets the number of OpenMP threads to use
// Parameters:
//   num_threads: The number of threads to use (defaults to the number of available cores)
/*void setNumThreads(int num_threads = 0) {
    if (num_threads <= 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
}*/

// Helper function to print a graph
void printGraph(const Graph& graph) {
    for (int v = 0; v < graph.size(); ++v) {
        std::cout << "Vertex " << v << ": ";
        for (const auto& [neighbor, flip] : graph[v]) {
            std::cout << "(" << neighbor << ", " << (flip ? "true" : "false") << ") ";
        }
        std::cout << std::endl;
    }
}

// Performs one round of edge contraction on the graph
// Parameters:
//   graph: The input graph to contract
// Returns:
//   A new graph after one round of edge contraction
//   The mapping from new vertex indices to original vertex indices
std::pair<Graph, std::vector<int>> contractGraph(const Graph& graph) {
    Graph result = graph;
    
    //#pragma omp parallel
    //{
    // Each thread gets its own RNG with a different seed
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<bool> coin(false, true);
    
    //#pragma omp for
    for (int v = 0; v < result.size(); ++v) {
        // For each neighbor
        for (auto& [neighbor, flip] : result[v]) {
            if (neighbor > v) {
                flip = coin(rng);
                result[neighbor][v] = flip;
            }
        }
    }
    //}
    printGraph(result);

    // Second phase: identify and contract matching pairs
    // initially empty int vector:
    std::vector<int> renaming;
    int current_index = 0;

    // mapping from vertex to 0, 1, 2. 0 is default value
    std::vector<int> vertex_to_label(result.size(), 0);
    
    for (int v = 0; v < result.size(); ++v) {
        if (vertex_to_label[v] == 0) {
            // haven't seen this vertex before
            renaming.push_back(current_index);

            int true_index = -1;
            bool found_two = false;
            for (auto& [neighbor, flip] : result[v]) {
                if (flip) {
                    if (true_index == -1) {
                        true_index = neighbor;
                    }
                    else {
                        found_two = true;
                    }
                }
            }
            if (true_index != -1 && !found_two && vertex_to_label[true_index] == 0) {
                // need to verify that true_index has only 1 true neighbor
                for (auto& [neighbor, flip] : result[true_index]) {
                    if (flip && neighbor != v) {
                        found_two = true;
                    }
                }
                if (!found_two) {
                    renaming[true_index] = 2; // true_index gets contracted away

                    // add all of true_index's neighbors to v's neighbors
                    for (auto& [neighbor, flip] : result[true_index]) {
                        result[v][neighbor] = flip;
                    }
                } else {
                    renaming[true_index] = 1; // true_index cannot be contracted
                }
            } else {
                vertex_to_label[v] = 1;
            }
        } else if (vertex_to_label[v] == 1) {
            renaming.push_back(current_index);
        }
        current_index++;
    }

    return {result, renaming};
}

int main() {
    // Set number of threads
    //setNumThreads(1);

    // Create a simple example graph
    // Graph with 4 vertices in a cycle: 0-1-2-3-0
    // Create a graph with 4 vertices (0 to 3)
    Graph graph(4);

    // Add edges to make a cycle: 0-1-2-3-0
    graph[0][1] = false;
    graph[0][3] = false;

    graph[1][0] = false;
    graph[1][2] = false;

    graph[2][1] = false;
    graph[2][3] = false;

    graph[3][0] = false;
    graph[3][2] = false;

    std::cout << "Original graph:" << std::endl;
    printGraph(graph);

    // Perform edge contraction
    auto [contracted_graph, renaming] = contractGraph(graph);

    std::cout << "\nContracted graph:" << std::endl;
    printGraph(contracted_graph);

    std::cout << "\nMapping:" << std::endl;
    for (int i = 0; i < renaming.size(); ++i) {
        std::cout << "Vertex " << i << " -> " << renaming[i] << std::endl;
    }

    return 0;
}