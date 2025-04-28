import osmnx as ox
import networkx as nx
import numpy as np
import json
import os

def load_graph(place_name="Pittsburgh, Pennsylvania, USA", network_type="drive"):
    """
    Load a road network from OSMnx
    """
    print(f"Loading graph for {place_name}...")
    G = ox.graph_from_place(place_name, network_type=network_type)
    
    # Add edge lengths only
    G = ox.distance.add_edge_lengths(G)
    
    print(f"Loaded graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G

def convert_to_contraction_format(G):
    """
    Convert OSMnx graph to format for contraction hierarchies
    Returns:
        nodes_map: dict mapping OSM node IDs to sequential integers
        edges: list of (from, to, weight) tuples
    """
    # Create a mapping from OSM node IDs to sequential integers
    nodes = list(G.nodes())
    nodes_map = {node_id: i for i, node_id in enumerate(nodes)}
    
    # Extract edges with weights (using distance in meters)
    edges = []
    for u, v, k, data in G.edges(data=True, keys=True):
        # Use length as weight
        if 'length' in data:
            weight = data['length']
        else:
            # Skip edges without length data
            continue
        
        # Add edge to list
        from_idx = nodes_map[u]
        to_idx = nodes_map[v]
        edges.append((from_idx, to_idx, weight))
    
    return nodes_map, edges

def write_to_file(nodes_map, edges, filename="pittsburgh_graph.txt"):
    """
    Write graph data to a text file in a format suitable for C++ input
    """
    with open(filename, 'w') as f:
        # Write number of nodes and edges
        f.write(f"{len(nodes_map)} {len(edges)}\n")
        
        # Write edges (from, to, weight)
        for from_node, to_node, weight in edges:
            f.write(f"{from_node} {to_node} {weight:.2f}\n")
    
    print(f"Graph written to {filename}")

def write_node_coordinates(G, nodes_map, filename="node_coordinates.txt"):
    """
    Write node coordinates to file for visualization
    """
    with open(filename, 'w') as f:
        for osm_id, idx in nodes_map.items():
            lon, lat = G.nodes[osm_id]['x'], G.nodes[osm_id]['y']
            f.write(f"{idx} {lat} {lon}\n")
    
    print(f"Node coordinates written to {filename}")

def create_c_plus_plus_input(place_name="Pittsburgh, Pennsylvania, USA", network_type="drive"):
    """
    Main function to generate input file for C++ contraction hierarchies
    """
    # Load graph
    G = load_graph(place_name, network_type)
    
    # Convert to contraction format
    nodes_map, edges = convert_to_contraction_format(G)
    
    # Write to files
    write_to_file(nodes_map, edges)
    write_node_coordinates(G, nodes_map)
    
    # Extract a small subgraph for testing
    print("\nCreating a smaller test graph...")
    test_nodes = 1000
    test_edges = [(u, v, w) for u, v, w in edges if u < test_nodes and v < test_nodes]
    write_to_file(
        {k: v for k, v in nodes_map.items() if v < test_nodes}, 
        test_edges, 
        "pittsburgh_small.txt"
    )
    
    return len(nodes_map), len(edges)

def generate_cpp_loader():
    """
    Generate C++ code to load the graph from the file
    """
    cpp_code = """
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
    """
    
    with open("graph_loader.cpp", "w") as f:
        f.write(cpp_code)
    
    print("C++ loader code written to graph_loader.cpp")

if __name__ == "__main__":
    try:
        num_nodes, num_edges = create_c_plus_plus_input()
        generate_cpp_loader()
        
        print(f"\nSummary:")
        print(f"Full graph: {num_nodes} nodes, {num_edges} edges")
        print(f"Graph files created: pittsburgh_graph.txt, pittsburgh_small.txt, node_coordinates.txt")
        print(f"C++ loader code: graph_loader.cpp")
        print(f"\nTo use in your contraction hierarchies code:")
        print(f"1. Copy graph_loader.cpp contents to your project")
        print(f"2. Load the graph with: SimpleGraph graph = loadGraphFromFile(\"pittsburgh_small.txt\");")
        print(f"3. Run makeContraction() on the loaded graph")
    except Exception as e:
        print(f"Error: {e}")
        print("\nIf you're having trouble with OSMnx, you might need to install or update it:")
        print("pip install -U osmnx") 