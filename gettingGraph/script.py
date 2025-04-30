import osmnx as ox
import networkx as nx

place_name = "Pittsburgh, Pennsylvania, USA"
G = ox.graph_from_place(place_name, network_type="drive")
G = ox.distance.add_edge_lengths(G)

# Print custom graph summary
print(f"Graph type: {type(G)}")
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")

for u in G.nodes():
    print(f"Node {u} has neighbors: {list(G.neighbors(u))}")

# Print some edges with distance
for u, v, data in list(G.edges(data=True))[:5]:
    print(f"From {u} to {v} — {data['length']:.1f} meters")

ox.plot_graph(G)
