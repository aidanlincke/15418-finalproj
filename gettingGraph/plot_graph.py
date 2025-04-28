import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

def load_node_coordinates(filename="node_coordinates.txt"):
    """
    Load node coordinates from file
    """
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        print("Make sure to run convert_to_ch.py first to generate coordinates.")
        sys.exit(1)
        
    nodes = {}  # Map from node id to (lat, lon)
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                node_id = int(parts[0])
                lat = float(parts[1])
                lon = float(parts[2])
                nodes[node_id] = (lat, lon)
    
    print(f"Loaded {len(nodes)} node coordinates.")
    return nodes

def load_edges(filename="pittsburgh_graph.txt"):
    """
    Load graph edges from file
    """
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        print("Make sure to run convert_to_ch.py first to generate the graph.")
        sys.exit(1)
        
    edges = []  # List of (from_id, to_id, weight)
    
    with open(filename, 'r') as f:
        # Skip first line (metadata)
        header = f.readline()
        num_nodes, num_edges = map(int, header.strip().split())
        
        # Read all edges
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                from_id = int(parts[0])
                to_id = int(parts[1])
                weight = float(parts[2])
                edges.append((from_id, to_id, weight))
    
    print(f"Loaded {len(edges)} edges.")
    return edges

def plot_graph(node_coords, edges, output_file=None, 
               title="Pittsburgh Road Network", 
               plot_nodes=True,
               color_by_weight=True,
               edge_alpha=0.4,
               node_alpha=0.3,
               node_size=2):
    """
    Plot the graph with edges and optionally nodes
    """
    # Create figure and axis
    plt.figure(figsize=(12, 10))
    
    # Extract all latitudes and longitudes for bounding box
    lats = [lat for lat, lon in node_coords.values()]
    lons = [lon for lat, lon in node_coords.values()]
    
    # Calculate bounding box with some padding
    lat_min, lat_max = min(lats), max(lats)
    lon_min, lon_max = min(lons), max(lons)
    
    # Add padding (10%)
    lat_padding = (lat_max - lat_min) * 0.1
    lon_padding = (lon_max - lon_min) * 0.1
    
    # Set axis limits
    plt.xlim(lon_min - lon_padding, lon_max + lon_padding)
    plt.ylim(lat_min - lat_padding, lat_max + lat_padding)
    
    # Create a colormap from blue to red
    cmap = plt.cm.coolwarm
    
    # Plot edges
    if color_by_weight:
        # Get edge weights for coloring
        weights = [weight for _, _, weight in edges]
        min_weight, max_weight = min(weights), max(weights)
        
        # Normalize weights to [0, 1] range for coloring
        if min_weight != max_weight:
            normalized_weights = [(w - min_weight) / (max_weight - min_weight) for w in weights]
        else:
            normalized_weights = [0.5] * len(weights)
        
        # Create line segments for edges
        lines = []
        colors = []
        for (from_id, to_id, weight), norm_weight in zip(edges, normalized_weights):
            if from_id in node_coords and to_id in node_coords:
                from_lat, from_lon = node_coords[from_id]
                to_lat, to_lon = node_coords[to_id]
                lines.append([(from_lon, from_lat), (to_lon, to_lat)])
                colors.append(norm_weight)
        
        # Create line collection with proper color array
        lc = LineCollection(lines, cmap=cmap, alpha=edge_alpha, linewidths=0.5)
        lc.set_array(np.array(colors))  # Use set_array instead of directly passing colors
        plt.gca().add_collection(lc)
        
        # Add colorbar
        plt.colorbar(lc, label="Edge Weight (meters)")
    else:
        for from_id, to_id, _ in edges:
            if from_id in node_coords and to_id in node_coords:
                from_lat, from_lon = node_coords[from_id]
                to_lat, to_lon = node_coords[to_id]
                plt.plot([from_lon, to_lon], [from_lat, to_lat], 'b-', alpha=edge_alpha, linewidth=0.5)
    
    # Plot nodes if requested
    if plot_nodes:
        node_lats = [lat for lat, lon in node_coords.values()]
        node_lons = [lon for lat, lon in node_coords.values()]
        plt.scatter(node_lons, node_lats, s=node_size, color='blue', alpha=node_alpha)
    
    # Add title and labels
    plt.title(title)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    
    # Add some explanatory text
    plt.text(lon_min, lat_min - lat_padding * 0.9, 
             f"Graph with {len(node_coords)} nodes and {len(edges)} edges",
             fontsize=10, ha='left')
    
    # Add background map (if cartopy is available)
    try:
        import cartopy.crs as ccrs
        from cartopy.io.img_tiles import Stamen
        
        # Replace the current axes with cartopy axes
        plt.close()
        
        # Create a new figure with cartopy axes
        fig = plt.figure(figsize=(12, 10))
        
        # Create a Stamen Terrain background
        tiles = Stamen('terrain-background')
        ax = fig.add_subplot(1, 1, 1, projection=tiles.crs)
        
        # Set the extent of the map
        ax.set_extent([lon_min - lon_padding, lon_max + lon_padding, 
                      lat_min - lat_padding, lat_max + lat_padding], 
                      crs=ccrs.PlateCarree())
        
        # Add the background image
        ax.add_image(tiles, 10)
        
        # Plot edges with cartopy transformation
        if color_by_weight:
            # Create line segments for edges
            lines = []
            for (from_id, to_id, _) in edges:
                if from_id in node_coords and to_id in node_coords:
                    from_lat, from_lon = node_coords[from_id]
                    to_lat, to_lon = node_coords[to_id]
                    lines.append([(from_lon, from_lat), (to_lon, to_lat)])
            
            # Create line collection with proper color array
            lc = LineCollection(lines, cmap=cmap, alpha=edge_alpha, 
                              linewidths=0.5, transform=ccrs.PlateCarree())
            lc.set_array(np.array(colors))  # Use set_array instead of directly passing colors
            ax.add_collection(lc)
            
            # Add colorbar
            plt.colorbar(lc, label="Edge Weight (meters)", shrink=0.7)
        else:
            for from_id, to_id, _ in edges:
                if from_id in node_coords and to_id in node_coords:
                    from_lat, from_lon = node_coords[from_id]
                    to_lat, to_lon = node_coords[to_id]
                    ax.plot([from_lon, to_lon], [from_lat, to_lat], 'b-', 
                           alpha=edge_alpha, linewidth=0.5, transform=ccrs.PlateCarree())
        
        # Plot nodes if requested
        if plot_nodes:
            node_lats = [lat for lat, lon in node_coords.values()]
            node_lons = [lon for lat, lon in node_coords.values()]
            ax.scatter(node_lons, node_lats, s=node_size, color='blue', 
                      alpha=node_alpha, transform=ccrs.PlateCarree())
        
        # Add title
        plt.title(title)
        
        # Add some explanatory text
        ax.text(lon_min, lat_min - lat_padding * 0.5, 
                f"Graph with {len(node_coords)} nodes and {len(edges)} edges",
                fontsize=10, ha='left', transform=ccrs.PlateCarree())
        
        print("Added map background using cartopy.")
    except ImportError:
        print("Cartopy not available. Using basic plot without map background.")
        print("To add a map background, install cartopy with: pip install cartopy")
    
    # Tight layout
    plt.tight_layout()
    
    # Save or display
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

def main():
    # Check for arguments (which file to use)
    graph_file = "pittsburgh_graph.txt"  # Default to the smaller graph
    if len(sys.argv) > 1 and sys.argv[1] == "full":
        graph_file = "pittsburgh_graph.txt"
        print("Using full Pittsburgh graph.")
    
    # Load data
    nodes = load_node_coordinates("node_coordinates.txt")
    edges = load_edges(graph_file)
    
    # Plot graph (with and without nodes, for comparison)
    plot_graph(nodes, edges, output_file="pittsburgh_road_network.png")
    plot_graph(nodes, edges, output_file="pittsburgh_road_network_edges_only.png", 
               plot_nodes=False, edge_alpha=0.6)
    
    print("Plots generated. Check pittsburgh_road_network.png and pittsburgh_road_network_edges_only.png")

if __name__ == "__main__":
    main() 