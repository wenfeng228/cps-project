import pandas as pd
import numpy as np
import sys
import os

def expand_network(nodes_df, edges_df, iterations=10):
    """
    Expands a network by recursively duplicating nodes and edges 
    and randomly rewiring them to preserve structure and connectivity.
    """
    # Ensure nodes are 0-indexed and contiguous
    nodes_df = nodes_df.sort_values('node').reset_index(drop=True)
    
    current_nodes = nodes_df.copy()
    current_edges = edges_df.copy()
    
    print(f"Initial: {len(current_nodes)} nodes, {len(current_edges)} edges")

    for i in range(iterations):
        n_count = len(current_nodes)
        
        # 1. Duplicate Nodes (Shift IDs by current count)
        new_nodes = current_nodes.copy()
        new_nodes['node'] = new_nodes['node'] + n_count
        
        # 2. Duplicate Edges (Shift IDs)
        new_edges = current_edges.copy()
        new_edges['origin'] = new_edges['origin'] + n_count
        new_edges['destination'] = new_edges['destination'] + n_count
        
        # 3. Mixing / Rewiring
        # With 50% probability, swap the destination of the edge 
        # between the original and the copy to connect the two components.
        mask = np.random.rand(len(current_edges)) < 0.5
        swap_indices = np.where(mask)[0]
        
        # Reset indices to ensure alignment
        current_edges = current_edges.reset_index(drop=True)
        new_edges = new_edges.reset_index(drop=True)
        
        # Perform the swap: (u, v) -> (u, v') and (u', v') -> (u', v)
        dest_old = current_edges.loc[swap_indices, 'destination'].values
        dest_new = new_edges.loc[swap_indices, 'destination'].values
        
        current_edges.loc[swap_indices, 'destination'] = dest_new
        new_edges.loc[swap_indices, 'destination'] = dest_old
        
        # 4. Merge
        current_nodes = pd.concat([current_nodes, new_nodes], ignore_index=True)
        current_edges = pd.concat([current_edges, new_edges], ignore_index=True)
        
        print(f"Iteration {i+1}: Nodes {len(current_nodes)}, Edges {len(current_edges)}")
        
    return current_nodes, current_edges

# --- Argument Control Added Below ---

if __name__ == "__main__":
    # Check if correct number of arguments are passed
    if len(sys.argv) < 4:
        print("Usage: python sampling.py <nodes_file> <edges_file> <iterations>")
        sys.exit(1)

    # 1. Parse Arguments
    node_file = sys.argv[1]
    edge_file = sys.argv[2]
    iterations = int(sys.argv[3])

    # Load original data
    nodes = pd.read_csv(node_file)
    edges = pd.read_csv(edge_file)

    # Run Expansion using the iterations argument
    sampled_nodes, sampled_edges = expand_network(nodes, edges, iterations=iterations)

    # 2. Create New Filenames
    # split "nodes.csv" into "nodes" and ".csv"
    node_name, node_ext = os.path.splitext(node_file)
    edge_name, edge_ext = os.path.splitext(edge_file)

    out_node_file = f"{node_name}_{iterations}{node_ext}"
    out_edge_file = f"{edge_name}_{iterations}{edge_ext}"

    # Save output with new names
    sampled_nodes.to_csv(out_node_file, index=False)
    sampled_edges.to_csv(out_edge_file, index=False)

    print(f"Saved to: {out_node_file} and {out_edge_file}")
    print("Expansion complete.")