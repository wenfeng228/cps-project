import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def visualize_network_directed(nodes_file, edges_file, output_file='plot_network.png'):
    # 1. Load Data
    try:
        nodes = pd.read_csv(nodes_file)
        edges = pd.read_csv(edges_file)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # 2. Create Directed Graph
    G = nx.from_pandas_edgelist(edges, 'origin', 'destination', edge_attr=['weight'], create_using=nx.DiGraph())
    
    if 'node' in nodes.columns:
        nx.set_node_attributes(G, nodes.set_index('node').to_dict('index'))

    # 3. Layout
    pos = nx.spring_layout(G, k=0.6, iterations=60, seed=42)

    # 4. Styling Config
    NODE_COLOR = '#333333'   # dark gray
    POS_EDGE   = '#1E90FF'   # blue
    NEG_EDGE   = '#DC143C'   # red
    
    weights = [G[u][v]['weight'] for u, v in G.edges]
    widths = [abs(w) * 2.0 for w in weights]
    edge_colors = [POS_EDGE if w > 0 else NEG_EDGE for w in weights]
    
    # Node sizing logic
    if 'initial' in nodes.columns:
        initial_vals = [G.nodes[n].get('initial', 1) for n in G.nodes]
        min_v, max_v = min(initial_vals), max(initial_vals)
        node_sizes = [300 + (x - min_v) * (900) / (max_v - min_v) if max_v > min_v else 500 for x in initial_vals]
    else:
        node_sizes = 500

    # 5. Plotting
    fig, ax = plt.subplots(figsize=(14, 12)) 

    # LAYER 1: The Connections (Curved Arrows)
    nx.draw_networkx_edges(G, pos,
                           width=widths,
                           edge_color=edge_colors,
                           alpha=0.7, 
                           arrows=True,
                           arrowstyle='-|>',
                           arrowsize=15,
                           connectionstyle='arc3, rad=0.1', 
                           ax=ax)

    # LAYER 2: The Nodes (Removed Shadows for simplicity)
    nx.draw_networkx_nodes(G, pos,
                           node_size=node_sizes,
                           node_color=NODE_COLOR,
                           edgecolors='black', # Changed to black for higher contrast
                           linewidths=1.0,
                           ax=ax)

    # LAYER 3: Labels
    nx.draw_networkx_labels(G, pos,
                            font_size=9,
                            font_family='sans-serif',
                            font_color='white',
                            font_weight='bold',
                            ax=ax)

    # Legend
    blue_line = mlines.Line2D([], [], color=POS_EDGE, linewidth=3, label='Positive Flow')
    red_line = mlines.Line2D([], [], color=NEG_EDGE, linewidth=3, label='Negative Flow')
    ax.legend(handles=[blue_line, red_line], loc='upper right', 
              frameon=True, fontsize=10)
 
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    visualize_network_directed('nodes.csv', 'edges.csv')