import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_allocation_trajectory(file_path):
    """
    Reads a simulation CSV and produces a stacked area plot of node allocations.
    Style matches previous scripts but with the grid removed.
    """
    # 1. Parse filename for title
    base_name = os.path.basename(file_path)
    # create a readable title from filename
    title_name = base_name.replace('.csv', '').replace('adjacency_list_', '').replace('_', ' ').title()
    
    print(f"Processing {base_name}...")

    # 2. Read Data
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return

    # 3. Extract Data
    # Filter columns that end with '_alloc'
    alloc_cols = [col for col in df.columns if col.endswith('_alloc')]
    
    if not alloc_cols:
        print(f"Warning: No allocation columns ('_alloc') found in {base_name}.")
        return

    # Clean labels (e.g., 'n52_alloc' -> 'Node 52')
    labels = [col.replace('_alloc', '').replace('n', 'Node ') for col in alloc_cols]
    
    # Prepare X and Y data
    steps = df['step']
    alloc_data = df[alloc_cols].T

    # 4. Plotting Configuration
    fig, ax = plt.subplots(figsize=(12, 8))

    # Color palette
    colors = plt.cm.tab10.colors
    if len(alloc_cols) > len(colors):
        import itertools
        colors = list(itertools.islice(itertools.cycle(colors), len(alloc_cols)))
    else:
        colors = colors[:len(alloc_cols)]

    # Stackplot
    ax.stackplot(steps, alloc_data, labels=labels, colors=colors, alpha=1, edgecolor='white', linewidth=0.5)

    # 5. Styling
    # Title and Labels
    ax.set_title(f'Trajectory of Allocations: {title_name}', fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Step', fontsize=12)
    ax.set_ylabel('Cumulative Allocation', fontsize=12)
    
    # Grid: Removed as requested ("cancel the grid inner plot")
    ax.grid(False)
    
    # Limits
    ax.set_xlim(left=0, right=steps.max())
    ax.set_ylim(bottom=0)

    # Legend
    handles, plt_labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], plt_labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), 
              borderaxespad=0, title="Top Spillovers")

    # 6. Save Output
    plt.tight_layout()
    output_filename = f"plot_{base_name.replace('.csv', '.png')}"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {output_filename}")
    
    plt.close(fig)

# --- Main Execution ---
if __name__ == "__main__":
    input_files = [
        "adj_list_global_trace_0.csv",
        "adj_list_local_trace_0.csv"
    ]
    
    for csv_file in input_files:
        plot_allocation_trajectory(csv_file)