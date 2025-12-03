import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def process_and_plot(file_path):
    file_name = os.path.basename(file_path)
    
    # 1. Determine scaling name and test type
    scaling_name = ""
    test_title = "Scaling Test"
    
    if "weak" in file_name.lower():
        scaling_name = "weak"
        test_title = "Weak Scaling"
    elif "strong" in file_name.lower():
        scaling_name = "strong"
        test_title = "Strong Scaling"
    elif "complexity" in file_name.lower():
        scaling_name = "complexity"
        test_title = "Complexity Scaling"
    else:
        scaling_name = "unknown"

    print(f"Processing {file_name} as {test_title}...")

    # 2. Read data
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {file_name}: {e}")
        return

    cols_of_interest = ['model', 'total_time', 'thread_count', 'num_nodes']
    if not all(col in df.columns for col in cols_of_interest):
        print(f"Missing columns in {file_name}")
        return
    
    # Clean data
    df = df[cols_of_interest].dropna(subset=['total_time'])

    # 3. Decide num of subplots
    distinct_threads = sorted(df['thread_count'].unique())
    has_baseline = 1 in distinct_threads
    
    num_subplots = 3 if has_baseline else 1
    
    if num_subplots == 1:
        figsize = (10, 6)
    else:
        figsize = (10, 15)

    # 4. Group by and Average
    df_agg = df.groupby(['model', 'thread_count', 'num_nodes'], as_index=False)['total_time'].mean()
    df_agg.rename(columns={'total_time': 'runtime'}, inplace=True)

    # 5. Calculate Metrics (Speedup/Efficiency) if baseline exists
    df_plot = df_agg.copy()
    if has_baseline:
        base_df = df_agg[df_agg['thread_count'] == 1].copy()
        base_df = base_df[['model', 'num_nodes', 'runtime']]
        base_df.rename(columns={'runtime': 'runtime_base'}, inplace=True)
        
        # Merge
        if scaling_name == "weak":
             # Weak scaling: compare against base thread=1
             base_df_weak = base_df.groupby('model', as_index=False)['runtime_base'].mean()
             df_merged = pd.merge(df_agg, base_df_weak, on=['model'], how='left')
              
             # Scaled Speedup = N * (T1 / TN)
             df_merged['speedup'] = df_merged['thread_count'] * (df_merged['runtime_base'] / df_merged['runtime'])
             
             # Weak Efficiency = T1 / TN (Should be 1.0 ideally)
             df_merged['efficiency'] = df_merged['runtime_base'] / df_merged['runtime']

        else:
             # Strong Scaling: compare same model same nodes
             df_merged = pd.merge(df_agg, base_df, on=['model', 'num_nodes'], how='left')
              
             # Speedup = T1 / TN
             df_merged['speedup'] = df_merged['runtime_base'] / df_merged['runtime']
             # Efficiency = Speedup / N
             df_merged['efficiency'] = df_merged['speedup'] / df_merged['thread_count']
        
        df_plot = df_merged

    # 6. Plotting
    # Determine X-axis
    if scaling_name == "complexity":
        x_col = 'num_nodes'
        x_label = 'Number of Nodes'
    else:
        x_col = 'thread_count'
        x_label = 'Thread Count'

    models = sorted(df_plot['model'].unique())
    
    fig, axes = plt.subplots(num_subplots, 1, figsize=figsize, sharex=True)
    if num_subplots == 1:
        axes = [axes] # Ensure iterable

    # -- Subplot 1: Runtime --
    ax = axes[0]
    for model in models:
        subset = df_plot[df_plot['model'] == model].sort_values(by=x_col)
        label_str = f"Model {model}"
        ax.plot(subset[x_col], subset['runtime'], marker='o', label=label_str)
    
    ax.set_ylabel('Runtime (s)')
    ax.set_title(f'{test_title}: Runtime vs {x_label}')
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.6)
    
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)

    if num_subplots == 3:
        # -- Subplot 2: Speedup --
        ax = axes[1]
        for model in models:
            subset = df_plot[df_plot['model'] == model].sort_values(by=x_col)
            label_str = f"Model {model}"
            ax.plot(subset[x_col], subset['speedup'], marker='o', linewidth=1.5, label=label_str)
        
        # Update label to indicate Scaled Speedup for Weak scaling
        y_label_speedup = 'Scaled Speedup' if scaling_name == 'weak' else 'Speedup'
        ax.set_ylabel(y_label_speedup)
        ax.set_title(f'{test_title}: {y_label_speedup} vs {x_label}')
        
        # Add ideal linear line for Weak Scaling Speedup
        if scaling_name == "weak":
            # Ideal weak speedup is linear (y=x) because work increases with threads
            ideal_x = sorted(df_plot[x_col].unique())
            ax.plot(ideal_x, ideal_x, 'k--', alpha=0.5, label='Ideal Linear')

        ax.legend()
        ax.grid(True, which="both", ls="--", alpha=0.6)
        ax.set_xscale('log', base=2)
        

        # -- Subplot 3: Efficiency --
        ax = axes[2]
        for model in models:
            subset = df_plot[df_plot['model'] == model].sort_values(by=x_col)
            label_str = f"Model {model}"
            ax.plot(subset[x_col], subset['efficiency'], marker='o', label=label_str)
        ax.set_ylabel('Efficiency')
        ax.set_title(f'{test_title}: Efficiency vs {x_label}')
         
        ax.set_ylim(0, 1.1) 
        
        ax.legend()
        ax.grid(True, which="both", ls="-", alpha=0.6)
        ax.set_xscale('log', base=2)

    plt.xlabel(x_label)
    plt.tight_layout()
    
    output_file = f'plot_{scaling_name}.png'
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Generated {output_file}")

# Usage
if __name__ == "__main__":
    files = ['weak_scaling.csv', 'strong_scaling.csv', 'complexity_scaling.csv']
    for f in files:
        if os.path.exists(f):
            process_and_plot(f)