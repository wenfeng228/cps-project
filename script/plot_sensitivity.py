import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def plot_sensitivity(file_path):
    # 1. Parse input data file name to get the parameter
    filename = os.path.basename(file_path)
    # Assumes format "sensitivity_{parameter}.csv"
    try:
        parameter = filename.split('_')[1].split('.')[0]
    except IndexError:
        print("Error: Filename must be in format 'sensitivity_{parameter}.csv'")
        return

    # 2. Read the data
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return

    # Check if parameter column exists
    if parameter not in df.columns:
        print(f"Error: Column '{parameter}' not found in the CSV.")
        return

    # Get distinct parameter values
    distinct_values = sorted(df[parameter].unique())
    
    # 6. Configure subplots (2,2)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    # 4. Generate subplot for each distinct parameter value
    for i, val in enumerate(distinct_values):
        if i >= 4:
            print("Warning: More than 4 distinct parameter values found. Only first 4 plotted.")
            break
            
        ax = axes[i]
        
        # Filter data by parameter value
        subset = df[df[parameter] == val]
        
        # Group by model to plot separate lines
        models = sorted(subset['model'].unique())
        
        for model in models:
            model_data = subset[subset['model'] == model]
            # 5. Plot: x=step, y=l2, label=model
            ax.plot(model_data['step'], model_data['l2'], label=f'Model {model}', linewidth=1.5)
        
        ax.axhline(y=1e-4, color='red', linestyle='--', linewidth=1.2, label='Tolerance ($10^{-4}$)')

        # Subplot configuration
        ax.set_title(f'{parameter} = {val}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Step')
        ax.set_ylabel('L2 Loss')
        ax.legend(title='Model')
        ax.grid(True, linestyle='--', alpha=0.6)

    # Clean up empty subplots if fewer than 4 values
    for j in range(len(distinct_values), 4):
        fig.delaxes(axes[j])

    plt.suptitle(f'Sensitivity Analysis: {parameter.capitalize()}', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to make room for suptitle

    # 7. Output the whole single plot
    output_filename = f"plot_{parameter}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Successfully generated '{output_filename}'")

if __name__ == "__main__":
    files_to_run = ["sensitivity_budget.csv", "sensitivity_gamma.csv", "sensitivity_lambda.csv"]
    for filename in files_to_run:
        if os.path.exists(filename):
            plot_sensitivity(filename)