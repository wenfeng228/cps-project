import pandas as pd
import os
import argparse

# helper
def get_param_name(dir_path): 
    folder = os.path.basename(os.path.normpath(dir_path))
    if folder.startswith('sensitivity_'):
        return folder.replace('sensitivity_', '')
    return folder

def parse_filename(filename, param_name): 
    parts = filename.split('_')
    model_id = None
    param_val = None

    for part in parts:
        # 1. Model ID
        if part.startswith('mod'):
            val = part.replace('mod', '')
            if val.isdigit():
                model_id = int(val)
        
        # 2. Parameter Value
        elif part.startswith(param_name):
            val = part.replace(param_name, '')
            try:
                # Handle int or float
                param_val = float(val)
                if param_val.is_integer():
                    param_val = int(param_val)
            except ValueError:
                continue

    return model_id, param_val

# mian
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='Directory containing csv files')
    args = parser.parse_args()

    # 1. Setup
    if not os.path.exists(args.input_dir):
        print(f"Error: {args.input_dir} not found.")
        return

    param_name = get_param_name(args.input_dir)
    output_filename = f"sensitivity_{param_name}.csv"
    
    print(f"Processing '{args.input_dir}'...")
    print(f"keeping columns: model, {param_name}, step, l2")

    all_data = []
    files = [f for f in os.listdir(args.input_dir) if f.endswith(".csv")]

    # 2. Process Files
    for f in files:
        # Get metadata from filename
        mod, val = parse_filename(f, param_name)

        if mod is not None and val is not None:
            path = os.path.join(args.input_dir, f)
            
            try:
                # Read CSV
                df = pd.read_csv(path)
                
                # Check if required columns exist before processing
                if 'step' in df.columns and 'l2' in df.columns: 
                    df = df[['step', 'l2']].copy()
                    
                    # Add metadata columns
                    df['model'] = mod
                    df[param_name] = val
                    
                    all_data.append(df)
                else:
                    print(f"Skipping {f}: Missing 'step' or 'l2' columns.")
                    
            except Exception as e:
                print(f"Error reading {f}: {e}")

    # 3. Save
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        
        # Reorder to: model, param, step, l2
        cols = ['model', param_name, 'step', 'l2']
        final_df = final_df[cols]
        
        # Sort
        final_df = final_df.sort_values(by=['model', param_name, 'step'])
        
        final_df.to_csv(output_filename, index=False)
        print(f"Success! Saved clean data to: {output_filename}")
    else:
        print("No valid data found.")

if __name__ == "__main__":
    main()