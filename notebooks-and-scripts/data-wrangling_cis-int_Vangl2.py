import pandas as pd
import glob
import os
#from pathlib import Path


def process_colocalization_data(df_celsr_vang: pd.DataFrame) -> pd.DataFrame:
    """
    Process colocalization data between Celsr+ puncta which have Vangl2.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing Celsr1 data with measurements for Vang
    
    Returns
    -------
    pd.DataFrame
        Merged DataFrame with counts grouped by variant, date, and roi
    """
    
    # Define column configurations
    column_configs = {
        'celsr1_base': ['Label', 'filename', 'date', 'Replicate', 'coverslip',
                        'Celsr1_Variant', 'image_number', 'NumberOfVoxels'],
        'intensity_patterns': {
            'Vangl2': ['IntensityMean_Celsr1-Vangl2', 'IntensityMax_Celsr1-Vangl2',]
        }
    }
    
    groupby_cols = ['Celsr1_Variant', 'Replicate', 'image_number','coverslip','filename']
    
    
    # Calculate number of Celsr puncta and % of those with Vang
    df_Celsr1_N = (df_celsr_vang.groupby(groupby_cols)['IntensityMax_Celsr1-Vangl2',]
          .size().reset_index(name='Celsr1_Counts')
          )
    
    df_Celsr1_N_with_Vang = (df_celsr_vang.groupby(groupby_cols)['IntensityMax_Celsr1-Vangl2']
          .apply(lambda x: (x > 0).sum()).reset_index(name='Celsr1_with_Vangl2_Counts')
          )

    df = pd.merge(df_Celsr1_N, df_Celsr1_N_with_Vang, how="outer", on=groupby_cols)
    df["FractionOfCelsr1_with_Vangl2"] = (df["Celsr1_with_Vangl2_Counts"] / df["Celsr1_Counts"])
    
    df.to_csv(r"C:/Users/wgiang/Documents/2026-04-07_Sarah-Vang-coloc/06_collated_results/COS7-Vangl2-coloc-merged.csv", index=False)

    # Calculate statistics
    stats_summary = df.groupby(['Celsr1_Variant', 'Replicate'])['FractionOfCelsr1_with_Vangl2'].agg(['mean', 'std']).reset_index()
    stats_summary.to_csv(r"C:/Users/wgiang/Documents/2026-04-07_Sarah-Vang-coloc/06_collated_results/COS7-Vangl2-coloc-stats_summary.csv", index=False)
    print("Statistics Summary:")
    print(stats_summary)

    return df



def ensure_dataframe(obj):
    """
    Converts a Pandas Series to a DataFrame if it is a Series, 
    otherwise returns the object unchanged.
    """
    if isinstance(obj, pd.Series):
        return obj.to_frame()
    else:
        return obj
    
def process_csv_files(directory_path):
    """
    Process CSV files for Celsr1+ puncta with measurements of masked Vangl2

    Parameters:
    directory_path (str): Path to directory containing CSV files
    
    Returns:
    pd.Dataframe 
    """
    
    # Find all CSV files with C2_V prefix and ending in Celsr1-Vangl2
    csv_pattern = os.path.join(directory_path, "C2_V*_Celsr1-Vangl2.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"No CSV files found with pattern C2_V*_Celsr1-Vangl2.csv in {directory_path}")
        return {}
    
    df_all = pd.DataFrame()
    
    # Process each CSV file
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        
        try:
            
            plasmid_mapping = {
                'p01':'Celsr1-WT',
                'p02':'Celsr-WT-mCherry',
                'p04':'Vang2-GFP',
                'p08':'Fzd6-mCherry',
                'p09':'Celsr1-Crsh-GFP',
                'p11':'Celsr1-Crsh-mCherry',
                }
            
            # Load CSV file 
            df = pd.read_csv(file_path,)
            
            # Add metadata columns

            df['filename'] = filename
            filename_parts = filename.split('.')[0].split('_')  # Remove extension and split by '_'
            # Extract metadata
            
            df['date'] = filename_parts[3]
            df['Celsr1_Plasmid'] = filename_parts[4][-6:-3] 
            df['Celsr1_Variant'] = df['Celsr1_Plasmid'].map(plasmid_mapping)
            df['Vangl2_Plasmid'] = filename_parts[4][-3:]
            df['Vangl2_Variant'] = df['Vangl2_Plasmid'].map(plasmid_mapping)
            df['Replicate'] = int(filename_parts[5])
            df['coverslip'] = int(filename_parts[6])
            df['image_number'] = filename_parts[7]
            df.rename(columns={"Max": "IntensityMax_Celsr1-Vangl2", "Mean": "IntensityMean_Celsr1-Vangl2"},inplace=True)
            df_all = pd.concat([df_all, df], ignore_index=True)
            
        except Exception as e:
            print(f"Error processing file {filename}: {str(e)}")
            continue

    
    return process_colocalization_data(df_all)

def plot_colocalization_results(df):
    """
    Plot colocalization results using Seaborn

    Parameters:
    df (pd.DataFrame): DataFrame containing colocalization results with columns 'Celsr1_Variant', 'FractionOfCelsr1_with_Vangl2'
    
    Returns:
    None
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    sns.set_theme(style="whitegrid")
    
    plt.figure(figsize=(10, 6))
    ax = sns.stripplot(x='Celsr1_Variant', y='FractionOfCelsr1_with_Vangl2',hue='Replicate', data=df)
    
    #ax.set_title('Fraction of Celsr1 Puncta with Vangl2 Colocalization by Variant')
    ax.set_xlabel('Celsr1 Variant')
    ax.set_ylabel('Fraction of Celsr1 Puncta with Vangl2')
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(r"C:/Users/wgiang/Documents/2026-04-07_Sarah-Vang-coloc/06_collated_results/COS7-Vangl2-coloc-stripplot.png")

# Example usage
def main():
    # Specify the directory containing your CSV files
    directory_path = r"C:/Users/wgiang/Documents/2026-04-07_Sarah-Vang-coloc/output/"  

    # Specify the directory where you want your output files to be saved
    output_path = r"C:/Users/wgiang/Documents/2026-04-07_Sarah-Vang-coloc/06_collated_results"

    # Process the files
    df = process_csv_files(directory_path)
    
    plot_colocalization_results(df)
if __name__ == "__main__":
    main()