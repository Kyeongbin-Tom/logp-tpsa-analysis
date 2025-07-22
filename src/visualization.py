import matplotlib.pyplot as plt
import seaborn as sns

def plot_tpsa_vs_logp(df, tpsa_col='TPSA', logp_col='Lipophilicity', save_path=None):
    """
    Plot scatter plot of TPSA vs LogP.
    
    Parameters:
    - df: DataFrame with TPSA and LogP values
    - tpsa_col: Column name for TPSA (default: 'TPSA')
    - logp_col: Column name for LogP (default: 'Lipophilicity')
    - save_path: If specified, saves the plot to the given path
    """
    plt.figure(figsize=(6, 4))
    sns.scatterplot(x=df[tpsa_col], y=df[logp_col], alpha=0.6)
    plt.title("TPSA vs LogP")
    plt.xlabel("TPSA")
    plt.ylabel("LogP")
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_distribution(df, col, title=None, bins=30, save_path=None):
    """
    Plot histogram of a single numerical column.
    
    Parameters:
    - df: DataFrame
    - col: Column name to plot
    - title: Plot title
    - bins: Number of histogram bins
    - save_path: If specified, saves the plot to the given path
    """
    plt.figure(figsize=(6, 4))
    sns.histplot(df[col], bins=bins, kde=True, color='skyblue')
    plt.xlabel(col)
    plt.ylabel("Frequency")
    plt.title(title if title else f"Distribution of {col}")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()
