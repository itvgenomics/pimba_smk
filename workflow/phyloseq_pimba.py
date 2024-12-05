# Import necessary libraries
import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from skbio import TreeNode
from skbio.stats.ordination import pcoa
import seaborn as sns
import skbio.diversity
from scipy.special import gammaln
from skbio import DistanceMatrix
from sklearn.manifold import MDS

def load_data(otu_file, tax_file, metadata_file):
    # Load OTU table
    otu_table_df = pd.read_csv(otu_file, sep='\t', index_col=0).sort_index()
    # Load and fix taxonomy table
    def is_numeric_string(s):
    # Check if the string represents a numeric value
        try:
            float(s)
            return True
        except ValueError:
            return False
    tax_table_df = pd.read_csv(tax_file, sep='\t', index_col=0).sort_index()
    tax_table_df = tax_table_df.drop(columns=['similarity', 'aux'])
    tax_table_df.replace(' ', 'Unassigned', inplace=True)
    tax_table_df = tax_table_df.applymap(lambda x: 'Unassigned' if isinstance(x, str) and is_numeric_string(x) else x)
    # Load metadata
    meta_table_df = pd.read_csv(metadata_file, sep=',', index_col=0).sort_index()
    return otu_table_df, tax_table_df, meta_table_df

def plot_taxonomic_composition(otu_table_df, tax_table_df, meta_table_df, rank, groupby, output_dir):
    # Copy and transform the OTU table
    ps_taxa = otu_table_df.copy()
    ps_taxa[rank] = tax_table_df[rank]
    ps_taxa = ps_taxa.groupby(rank).sum()
    ps_taxa = ps_taxa / ps_taxa.sum()  # Transform to relative abundance
    ps_taxa = ps_taxa.reset_index()
    ps_taxa = pd.melt(ps_taxa, id_vars=rank, var_name='Sample', value_name='Abundance')
    ps_taxa = ps_taxa[ps_taxa['Abundance'] > 0.01]  # Filter out low abundance taxa
    ps_taxa = ps_taxa.sort_values(by=rank)  # Sort data frame alphabetically by the taxonomic rank
    
    # Merge metadata into the transformed OTU table
    ps_taxa = ps_taxa.merge(meta_table_df, left_on='Sample', right_index=True, how='left')
    
    # Define the number of unique taxa and generate the color palette
    n = len(ps_taxa[rank].unique()) + 1
    palette = sns.color_palette("tab20", n)
    
    # Set plot size and font scale
    sns.set_context("talk", font_scale=1.2)
    
    # Plotting with seaborn FacetGrid for faceting
    if groupby:
        # Pivot the dataframe to get the format needed for stacked bars
        pivot_df = ps_taxa.pivot_table(index=['Sample', groupby], columns=rank, values='Abundance').fillna(0).reset_index()
        
        # Calculate low abundance taxa per sample and add to pivot_df
        pivot_df['Low abundance taxa'] = 1 - pivot_df.drop(columns=['Sample', groupby]).sum(axis=1)
        
        # Save pivot_df to a CSV file
        pivot_df.to_csv(os.path.join(output_dir, f'{rank}_pivot_table.tsv'), sep='\t')
        
        # Initialize the FacetGrid
        g = sns.FacetGrid(pivot_df, col=groupby, col_wrap=3, sharey=False, sharex=False, height=8, aspect=0.5)
        
        # Plot the stacked bar plot within each facet
        for ax, (name, df) in zip(g.axes.flat, pivot_df.groupby(groupby)):
            df.set_index('Sample').drop(columns=[groupby]).plot(kind='bar', stacked=True, color=palette, ax=ax, legend=False)
            ax.set_title(name, fontsize=18)
            ax.set_xlabel("Sample", fontsize=14)
            ax.set_ylabel(f"Relative Abundance ({rank.capitalize()})", fontsize=14)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=12)
            ax.set_yticks(ax.get_yticks())  # Ensure that ticks are set before setting labels
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
            ax.set_ylim(0, 1.0)  # Set the y-axis limits for each subplot

        # Manually create the legend outside the loop
        handles, labels = ax.get_legend_handles_labels()  # Get handles and labels from the last ax
        g.fig.legend(handles, labels, title=rank.capitalize(), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='16', title_fontsize='18')
        plt.subplots_adjust(top=0.9, hspace=0.3)  # Adjust the top to make room for the title
        g.fig.suptitle(f'{rank.capitalize()} Stacked Barplot', fontsize=18)
        g.fig.savefig(os.path.join(output_dir, f'{rank}_barplots.svg'), bbox_inches='tight')
    else:
        # Pivot the dataframe to get the format needed for stacked bars
        pivot_df = ps_taxa.pivot(index='Sample', columns=rank, values='Abundance').fillna(0)
        
        # Calculate low abundance taxa per sample and add to pivot_df
        pivot_df['Low abundance taxa'] = 1 - pivot_df.sum(axis=1)
        
        # Save pivot_df to a CSV file
        pivot_df.to_csv(os.path.join(output_dir, f'{rank}_pivot_table.tsv'), sep='\t')
        
        # Plotting the stacked bar plot
        pivot_df.plot(kind='bar', stacked=True, color=palette, figsize=(20, 10))
        plt.ylabel(f'Relative Abundance ({rank.capitalize()})', fontsize=18)
        plt.xlabel('Sample', fontsize=18)
        plt.xticks(rotation=30, ha='right', fontsize=16)
        plt.yticks(fontsize=16)
        plt.legend(title=rank.capitalize(), bbox_to_anchor=(1.05, 1), loc='upper left', title_fontsize='18', fontsize='14')
        plt.title(f'{rank.capitalize()} Stacked Barplot', fontsize=22)
        plt.savefig(os.path.join(output_dir, f'{rank}_barplots.svg'), bbox_inches='tight')

def plot_alpha_diversity(otu_table_df, output_dir):
    def calculate_alpha_diversity(otu_table):
        # Calculate observed species
        observed = otu_table.apply(lambda x: x[x > 0].count(), axis=0)
        observed_se = otu_table.apply(lambda x: x.std() / np.sqrt(len(x)), axis=0)

        # Calculate Chao1
        chao1 = otu_table.apply(lambda x: skbio.diversity.alpha.chao1(x), axis=0)
        chao1_se = np.sqrt((chao1 * (1 / len(otu_table))))  # Adjust standard error calculation

        # Calculate Shannon index
        shannon = otu_table.apply(lambda x: skbio.diversity.alpha.shannon(x), axis=0)
        shannon_se = np.sqrt((shannon * (1 / len(otu_table))))  # Adjust standard error calculation

        return (observed, observed_se), (chao1, chao1_se), (shannon, shannon_se)

    # Calculate the alpha diversity measures and their standard errors
    (observed, observed_se), (chao1, chao1_se), (shannon, shannon_se) = calculate_alpha_diversity(otu_table_df)

    # Create a DataFrame for plotting
    alpha_diversity_df = pd.DataFrame({
        'Sample': otu_table_df.columns,  # Using the sample names from the columns
        'Observed': observed,
        'Observed_se': observed_se,
        'Chao1': chao1,
        'Chao1_se': chao1_se,
        'Shannon': shannon,
        'Shannon_se': shannon_se,
    })

    # Set the aesthetics for the plots
    sns.set(style='whitegrid')

    # Create a figure with three subplots
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(24, 6))

    # Create a boxplot for each alpha diversity metric with standard error bars
    metrics = ['Observed', 'Chao1', 'Shannon']
    for i, metric in enumerate(metrics):
        ax = axes[i]
        sns.boxplot(data=alpha_diversity_df, x='Sample', y=metric, color='lightgray', boxprops=dict(alpha=.5), ax=ax)
        
        # Overlay the standard error bars
        ax.errorbar(x=np.arange(len(alpha_diversity_df['Sample'].unique())), 
                    y=alpha_diversity_df.groupby('Sample')[metric].mean(),
                    yerr=alpha_diversity_df.groupby('Sample')[f'{metric}_se'].mean(),
                    fmt='o', color='red', capsize=5, label='Standard Error', markersize=5)
        
        ax.set_title(f'Boxplot of {metric} Diversity Metric with Standard Error')
        ax.set_ylabel(metric)
        ax.grid(True)
        ax.legend()
        # Set the ticks and rotate x-axis labels
        ax.set_xticks(np.arange(len(alpha_diversity_df['Sample'].unique())))
        ax.set_xticklabels(alpha_diversity_df['Sample'].unique(), rotation=90)

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the figure to a file
    plt.savefig(os.path.join(output_dir, 'alpha_diversity_dotplot.svg'))

def plot_rarefaction_curves(otu_table_df, output_dir):
    # Transpose the OTU table
    otu_dataFC = otu_table_df.T

    # Calculate the number of species (S)
    S = (otu_dataFC > 0).sum(axis=1)

    # Calculate the minimum row sum (raremax)
    raremax = otu_dataFC.sum(axis=1).min()

    # Define the rarefy function
    def rarefy(x, sample):
        x = x[x > 0]
        n = x.sum()
        if sample > n:
            print(f"Warning: specified 'sample' size is larger than the observed number of individuals in one or more sites")
            sample = n
        p = sample / n
        log_comb = lambda n, k: gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        rare = sum(1 - np.exp(log_comb(n - x, sample) - log_comb(n, sample)))
        return rare

    # Apply rarefy function to each row of the OTU table
    Srare = otu_dataFC.apply(lambda row: rarefy(row, raremax), axis=1)

    # Plotting parameters
    lty = ['-', '--', '-.', ':']
    palette = plt.cm.tab20.colors

    # Generate the parameter grid
    pars = [{'col': palette[i % len(palette)], 'lty': lty[i % len(lty)]} for i in range(len(otu_dataFC))]

    # Create the plot
    plt.figure(figsize=(15, 8))
    for i, (idx, row) in enumerate(otu_dataFC.iterrows()):
        col = pars[i]['col']
        linestyle = pars[i]['lty']
        rare_curve = [rarefy(row, n) for n in range(1, int(raremax) + 1, 50)]
        plt.plot(range(1, int(raremax) + 1, 50), rare_curve, label=idx, color=col, linestyle=linestyle)

    # Adding labels and legend
    plt.xlabel('Sample Size')
    plt.ylabel('Species')
    plt.title('Rarefaction Curves')
    plt.legend(loc='best')

    # Save the plot as an SVG file
    plt.savefig(os.path.join(output_dir, 'rarefaction_curve2.svg'))

def plot_cluster_dendrogram(otu_table_df, meta_table_df, output_dir, groupby):
    # Transpose OTU table to have samples as rows
    otu_data = otu_table_df.T

    # Function to calculate Bray-Curtis distance
    def bray_curtis_distance(otu_table):
        """Calculate the Bray-Curtis distance matrix."""
        otu_df = pd.DataFrame(otu_table)
        # Calculate the total for each row and create a normalized DataFrame
        row_sums = otu_df.sum(axis=1)
        normalized_df = otu_df.div(row_sums, axis=0)  # Normalize by row sums
        # Calculate the Bray-Curtis distance matrix
        return 1 - normalized_df.dot(normalized_df.T)

    # Compute the distance matrix
    distance_matrix = bray_curtis_distance(otu_data)

    # Hierarchical clustering using average linkage
    linked = linkage(distance_matrix, method='average')

    # Reorder samples and append group information if groupby is specified
    if groupby and groupby in meta_table_df.columns:
        order = dendrogram(linked, no_plot=True)['leaves']
        reordered_samples = [otu_table_df.columns[i] for i in order]
        reordered_groups = meta_table_df.loc[reordered_samples, groupby]
        labels = [f"{group}-{sample}" for group, sample in zip(reordered_groups, reordered_samples)]
    else:
        labels = otu_table_df.columns

    # Plotting the dendrogram with new labels
    plt.figure(figsize=(15, 8))
    dendrogram(linked, labels=labels, orientation='top')
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Samples')
    plt.xticks(rotation=90)
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cluster_dendrogram.svg"))

def perform_pvclust(otu_table_df, meta_table_df, output_dir, groupby, method_dist='correlation', method_hclust='ward', nboot=1000):
    # Function to perform hierarchical clustering with bootstrap resampling
    def pvclust(data, method_dist, method_hclust, nboot):
        def hierarchical_clustering(data, method_dist, method_hclust):
            if method_dist == 'correlation':
                dist = pdist(data.T, metric='correlation')
            else:
                dist = pdist(data.T, metric=method_dist)
            hc = linkage(dist, method=method_hclust)
            return hc

        def bootstrap(data, method_dist, method_hclust, nboot):
            boot_results = []
            n = data.shape[1]
            for _ in range(nboot):
                indices = np.random.choice(n, n, replace=True)
                boot_data = data[:, indices]
                hc = hierarchical_clustering(boot_data, method_dist, method_hclust)
                boot_results.append(hc)
            return boot_results

        hc = hierarchical_clustering(data, method_dist, method_hclust)
        boot_results = bootstrap(data, method_dist, method_hclust, nboot)
        
        return hc, boot_results

    # Function to calculate AU (Approximately Unbiased) p-values
    def calculate_au_pvalues(hc, boot_results, data):
        n = data.shape[1]
        au_pvalues = np.zeros(len(hc))

        for i in range(len(hc)):
            original_cluster = fcluster(hc, t=hc[i, 2], criterion='distance')
            count_au = 0
            for boot_hc in boot_results:
                boot_cluster = fcluster(boot_hc, t=hc[i, 2], criterion='distance')
                if len(np.unique(original_cluster)) == len(np.unique(boot_cluster)):
                    count_au += 1
            au_pvalues[i] = count_au / len(boot_results) * 100

        return au_pvalues

    # Function to calculate bootstrap values
    def calculate_bootstrap_values(hc, boot_results):
        n = len(hc) + 1  # number of original observations
        bootstrap_values = np.zeros(n - 1)

        for i in range(n - 1):
            original_cluster = fcluster(hc, t=hc[i, 2], criterion='distance')
            count = 0
            for boot_hc in boot_results:
                boot_cluster = fcluster(boot_hc, t=hc[i, 2], criterion='distance')
                if np.array_equal(np.sort(original_cluster), np.sort(boot_cluster)):
                    count += 1
            bootstrap_values[i] = count / len(boot_results) * 100

        return bootstrap_values

    # Function to plot the dendrogram with both AU p-values and bootstrap values
    def plot_pvclust(hc, boot_results, data, labels, output_dir, alpha=0.95):
        plt.figure(figsize=(15, 8))
        dend = dendrogram(hc, labels=labels, leaf_rotation=90)
        
        # Calculate AU p-values and bootstrap values
        au_pvalues = calculate_au_pvalues(hc, boot_results, data)
        bootstrap_values = calculate_bootstrap_values(hc, boot_results)
        
        # Add AU p-values and bootstrap values to the dendrogram
        for i, (x_coords, y_coords) in enumerate(zip(dend['icoord'], dend['dcoord'])):
            x = 0.5 * sum(x_coords[1:3])
            y = y_coords[1]
            
            # Display AU p-value
            plt.text(x, y, f'AU: {au_pvalues[i]:.2f}', va='bottom', ha='center', fontsize=10, color='blue')

            # Display Bootstrap value
            plt.text(x, y - 0.05, f'Bootstrap: {bootstrap_values[i]:.2f}', va='bottom', ha='center', fontsize=10, color='red')

        plt.title("Cluster analysis with AU and Bootstrap values")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "cluster_bootstrap1000.svg"))

    # Convert OTU table to matrix and ensure all values are finite
    taxon_matrix = otu_table_df.values

    # Perform clustering
    hc, boot_results = pvclust(taxon_matrix, method_dist, method_hclust, nboot)

    # Optionally reorder labels based on a grouping variable
    if groupby:
        order = [i[0] for i in sorted(enumerate(hc[:, 2]), key=lambda x: x[1])]
        hell_tip_labels = meta_table_df.loc[otu_table_df.columns, groupby].values
        labels = [f"{hell_tip_labels[i]}-{meta_table_df.loc[otu_table_df.columns[i], 'SampleName']}" for i in range(len(hell_tip_labels))]
    else:
        labels = otu_table_df.columns

    # Plot the result
    plot_pvclust(hc, boot_results, taxon_matrix, labels, output_dir)

def perform_pcoa_unweighted_unifrac(otu_table_df, meta_table_df, output_dir, groupby):
    # Function to create a random tree using a simple balanced tree approach
    def create_random_tree(tips):
        tree_str = "(" + ",".join(tips) + ");"  # Simple star-like tree
        tree = TreeNode.read([tree_str])
        return tree

    # Create a random tree using OTU names
    random_tree = create_random_tree(otu_table_df.index.tolist())

    # Function to calculate unweighted UniFrac distance manually
    def calculate_unweighted_unifrac(otu_table, tree):
        # Normalize the OTU table
        otu_table_normalized = otu_table.div(otu_table.sum(axis=0), axis=1)

        # Prepare the UniFrac distance matrix
        n_samples = otu_table_normalized.shape[1]
        unifrac_dist_matrix = np.zeros((n_samples, n_samples))

        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                # Find the shared and unique OTUs
                otu_i = otu_table_normalized.iloc[:, i]
                otu_j = otu_table_normalized.iloc[:, j]
                shared_otus = (otu_i > 0) & (otu_j > 0)

                # Calculate the UniFrac distance
                unique_i = otu_i[~shared_otus].sum()
                unique_j = otu_j[~shared_otus].sum()
                unifrac_dist = (unique_i + unique_j) / (otu_i.sum() + otu_j.sum())

                # Fill in the distance matrix
                unifrac_dist_matrix[i, j] = unifrac_dist
                unifrac_dist_matrix[j, i] = unifrac_dist  # Symmetric matrix

        return DistanceMatrix(unifrac_dist_matrix, otu_table.columns.tolist())

    # Calculate the distance matrix
    unifrac_dist = calculate_unweighted_unifrac(otu_table_df, random_tree)

    # Perform PCoA
    ordination_results = pcoa(unifrac_dist)

    # Create a DataFrame for the results
    pcoa_df = pd.DataFrame(ordination_results.samples, index=otu_table_df.columns.tolist())
    pcoa_df['SampleName'] = otu_table_df.columns.tolist()  # Add sample names for labeling

    # Add group information if groupby is specified
    if groupby and groupby in meta_table_df.columns:
        pcoa_df[groupby] = meta_table_df.loc[pcoa_df['SampleName'], groupby]

    # Plotting the PCoA results
    plt.figure(figsize=(15, 8))
    sns.scatterplot(data=pcoa_df, x='PC1', y='PC2', hue=groupby if groupby and groupby in pcoa_df.columns else None, palette='Set2', s=300)

    # Annotate each point with the sample name
    for i, sample_name in enumerate(pcoa_df['SampleName']):
        plt.text(pcoa_df['PC1'].iloc[i], pcoa_df['PC2'].iloc[i], sample_name,
                 horizontalalignment='left', size='medium', color='black', weight='semibold')

    plt.title('PCoA with Unweighted UniFrac Distance')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.axhline(0, color='grey', lw=0.5, ls='--')
    plt.axvline(0, color='grey', lw=0.5, ls='--')
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "PCoA_unweighted_unifrac.svg"))

def perform_nmds_bray_curtis(otu_table_df, meta_table_df, output_dir, groupby):
    # Function to calculate Bray-Curtis distance
    def calculate_bray_curtis(otu_table):
        # Normalize the OTU table
        otu_table_normalized = otu_table.div(otu_table.sum(axis=0), axis=1)

        # Create a distance matrix
        n_samples = otu_table_normalized.shape[1]
        bray_curtis_matrix = np.zeros((n_samples, n_samples))

        for i in range(n_samples):
            for j in range(n_samples):
                # Calculate Bray-Curtis distance
                if i != j:
                    bray_curtis_dist = (otu_table_normalized.iloc[:, i] - otu_table_normalized.iloc[:, j]).abs().sum() / (otu_table_normalized.iloc[:, i] + otu_table_normalized.iloc[:, j]).sum()
                    bray_curtis_matrix[i, j] = bray_curtis_dist
                else:
                    bray_curtis_matrix[i, j] = 0.0  # Distance to itself is zero

        return DistanceMatrix(bray_curtis_matrix, otu_table.columns.tolist())

    # Calculate the distance matrix
    bray_curtis_dist = calculate_bray_curtis(otu_table_df)

    # Perform NMDS using scikit-learn
    nmds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    nmds_results = nmds.fit_transform(bray_curtis_dist.data)

    # Create a DataFrame for the results
    nmds_df = pd.DataFrame(nmds_results, index=otu_table_df.columns.tolist(), columns=['NMDS1', 'NMDS2'])
    nmds_df['SampleName'] = otu_table_df.columns.tolist()  # Add sample names for labeling

    # Add group information if groupby is specified
    if groupby and groupby in meta_table_df.columns:
        nmds_df[groupby] = meta_table_df.loc[nmds_df['SampleName'], groupby]

    # Plotting the NMDS results
    plt.figure(figsize=(15, 8))
    sns.scatterplot(data=nmds_df, x='NMDS1', y='NMDS2', hue=groupby if groupby and groupby in nmds_df.columns else None, palette='Set2', s=300)

    # Annotate each point with the sample name
    for i, sample_name in enumerate(nmds_df['SampleName']):
        plt.text(nmds_df['NMDS1'].iloc[i], nmds_df['NMDS2'].iloc[i], sample_name,
                 horizontalalignment='left', size='medium', color='black', weight='semibold')

    # Title and labels
    plt.title('NMDS with Bray-Curtis Distance')
    plt.xlabel('NMDS1')
    plt.ylabel('NMDS2')
    plt.axhline(0, color='grey', lw=0.5, ls='--')
    plt.axvline(0, color='grey', lw=0.5, ls='--')
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "NMDS_bray.svg"))

def main():
    parser = argparse.ArgumentParser(description='Process OTU table, taxonomic assignments, and metadata.')
    parser.add_argument('-o', '--otu', required=True, help='OTU table input file')
    parser.add_argument('-t', '--tax', required=True, help='Taxonomic assignments input file')
    parser.add_argument('-m', '--metadata', required=True, help='Metadata input file')
    parser.add_argument('-g', '--group_by', required=True, help='Group by flag (e.g., Species, SampleName). If you do not want to group samples, set this variable to False')
    args = parser.parse_args()

    otu_table_df, tax_table_df, meta_table_df = load_data(args.otu, args.tax, args.metadata)

    if args.group_by == 'False':
        args.group_by = False
    else:
        if args.group_by not in meta_table_df.columns:
            raise ValueError(f"Group by column '{args.group_by}' not found in metadata")

    output_dir = 'plots'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plot_taxonomic_composition(otu_table_df, tax_table_df, meta_table_df, 'phylum', args.group_by, output_dir)
    plot_taxonomic_composition(otu_table_df, tax_table_df, meta_table_df, 'class', args.group_by, output_dir)
    plot_taxonomic_composition(otu_table_df, tax_table_df, meta_table_df, 'order', args.group_by, output_dir)
    plot_taxonomic_composition(otu_table_df, tax_table_df, meta_table_df, 'family', args.group_by, output_dir)
    plot_taxonomic_composition(otu_table_df, tax_table_df, meta_table_df, 'genus', args.group_by, output_dir)

    plot_alpha_diversity(otu_table_df, output_dir)

    plot_rarefaction_curves(otu_table_df, output_dir)

    plot_cluster_dendrogram(otu_table_df, meta_table_df, output_dir, args.group_by)

    perform_pvclust(otu_table_df, meta_table_df, output_dir, args.group_by)

    perform_pcoa_unweighted_unifrac(otu_table_df, meta_table_df, output_dir, args.group_by)

    perform_nmds_bray_curtis(otu_table_df, meta_table_df, output_dir, args.group_by)


if __name__ == '__main__':
    main()
