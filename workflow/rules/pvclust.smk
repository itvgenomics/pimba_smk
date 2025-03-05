rule run_pvclust:
    input:
        otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt'),
        prepared_tax_assignments = os.path.join(output_dir, 'output', 'plot_' + file_name_raw + '_tax_assignments.txt'),
        metadata = config['metadata']
    output:
        os.path.join(current_path, 'results', '02-plot', 'plots', 'cluster_bootstrap1000.svg')
    log:
        os.path.join(current_path, 'results', '02-plot', 'logs', 'pvclust', file_name_raw + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_phyloseq:latest"
    params:
        group_by = config['group_by']
    shell:
        """
        mkdir -p {current_path}/results/02-plot/
        cd {current_path}/results/02-plot/plots
        echo "Running the pvclust container: "
        Rscript /data/pvclust.R {input.otu_table} {input.prepared_tax_assignments} {input.metadata} {params.group_by} \
        > {log} 2>&1
        """