rule run_plot_graphs:
    input:
        otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt'),
        prepared_tax_assignments = os.path.join(output_dir, 'output', 'plot_' + file_name_raw + '_tax_assignments.txt'),
        metadata = config['metadata']
    output:
        os.path.join(current_path, 'results', '02-plot', 'plots', 'alpha_diversity_dotplot.svg')
    log:
        os.path.join(current_path, 'results', '02-plot', 'logs', 'plot_graphs', file_name_raw + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_python_plot:latest"
    params:
        group_by = config['group_by']
    shell:
        """
        mkdir -p {current_path}/results/02-plot/
        cd {current_path}/results/02-plot/
        echo "Running the plot python container: "
        eval "$(conda shell.bash hook)"
        conda activate plot_env
        python /app/pimba_plot_graphs.py -o {input.otu_table} -t {input.prepared_tax_assignments} -m {input.metadata} -g {params.group_by} \
        > {log} 2>&1
        """