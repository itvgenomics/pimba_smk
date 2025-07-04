rule run_vsearch_map:
    input:
        renamed = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta'),
        rawdata = rawdata
    output:
        mapfile = os.path.join(output_dir, file_name_raw + '_map.uc'),
        table_txt = os.path.join(output_dir, file_name_raw + '_pre_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_pre_asv_table.txt')
    params:
        similarity = config['otu_similarity'],
        threads = config['num_threads']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_map", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_map.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.29.1"
    shell:
        """
        echo "Running the VSEARCH Container - --usearch_global: "
        vsearch --usearch_global {input.rawdata} --db {input.renamed} --strand both \
            --id {params.similarity} --uc {output.mapfile} \
            --otutabout {output.table_txt}.raw --threads {params.threads} > {log} 2>&1
        awk 'NR==1 {{print; next}} {{sum=0; for (i=2; i<=NF; i++) sum+=$i; if (sum>0) print}}' \
            {output.table_txt}.raw > {output.table_txt}
        rm -f {output.table_txt}.raw
        """