rule run_lulu:
    input:
        table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt'),
        match_list = os.path.join(output_dir, 'lulu_output', 'match_list.txt')
    output:
        lulu_table = os.path.join(output_dir, 'lulu_output', file_name_raw + '_otu_table_lulu.txt') if strategy == "otu" else os.path.join(output_dir, 'lulu_output', file_name_raw + '_asv_table_lulu.txt'),
    params:
        old_table = os.path.join(output_dir, file_name_raw + '_otu_table_old.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_old.txt')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "lulu", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_lulu.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_r:latest"
    shell:
        """
        echo "Running the R Container - lulu: "
        Rscript /data/file.R {input.table} {input.match_list} {output.lulu_table} > {log} 2>&1

        sed -i -e 's/"X//g' {output.lulu_table}
        sed -i -e 's/"//g' {output.lulu_table}
        sed -i '1s/^/OTUId\t/' {output.lulu_table}

        cp {input.table} {params.old_table}
        cp {output.lulu_table} {input.table}
        """