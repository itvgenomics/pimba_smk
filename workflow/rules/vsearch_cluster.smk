rule run_vsearch_otu_clustering:
    input:
        trimmed = os.path.join(output_dir, file_name_raw + '_trimmed.fasta')
    output:
        otus = os.path.join(output_dir, file_name_raw + '_otus1.fasta')
    params:
        similarity = config['otu_similarity'],
        threads = config['num_threads']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_otu_clustering", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_otu_clustering.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
    shell:
        """
        echo "Running the VSEARCH Container - --cluster_size: "
        vsearch --cluster_size {input.trimmed} --consout {output.otus} \
        --id {params.similarity} --threads {params.threads} > {log} 2>&1
        """