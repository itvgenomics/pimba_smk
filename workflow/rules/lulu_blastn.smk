rule run_makeblastdb_and_blastn:
    input:
        fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
    output:
        blast_db = os.path.join(output_dir, 'lulu_output', file_name_raw + '_otus.fasta.nsq'),
        match_list = os.path.join(output_dir, 'lulu_output', 'match_list.txt')
    params:
        threads = config['num_threads']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "makeblastdb_blastn", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_makeblastdb_blastn.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_blast:latest"
    shell:
        """
        mkdir -p {output_dir}/lulu_output
        echo "Running the BLAST Container - makeblastdb: "
        makeblastdb -in {input.fasta_file} -out {output_dir}/lulu_output/{file_name_raw}_otus.fasta \
        -parse_seqids -dbtype nucl > {log} 2>&1

        echo "Running the BLAST Container - blastn: "
        blastn -db {output_dir}/lulu_output/{file_name_raw}_otus.fasta -outfmt "6 qseqid sseqid pident" \
        -out {output_dir}/lulu_output/match_list.txt -num_threads {params.threads} \
        -qcov_hsp_perc 80 -perc_identity 84 -query {input.fasta_file} >> {log} 2>&1
        """