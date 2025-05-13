rule run_qiime_assign_taxonomy:
    input:
        fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
    output:
        taxonomy = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    params:
        db_dir = config['ITS-FUNGI-UNITE-DB'],
        similarity_int = similarity_int,
        similarity_assign = config['assign_similarity']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "qiime_assign_taxonomy", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_qiime_assign_taxonomy.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiime:latest"
    shell:
        """
        echo "Running the Qiime Container - assign_taxonomy.py: "
        assign_taxonomy.py -i {input.fasta_file} -o {output_dir}/output/ \
        -t {params.db_dir}/*{params.similarity_int}*.txt \
        -r {params.db_dir}/*{params.similarity_int}*.fasta \
        --similarity={params.similarity_assign} > {log} 2>&1
        """
rule run_create_abundance_file:
    input:
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
        otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
    output:
        abundance_file_dir = directory(os.path.join(output_dir, 'diversity_by_sample'))
    log:
        os.path.join(current_path, "results", "01-run", "logs", "create_abundance_file", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_create_abundance_file.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        mkdir -p {output.abundance_file_dir}
        cd {output.abundance_file_dir}
        echo "Running the QiimePipe Container - createAbundanceFile.py: "
        python /qiimepipe/createAbundanceFile.py {input.tax_assignments} {input.otu_table} > {log} 2>&1
        """