if marker_gene == "ITS-FUNGI-NCBI" or marker_gene == "ITS-PLANTS-NCBI":
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
        singularity:
            "docker://itvdsbioinfo/pimba_vsearch:v2.29.1"
        shell:
            """
            echo "Running the VSEARCH Container - --usearch_global: "
            vsearch --usearch_global {input.rawdata} --db {input.renamed} --strand both \
            --id {params.similarity} --uc {output.mapfile} --otutabout {output.table_txt} --threads {params.threads} > {log} 2>&1
            """
else:
    rule run_vsearch_map:
        input:
            renamed = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta'),
            rawdata = rawdata
        output:
            mapfile = os.path.join(output_dir, file_name_raw + '_map.uc'),
            table_txt = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
        params:
            similarity = config['otu_similarity'],
            threads = config['num_threads']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "vsearch_map", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_vsearch:v2.29.1"
        shell:
            """
            echo "Running the VSEARCH Container - --usearch_global: "
            vsearch --usearch_global {input.rawdata} --db {input.renamed} --strand both \
            --id {params.similarity} --uc {output.mapfile} --otutabout {output.table_txt} --threads {params.threads} > {log} 2>&1
            """

if lulu == "yes":
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
    if marker_gene == "ITS-FUNGI-NCBI" or marker_gene == "ITS-PLANTS-NCBI":
        rule run_lulu:
            input:
                table = os.path.join(output_dir, file_name_raw + '_pre_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_pre_asv_table.txt'),
                match_list = os.path.join(output_dir, 'lulu_output', 'match_list.txt')
            output:
                lulu_table = os.path.join(output_dir, 'lulu_output', file_name_raw + '_otu_table_lulu.txt') if strategy == "otu" else os.path.join(output_dir, 'lulu_output', file_name_raw + '_asv_table_lulu.txt'),
            params:
                old_table = os.path.join(output_dir, file_name_raw + '_otu_table_old.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_old.txt')
            log:
                os.path.join(current_path, "results", "01-run", "logs", "lulu", file_name_raw + ".log")
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
    else:
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