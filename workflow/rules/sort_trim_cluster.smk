rule run_vsearch_sort:
    input:
        derep = os.path.join(output_dir, file_name_raw + '_derep.fasta')
    output:
        sorted = os.path.join(output_dir, file_name_raw + '_sorted.fasta')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_sort", file_name_raw + ".log")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
    shell:
        """
        echo "Running the VSEARCH Container - --sortbysize: "
        vsearch --sortbysize {input.derep} --output {output.sorted} --minsize 2 > {log} 2>&1
        """

rule run_vsearch_trim:
    input:
        sorted = os.path.join(output_dir, file_name_raw + '_sorted.fasta')
    output:
        trimmed = os.path.join(output_dir, file_name_raw + '_trimmed.fasta')
    params:
        otu_length = config['otu_length']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_trim", file_name_raw + ".log")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
    shell:
        """
        if [ {params.otu_length} != "0" ]; then
            echo "Running the VSEARCH Container - --fastx_filter: "
            vsearch --fastx_filter {input.sorted} --fastq_trunclen {params.otu_length} \
            --fastaout {output.trimmed} > {log} 2>&1
        else
            cp {input.sorted} {output.trimmed}
        fi
        """

if strategy == "otu":
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
        singularity:
            "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
        shell:
            """
            echo "Running the VSEARCH Container - --cluster_size: "
            vsearch --cluster_size {input.trimmed} --consout {output.otus} \
            --id {params.similarity} --threads {params.threads} > {log} 2>&1
            """

if strategy == "asv":
    rule run_vsearch_asv_fastx:
        input:
            trimmed = os.path.join(output_dir, file_name_raw + '_trimmed.fasta')
        output:
            trimmed_noN = os.path.join(output_dir, file_name_raw + '_trimmed_noN.fasta')
        params:
            threads = config['num_threads']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "vsearch_asv_fastx", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
        shell:
            """
            echo "Running the VSEARCH Container - --fastx_filter: "
            vsearch --fastx_filter {input.trimmed} --fastq_maxns 0 --fastaout {output.trimmed_noN} > {log} 2>&1
            """
    rule swarm_derep:
        input:
            trimmed_noN = os.path.join(output_dir, file_name_raw + '_trimmed_noN.fasta')
        output:
            derep = os.path.join(output_dir, file_name_raw + '_trimmed_noN_derep.fasta')
        params:
            threads = config['num_threads']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "swarm_derep", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_swarm:v3.1.0"
        shell:
            """
            echo "Running the SWARM Container: "
            swarm -d 0 -z -t {params.threads} -w {output.derep} -o /dev/null {input.trimmed_noN} > {log} 2>&1
            """
    rule find_asvs_with_swarm:
        input:
            derep = os.path.join(output_dir, file_name_raw + '_trimmed_noN_derep.fasta')
        output:
            asvs = os.path.join(output_dir, file_name_raw + '_asvs1.fasta')
        params:
            threads = config['num_threads']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "find_asvs_with_swarm", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_swarm:v3.1.0"
        shell:
            """
            swarm -f -t {params.threads} -w {output.asvs} -u {output_dir}/uclust_format.file -z \
            -d 1 {input.derep} > {log} 2>&1
            """