rule run_vsearch_sort:
    input:
        derep = os.path.join(output_dir, file_name_raw + '_derep.fasta')
    output:
        sorted = os.path.join(output_dir, file_name_raw + '_sorted.fasta')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_sort", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_sort.txt")
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
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_trim.txt")
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