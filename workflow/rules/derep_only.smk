rule run_vsearch_derep_fulllength:
    input:
        rawdata = rawdata
    output:
        derep = os.path.join(output_dir, file_name_raw + '_derep.fasta')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_derep", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_derep.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
    shell:
        """
        echo "Running the VSEARCH Container - --derep_fulllength:"
        mkdir -p {output_dir}
        vsearch --derep_fulllength {input.rawdata} --output {output.derep} --sizeout > {log} 2>&1
        """