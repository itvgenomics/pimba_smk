rule run_itsx:
    input:
        rawdata = rawdata
    output:
        itsx_fasta = os.path.join(output_dir, file_name_raw + '_itsx.fasta')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'itsx', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_itsx.txt")
    params:
        threads = config['num_threads']
    singularity:
        "docker://metashot/itsx:1.1.2-1"
    shell:
        """
        echo "Running the ITSx Container: "
        mkdir -p {output_dir}
        ITSx -i {input.rawdata} -o {output_dir}/{file_name_raw}_itsx -t f \
        --cpu {params.threads} --graphical F > {log} 2>&1
        
        cat {output_dir}/*_itsx.ITS1.fasta {output_dir}/*_itsx.ITS2.fasta > {output.itsx_fasta}
        """
rule run_vsearch_derep_fulllength:
    input:
        itsx_fasta = os.path.join(output_dir, file_name_raw + '_itsx.fasta')
    output:
        derep_fasta = os.path.join(output_dir, file_name_raw + '_derep.fasta')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'vsearch_derep', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_derep.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
    shell:
        """
        echo "Creating a VSEARCH Container: "
        vsearch --derep_fulllength {input.itsx_fasta} --output {output.derep_fasta} --sizeout > {log} 2>&1
        """