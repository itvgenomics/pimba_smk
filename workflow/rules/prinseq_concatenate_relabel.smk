rule run_prinseq:
    input:
        assembled = os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.assembled.fastq")
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}.assembled.fasta"),
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "prinseq", "{sample}.log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "{sample}_prepare_prinseq.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_prinseq:v0.20.4"
    shell:
        """
        export LANG=C
        export LC_ALL=C
        prinseq-lite.pl -fastq {input.assembled} -out_format 1 -seq_id Seq -out_good {current_path}/results/00-prepare/assemblies/pear/{wildcards.sample}.assembled \
        > {log} 2>&1
        """

rule run_relabel:
    input:
        assembled_fasta = os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}.assembled.fasta")
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}_relabel_notSingleton.fasta")
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "qiimepipe", "{sample}.log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "{sample}_prepare_relabel.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        cd {current_path}/results/00-prepare/assemblies/pear/
        python3.6 /qiimepipe/relabelReads-v2.py {wildcards.sample}.assembled.fasta . \
        > {log} 2>&1
        mv {wildcards.sample}_relabel.fasta {wildcards.sample}_relabel_notSingleton.fasta
        """

rule concatenate_and_rename:
    input:
        expand(os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}_relabel_notSingleton.fasta"), sample=SAMPLES)
    output:
        os.path.join(current_path, "results", "00-prepare", config["outputprepare"] + ".fasta")
    params:
        outputname = config["outputprepare"]
    shell:
        """
        cd {current_path}/results/00-prepare/assemblies/pear/
        cat *relabel_notSingleton.fasta > {params.outputname}.fasta
        mv {params.outputname}.fasta {current_path}/results/00-prepare/
        """