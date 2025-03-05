rule concatenate_singletons:
    input:
        assembled_fastq = os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.assembled.fastq"),
        unassembled_forward = os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.unassembled.forward.fastq"),
        unassembled_reverse = os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.unassembled.reverse.fastq"),
        singleton_truncated = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.singleton.truncated")
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}_withSingleton.fastq")
    shell:
        """
        cat {input.assembled_fastq} {input.unassembled_forward} {input.unassembled_reverse} {input.singleton_truncated} > {output}
        sed -i 's/ /_/g' {output}
        """

rule run_prinseq_singleton:
    input:
        assembled_singleton = os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}_withSingleton.fastq")
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}.assembled.withSingleton.fasta"),
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "prinseq_singleton", "{sample}.log")
    singularity:
        "docker://itvdsbioinfo/pimba_prinseq:v0.20.4"
    shell:
        """
        export LANG=C
        export LC_ALL=C
        prinseq-lite.pl -fastq {input.assembled_singleton} -out_format 1 -seq_id Seq -out_good {current_path}/results/00-prepare/assemblies/pear/{wildcards.sample}.assembled.withSingleton \
        > {log} 2>&1
        """

rule run_relabel_singleton:
    input:
        assembled_fasta_singleton = expand(os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}.assembled.withSingleton.fasta"), sample=SAMPLES)
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}_relabel_withSingleton.fasta")
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "qiimepipe_singleton", "{sample}.log")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        cd {current_path}/results/00-prepare/assemblies/pear/
        python3.6 /qiimepipe/relabelReads-v2.py {wildcards.sample}.assembled.withSingleton.fasta . \
        > {log} 2>&1
        mv {wildcards.sample}_relabel.fasta {wildcards.sample}_relabel_withSingleton.fasta
        """

rule concatenate_and_rename_singletons:
    input:
        expand(os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}_relabel_withSingleton.fasta"), sample=SAMPLES)
    output:
        os.path.join(current_path, "results", "00-prepare", config["outputprepare"] + "_withSingleton.fasta")
    params:
        outputname = config["outputprepare"]
    shell:
        """
        cd {current_path}/results/00-prepare/assemblies/pear/
        cat *_relabel_withSingleton.fasta > {params.outputname}_withSingleton.fasta
        mv {params.outputname}_withSingleton.fasta {current_path}/results/00-prepare/
        rm -r {current_path}/results/00-prepare/assemblies/
        """