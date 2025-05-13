rule run_overlapper:
    input:
        pair1 = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair1.truncated"),
        pair2 = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair2.truncated"),
    output:
        assembled = os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.assembled.fastq"),
        un_forward = os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.unassembled.forward.fastq"),
        un_reverse = os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.unassembled.reverse.fastq")
    params:
        minoverlap = config["minoverlap"],
        minsim = config["minsim"]
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "overlapper", "{sample}.log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "{sample}_prepare_overlapper.txt")
    shell:
        """
        mkdir -p {current_path}/results/00-prepare/assemblies/{wildcards.sample}
        cd {current_path}/results/00-prepare/assemblies/{wildcards.sample}
        python {current_path}/workflow/scripts/overlapper.py -f {input.pair1} -r {input.pair2} --mo {params.minoverlap} --ms {params.minsim} > {log} 2>&1
        mv overlapped.fastq {output.assembled}
        mv notAssembled-1.fastq {output.un_forward}
        mv notAssembled-2.fastq {output.un_reverse}
        for i in {current_path}/results/00-prepare/assemblies/{wildcards.sample}/*.assembled.fastq; do sed -i "s/ /_/g" $i; done
        """