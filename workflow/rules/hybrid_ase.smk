import os


rule hisat2_index_concat:
    input:
        fasta=config["ref_concat"]["genome"]["fasta"]
    output:
        expand(
            config["ref_concat"]["hisat2_index"]["prefix"] + ".{i}.ht2",
            i=range(1, 9)
        )
    params:
        prefix=config["ref_concat"]["hisat2_index"]["prefix"],
        index_dir=lambda wc: os.path.dirname(
            config["ref_concat"]["hisat2_index"]["prefix"]
        )
    log:
        "logs/hybrid_ase/hisat2_index_concat.log"
    threads: 8
    shell:
        r"""
        mkdir -p {params.index_dir} logs/hybrid_ase

        hisat2-build \
            -p {threads} \
            {input.fasta} \
            {params.prefix} \
            > {log} 2>&1
        """
