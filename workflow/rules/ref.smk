rule get_ref_genome:
    output:
        temp("resources/ref_genome.fasta")
    log:
        "logs/get_ref_genome.log"
    conda:
        "../envs/curl.yaml"
    params:
        link=config["ref_genome"]["link"]
    cache: True
    shell:
        r"""
        curl {params.link} > {output}.gz 2> {log}
        unpigz {output}.gz 2>> {log}
        """


rule get_genome_annotation:
    output:
        "resources/dm6.ncbiRefSeq.gtf"
    log:
        "logs/get_genome_annotation.log"
    conda:
        "../envs/curl.yaml"
    params:
        link=config["ref_annotation"]["link"]
    cache: True
    shell:
        r"""
        curl -L {params.link} 2> {log} | gunzip -c > {output}
        """


if config["use_spikeIn"]:

    rule get_spikeIn_genome:
        output:
            temp("resources/spikeIn_genome.fasta")
        log:
            "logs/get_spikeIn_genome.log"
        conda:
            "../envs/curl.yaml"
        params:
            link=config["spikeIn_genome"]["link"]
        cache: True
        shell:
            r"""
            curl {params.link} > {output}.gz 2> {log}
            unpigz {output}.gz 2>> {log}
            """

    rule combine_genomes:
        input:
            ref="resources/ref_genome.fasta",
            spikeIn="resources/spikeIn_genome.fasta"
        output:
            temp("resources/genome.fasta.gz")
        log:
            "logs/combine_genomes.log"
        conda:
            "../envs/seqkit.yaml"
        cache: True
        shell:
            r"""
            # Add a prefix to spike-in contig names, then concatenate with ref genome
            seqkit replace -p "(.+)" -r "spikeIn_$1" -o resources/tmp_spikeIn.fasta {input.spikeIn} 2> {log}

            # IMPORTANT: actually gzip the combined fasta to match the .gz filename
            cat {input.ref} resources/tmp_spikeIn.fasta | gzip -c > {output} 2>> {log}

            rm resources/tmp_spikeIn.fasta
            """

else:

    rule rename_genome:
        input:
            "resources/ref_genome.fasta"
        output:
            temp("resources/genome.fasta.gz")
        log:
            "logs/rename_genome.log"
        cache: True
        shell:
            r"""
            mv {input} resources/genome.fasta 2> {log}
            gzip -f resources/genome.fasta 2>> {log}
            """


if config["filter_chroms"]:

    rule define_keep_chroms:
        input:
            genome="resources/genome.fasta.gz",
            keep_chroms=config["keep_chroms"]
        output:
            "resources/keep_chroms.bed"
        log:
            "logs/define_keep_chroms.log"
        conda:
            "../envs/seqkit.yaml"
        cache: True
        shell:
            r"""
            seqkit grep -f {input.keep_chroms} {input.genome} \
              | seqkit fx2tab -nil \
              | awk -v OFS='\t' '{{print $1, 1, $2}}' > {output}
            """


rule hisat2_index:
    input:
        fasta="resources/genome.fasta.gz"
    output:
        "resources/hisat2_index.1.ht2",
        "resources/hisat2_index.2.ht2",
        "resources/hisat2_index.3.ht2",
        "resources/hisat2_index.4.ht2",
        "resources/hisat2_index.5.ht2",
        "resources/hisat2_index.6.ht2",
        "resources/hisat2_index.7.ht2",
        "resources/hisat2_index.8.ht2"
    log:
        "logs/hisat2_index_genome.log"
    conda:
        "../envs/hisat2.yaml"
    params:
        prefix="resources/hisat2_index",
        extra=config["params"]["hisat2_index"]
    threads: 8
    shell:
        r"""
        gunzip -c {input.fasta} > resources/genome.fasta
        hisat2-build {params.extra} -p {threads} resources/genome.fasta {params.prefix} > {log} 2>&1
        """
