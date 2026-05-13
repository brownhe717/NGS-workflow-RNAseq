import os


HYBRID_ASE_ENV = "../envs/hybrid_ase.yaml"


rule hisat2_index_concat:
    input:
        fasta=config["ref_concat"]["genome"]["fasta"]
    output:
        expand(
            config["ref_concat"]["hisat2_index"]["prefix"] + ".{i}.ht2",
            i=range(1, 9)
        )
    conda:
        HYBRID_ASE_ENV
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


rule hisat2_align_concat:
    input:
        idx=expand(
            config["ref_concat"]["hisat2_index"]["prefix"] + ".{i}.ht2",
            i=range(1, 9)
        ),
        fq1=lambda wc: units.loc[wc.sample_name, "fq1"],
        fq2=lambda wc: units.loc[wc.sample_name, "fq2"]
    output:
        bam="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
        bai="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam.bai"
    conda:
        HYBRID_ASE_ENV
    params:
        idx_prefix=config["ref_concat"]["hisat2_index"]["prefix"],
        extra=lambda wc: (
            config["params"]["hisat2_align"]
            + f" -k {config['hybrid_ase']['hisat2_k']}"
        )
    log:
        "logs/hybrid_ase/hisat2_align/{sample_name}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/aligned logs/hybrid_ase/hisat2_align

        hisat2 \
            -x {params.idx_prefix} \
            -1 {input.fq1} \
            -2 {input.fq2} \
            -p {threads} \
            {params.extra} \
            2> {log} | \
        samtools sort \
            -@ {threads} \
            -o {output.bam}

        samtools index {output.bam}
        """

rule filter_unique_concat_bam:
    input:
        bam="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
        bai="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam.bai"
    output:
        bam="results/hybrid_ase/aligned_unique/{sample_name}.concat.unique.sorted.bam",
        bai="results/hybrid_ase/aligned_unique/{sample_name}.concat.unique.sorted.bam.bai"
    conda:
        HYBRID_ASE_ENV
    params:
        min_mapq=config["hybrid_ase"]["unique_bam_min_mapq"]
    log:
        "logs/hybrid_ase/filter_unique/{sample_name}.log"
    threads: 4
    shell:
        r"""
        mkdir -p results/hybrid_ase/aligned_unique logs/hybrid_ase/filter_unique

        samtools view \
            -h \
            -q {params.min_mapq} \
            {input.bam} \
            | awk '$0 ~ /^@/ || $0 ~ /NH:i:1/' \
            | samtools sort \
                -@ {threads} \
                -o {output.bam}

        samtools index {output.bam}

        echo "Created unique-only BAM from {input.bam}" > {log}
        """

rule featurecounts_concat:
    input:
        bam=expand(
            "results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
            sample_name=samples.index
        ),
        gtf=config["ref_concat"]["annotation"]["gtf"]
    output:
        counts="results/hybrid_ase/counts/concat_featureCounts.txt"
    conda:
        HYBRID_ASE_ENV
    params:
        extra=config["params"]["featurecounts"]
    log:
        "logs/hybrid_ase/featureCounts_concat.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/counts logs/hybrid_ase

        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            {params.extra} \
            {input.bam} \
            > {log} 2>&1
        """


rule call_parental_snps_concat:
    input:
        mel=expand(
            "results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
            sample_name=samples.query("ase_role == 'mel_parent'").index
        ),
        sim=expand(
            "results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
            sample_name=samples.query("ase_role == 'sim_parent'").index
        ),
        ref=config["ref_concat"]["genome"]["fasta"]
    output:
        vcf="results/hybrid_ase/snps/parental_raw.vcf.gz",
        tbi="results/hybrid_ase/snps/parental_raw.vcf.gz.tbi"
    conda:
        HYBRID_ASE_ENV
    log:
        "logs/hybrid_ase/call_parental_snps_concat.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/snps logs/hybrid_ase

        bcftools mpileup \
            -Ou \
            -f {input.ref} \
            -a FORMAT/DP,FORMAT/AD \
            {input.mel} \
            {input.sim} \
            2> {log} | \
        bcftools call \
            -mv \
            -Oz \
            -o {output.vcf}

        tabix -p vcf {output.vcf}
        """


rule filter_diagnostic_snps:
    input:
        vcf="results/hybrid_ase/snps/parental_raw.vcf.gz"
    output:
        vcf="results/hybrid_ase/snps/diagnostic_mel_sim_snps.vcf.gz",
        tbi="results/hybrid_ase/snps/diagnostic_mel_sim_snps.vcf.gz.tbi"
    conda:
        HYBRID_ASE_ENV
    params:
        min_depth=config["hybrid_ase"]["min_snp_depth"]
    log:
        "logs/hybrid_ase/filter_diagnostic_snps.log"
    shell:
        r"""
        mkdir -p results/hybrid_ase/snps logs/hybrid_ase

        bcftools view \
            -v snps \
            -m2 -M2 \
            -i 'FORMAT/DP[0]>={params.min_depth} && FORMAT/DP[1]>={params.min_depth} && ((GT[0]="0/0" && GT[1]="1/1") || (GT[0]="1/1" && GT[1]="0/0"))' \
            {input.vcf} \
            -Oz \
            -o {output.vcf} \
            2> {log}

        tabix -p vcf {output.vcf}
        """


rule hybrid_snp_pileup:
    input:
        bam="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
        vcf="results/hybrid_ase/snps/diagnostic_mel_sim_snps.vcf.gz",
        ref=config["ref_concat"]["genome"]["fasta"]
    output:
        bcf="results/hybrid_ase/allele_counts/{sample_name}.diagnostic_sites.bcf"
    conda:
        HYBRID_ASE_ENV
    log:
        "logs/hybrid_ase/allele_counts/{sample_name}.log"
    shell:
        r"""
        mkdir -p results/hybrid_ase/allele_counts logs/hybrid_ase/allele_counts

        bcftools mpileup \
            -Ou \
            -f {input.ref} \
            -R {input.vcf} \
            -a FORMAT/DP,FORMAT/AD \
            {input.bam} \
            2> {log} | \
        bcftools call \
            -m \
            -Ob \
            -o {output.bcf}
        """
