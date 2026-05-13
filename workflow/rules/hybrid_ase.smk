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

        samtools index {output.bam} 2>> {log}
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
            2> {log} \
            | awk '$0 ~ /^@/ || $0 ~ /NH:i:1/' \
            | samtools sort \
                -@ {threads} \
                -o {output.bam} \
                2>> {log}

        samtools index {output.bam} 2>> {log}

        echo "Created unique-only BAM from {input.bam}" >> {log}
        """


rule featurecounts_concat_multi:
    input:
        bam=expand(
            "results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
            sample_name=samples.index
        ),
        gtf=config["ref_concat"]["annotation"]["gtf"]
    output:
        counts="results/hybrid_ase/counts/concat_featureCounts_multi_fraction.txt"
    conda:
        HYBRID_ASE_ENV
    params:
        extra=config["params"]["featurecounts_multi"]
    log:
        "logs/hybrid_ase/featureCounts_concat_multi.log"
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


rule featurecounts_concat_unique:
    input:
        bam=expand(
            "results/hybrid_ase/aligned_unique/{sample_name}.concat.unique.sorted.bam",
            sample_name=samples.index
        ),
        gtf=config["ref_concat"]["annotation"]["gtf"]
    output:
        counts="results/hybrid_ase/counts/concat_featureCounts_unique.txt"
    conda:
        HYBRID_ASE_ENV
    params:
        extra=config["params"]["featurecounts_unique"]
    log:
        "logs/hybrid_ase/featureCounts_concat_unique.log"
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


'''
SNP-based ASE rules are paused for now.

rule call_parental_snps_concat:
    ...

rule filter_diagnostic_snps:
    ...

rule hybrid_snp_pileup:
    ...
'''


rule make_bigwigs_hybrid_multi:
    input:
        bam="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
        bai="results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam.bai"
    output:
        "results/hybrid_ase/bigwigs/multimapper_inclusive/{sample_name}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_ind"]
    log:
        "logs/hybrid_ase/bigwigs_multi/{sample_name}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/bigwigs/multimapper_inclusive logs/hybrid_ase/bigwigs_multi

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule make_bigwigs_hybrid_unique:
    input:
        bam="results/hybrid_ase/aligned_unique/{sample_name}.concat.unique.sorted.bam",
        bai="results/hybrid_ase/aligned_unique/{sample_name}.concat.unique.sorted.bam.bai"
    output:
        "results/hybrid_ase/bigwigs/unique_only/{sample_name}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_ind"]
    log:
        "logs/hybrid_ase/bigwigs_unique/{sample_name}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/bigwigs/unique_only logs/hybrid_ase/bigwigs_unique

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


HYBRID_SAMPLE_GROUPS = sorted(samples["condition"].unique())


def get_hybrid_multi_bams_by_group(wildcards):
    group_samples = samples.query("condition == @wildcards.sample_group").index
    return expand(
        "results/hybrid_ase/aligned/{sample_name}.concat.sorted.bam",
        sample_name=group_samples
    )


def get_hybrid_unique_bams_by_group(wildcards):
    group_samples = samples.query("condition == @wildcards.sample_group").index
    return expand(
        "results/hybrid_ase/aligned_unique/{sample_name}.concat.unique.sorted.bam",
        sample_name=group_samples
    )


rule merge_hybrid_multi_bam:
    input:
        get_hybrid_multi_bams_by_group
    output:
        bam="results/hybrid_ase/aligned_merged/multimapper_inclusive/{sample_group}.bam",
        bai="results/hybrid_ase/aligned_merged/multimapper_inclusive/{sample_group}.bam.bai"
    conda:
        HYBRID_ASE_ENV
    log:
        "logs/hybrid_ase/merge_multi/{sample_group}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/aligned_merged/multimapper_inclusive logs/hybrid_ase/merge_multi

        samtools merge \
            -@ {threads} \
            -f \
            {output.bam} \
            {input} \
            2> {log}

        samtools index {output.bam} 2>> {log}
        """


rule merge_hybrid_unique_bam:
    input:
        get_hybrid_unique_bams_by_group
    output:
        bam="results/hybrid_ase/aligned_merged/unique_only/{sample_group}.bam",
        bai="results/hybrid_ase/aligned_merged/unique_only/{sample_group}.bam.bai"
    conda:
        HYBRID_ASE_ENV
    log:
        "logs/hybrid_ase/merge_unique/{sample_group}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/aligned_merged/unique_only logs/hybrid_ase/merge_unique

        samtools merge \
            -@ {threads} \
            -f \
            {output.bam} \
            {input} \
            2> {log}

        samtools index {output.bam} 2>> {log}
        """


rule make_bigwigs_hybrid_multi_merged:
    input:
        bam="results/hybrid_ase/aligned_merged/multimapper_inclusive/{sample_group}.bam",
        bai="results/hybrid_ase/aligned_merged/multimapper_inclusive/{sample_group}.bam.bai"
    output:
        "results/hybrid_ase/bigwigs/multimapper_inclusive_merged/{sample_group}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_merged"]
    log:
        "logs/hybrid_ase/bigwigs_multi_merged/{sample_group}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/bigwigs/multimapper_inclusive_merged logs/hybrid_ase/bigwigs_multi_merged

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule make_bigwigs_hybrid_unique_merged:
    input:
        bam="results/hybrid_ase/aligned_merged/unique_only/{sample_group}.bam",
        bai="results/hybrid_ase/aligned_merged/unique_only/{sample_group}.bam.bai"
    output:
        "results/hybrid_ase/bigwigs/unique_only_merged/{sample_group}.bw"
    conda:
        "../envs/deeptools.yaml"
    params:
        extra=config["params"]["bigwigs_merged"]
    log:
        "logs/hybrid_ase/bigwigs_unique_merged/{sample_group}.log"
    threads: 8
    shell:
        r"""
        mkdir -p results/hybrid_ase/bigwigs/unique_only_merged logs/hybrid_ase/bigwigs_unique_merged

        bamCoverage \
            --bam {input.bam} \
            -o {output} \
            -p {threads} \
            {params.extra} \
            > {log} 2>&1
        """


rule zscore_normalize_hybrid_multi_ind_bigwigs:
    input:
        "results/hybrid_ase/bigwigs/multimapper_inclusive/{sample_name}.bw"
    output:
        "results/hybrid_ase/bigwigs_zscore/multimapper_inclusive/individual/{sample_name}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"


rule zscore_normalize_hybrid_unique_ind_bigwigs:
    input:
        "results/hybrid_ase/bigwigs/unique_only/{sample_name}.bw"
    output:
        "results/hybrid_ase/bigwigs_zscore/unique_only/individual/{sample_name}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"


rule zscore_normalize_hybrid_multi_merged_bigwigs:
    input:
        "results/hybrid_ase/bigwigs/multimapper_inclusive_merged/{sample_group}.bw"
    output:
        "results/hybrid_ase/bigwigs_zscore/multimapper_inclusive/merged/{sample_group}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"


rule zscore_normalize_hybrid_unique_merged_bigwigs:
    input:
        "results/hybrid_ase/bigwigs/unique_only_merged/{sample_group}.bw"
    output:
        "results/hybrid_ase/bigwigs_zscore/unique_only/merged/{sample_group}.bw"
    resources:
        mem_mb=32000,
        high_mem=1
    conda:
        "../envs/zscore_normalize_bw.yaml"
    script:
        "../scripts/zscore_normalize_bw.R"