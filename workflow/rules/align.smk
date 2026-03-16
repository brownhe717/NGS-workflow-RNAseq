rule get_sra_se:
	output:
		temp("data/sra/se/{accession}.fastq.gz"),
	conda:
		"../envs/sratools.yaml"
	log:
		"logs/get_sra/{accession}.log",
	threads: 2
	script:
			"../scripts/fasterq-dump.py"
                
rule get_sra_pe:
	output:
		temp("data/sra/pe/{accession}_1.fastq.gz"),
		temp("data/sra/pe/{accession}_2.fastq.gz"),
	conda:
		"../envs/sratools.yaml"
	log:
		"logs/get_sra/{accession}.log",
	threads: 2
	script:
			"../scripts/fasterq-dump.py"
        
rule merge_fastqs:
    input:
        get_fq_merge,
    output:
        temp("data/merged/{sample}_{read}.fastq.gz"),
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|1|2",
    shell:
        "cat {input} > {output} 2> {log}"
        
rule hisat2_align:
    input:
        reads=get_hisat2_input,
        idx=rules.hisat2_index.output
    output:
        temp("results/aligned_reads/mapped/{sample}.bam")
    log:
        "logs/hisat2_align/{sample}.log"
    params:
        idx="resources/hisat2_index",
        extra=config["params"]["hisat2_align"]
    threads: 8
    wrapper:
        "v1.3.2/bio/hisat2/align"

rule samtools_sort:
    input:
       "results/aligned_reads/mapped/{sample}.bam"
    output:
        temp("results/aligned_reads/sorted/{sample}.bam")
    log:
        "logs/samtools_sort/{sample}.log"
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v1.1.0/bio/samtools/sort"
        
rule samtools_index_aligned:
    input:
        "results/aligned_reads/sorted/{sample}.bam"
    output:
        temp("results/aligned_reads/sorted/{sample}.bam.bai")
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"
