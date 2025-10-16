

rule bwaindex:
    input:
        genome=genome_path,
    output:
        idx=idx,
    params:
        bwa=config["map"]["mapper"],
    threads: 1  #Only affects bwa-meme
    log:
        f"logs/bwa-memx_index/{assembly}.log",
    cache: True
    wrapper:
        "v4.6.0/bio/bwa-memx/index"


rule chromap_index:
    input:
        genome=genome_path,
    output:
        idx=multiext(genome_path, ".chromap.index"),
    log:
        f"logs/chromap_index/{assembly}.log",
    conda:
        "envs/chromap.yml"
    shell:
        r"chromap -i -r {input.genome} -o {output.idx} >{log} 2>&1"



rule map_chunks_bwa:
    input:
        reads=lambda wildcards: (
            [
                f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}_trimmed.fastq.gz",
                f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}_trimmed.fastq.gz",
            ]
            if config["map"]["trim_options"]
            else [
                f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}.fastq.gz",
                f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}.fastq.gz",
            ]
        ),
        reference=genome_path,
        idx=idx,
    params:
        bwa=config["map"]["mapper"],
        extra="-SP5M",
        sort="none",
        dedup="none",
    threads: 4
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.bam",
    log:
        "logs/bwa_memx/{library}.{run}.{chunk_id}.log",
    benchmark:
        "benchmarks/bwa_memx/{library}.{run}.{chunk_id}.tsv"
    wrapper:
        "v4.6.0/bio/bwa-memx/mem"


rule map_chunks_chromap:
    input:
        reads=lambda wildcards: (
            [
                f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}_trimmed.fastq.gz",
                f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}_trimmed.fastq.gz",
            ]
            if config["map"]["trim_options"]
            else [
                f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}.fastq.gz",
                f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}.fastq.gz",
            ]
        ),
        reference=genome_path,
        idx=multiext(genome_path, ".chromap.index"),
    params:
        extra=config["map"].get("mapping_options", ""),
    threads: 8
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.{assembly}.pairs.gz",
    log:
        "logs/chromap/{library}.{run}.{chunk_id}.log",
    benchmark:
        "benchmarks/chromap/{library}.{run}.{chunk_id}.tsv"
    conda:
        "envs/chromap.yml"
    shell:
        # chromap doesn't output gzip files, so we need to pipe it to bgzip
        # It doesn't work with low memory mode, so we can't use the hic preset
        # Hence I provide all arguments manually except for --low-mem
        r"""
        chromap -e 4 -q 1 --split-alignment --pairs -x {input.idx} -r {input.reference} \
        -t {threads} {params.extra} \
        -1 {input.reads[0]} -2 {input.reads[1]} -o /dev/stdout 2>{log} | \
        bgzip > {output} \
        """

