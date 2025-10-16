
rule download_fastqs:
    output:
        fastq1=f"{downloaded_fastqs_folder}/{{library}}.{{run}}.1.fastq.gz",
        fastq2=f"{downloaded_fastqs_folder}/{{library}}.{{run}}.2.fastq.gz",
    log:
        "logs/download_fastqs/{library}.{run}.log",
    threads: 8
    conda:
        "../envs/download_fastqs.yml"
    params:
        bgzip_threads=lambda wildcards, threads: max(1, (threads - 2) / 2),
        fastq_files=lambda wildcards: LIBRARY_RUN_FASTQS[wildcards.library][
            wildcards.run
        ],
    script:
        "../scripts/download_split_fastq.py"


# Chunking command, depending on whether chunking is required - simply copying the files
# if no chunking needed, or splitting them into chunks of required size
if config["map"]["chunksize"] > 0:
    chunk_command = """
        mkdir -p {output};
        zcat {input.fastq1} | split -l {params.chunksize} -d \
            --filter 'bgzip -c -@ {threads} > $FILE.fastq.gz' - \
            {output}/1.
        zcat {input.fastq2} | split -l {params.chunksize} -d \
            --filter 'bgzip -c -@ {threads} > $FILE.fastq.gz' - \
            {output}/2.
        """
else:
    chunk_command = """
        mkdir -p {output}
        cp {input.fastq1} {output}/1.00.fastq.gz
        cp {input.fastq2} {output}/2.00.fastq.gz
        """


rule trim:
    input:
        sample=[
            f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}.fastq.gz",
            f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}.fastq.gz",
        ],
    log:
        "logs/fastp/{library}.{run}.{chunk_id}.log",
    params:
        extra=config["map"]["trim_options"],
    output:
        # Would be better with a pipe, but it causes weird freezing of the pipeline
        trimmed=[
            temp(
                f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}_trimmed.fastq.gz"
            ),
            temp(
                f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}_trimmed.fastq.gz"
            ),
        ],
        json=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.fastp.json",
        html=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.fastp.html",
    wrapper:
        "v4.6.0/bio/fastp"



rule fastqc:
    input:
        lambda wildcards: (
            f"downloaded_fastqs/{wildcards.library}.{wildcards.run}.{wildcards.side}.fastq.gz"
            if needs_downloading(
                LIBRARY_RUN_FASTQS[wildcards.library][wildcards.run],
                int(wildcards.side) - 1,
            )
            else LIBRARY_RUN_FASTQS[wildcards.library][wildcards.run][
                int(wildcards.side) - 1
            ]
        ),
    output:
        html=f"{fastqc_folder}/{{library}}.{{run}}.{{side}}_fastqc.html",
        zip=f"{fastqc_folder}/{{library}}.{{run}}.{{side}}_fastqc.zip",
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/{library}.{run}.{side}.log",
    resources:
        mem_mb = 1024,
    benchmark:
        "benchmarks/fastqc/{library}.{run}.{side}.tsv",
    conda:
        "../envs/fastqc.yml",
    threads: 1
    script:
        "../scripts/fastqc.py"