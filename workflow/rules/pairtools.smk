

rule parse_sort_chunks:
    input:
        bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.bam",
        chromsizes=chromsizes_path,
    threads: 4
    params:
        # keep_bams_command=f"| tee >(samtools view -bS > {mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.{assembly}.bam)"
        # if config["parse"]["keep_unparsed_bams"]
        # else "",
        dropsam_flag="" if config["parse"].get("make_pairsam", False) else "--drop-sam",
        dropreadid_flag=(
            "--drop-readid" if config["parse"].get("drop_readid", False) else ""
        ),
        dropseq_flag="--drop-seq" if config["parse"].get("drop_seq", True) else "",
        parsing_options=config["parse"].get("parsing_options", ""),
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.{assembly}.pairs.gz",
    benchmark:
        "benchmarks/parse_sort_chunks/{library}.{run}.{chunk_id}.tsv"
    log:
        "logs/parse_sort_chunks/{library}.{run}.{chunk_id}.log",
    shell:
        r"""
        pairtools parse {params.dropsam_flag} {params.dropreadid_flag} {params.dropseq_flag} \
        {params.parsing_options} \
        -c {input.chromsizes} {input.bam} \
        | pairtools sort --nproc {threads} -o {output} \
        >{log[0]} 2>&1
        """



# Find out the chunk ids for each library and run - since we don't know them beforehand
def get_pair_chunks(wildcards):
    checkpoint_output = checkpoints.chunk_fastq.get(**wildcards).output
    chunk_ids = glob_wildcards(
        f"{processed_fastqs_folder}/{wildcards.library}/{wildcards.run}/2.{{chunk_id}}.fastq.gz",
    ).chunk_id
    # Normalize chunk_ids: strip "_trimmed" suffix if present
    # This ensures consistent chunk_id resolution regardless of trimming state
    chunk_ids = [cid.replace("_trimmed", "") for cid in chunk_ids]
    paths = expand(
        f"{mapped_parsed_sorted_chunks_folder}/{wildcards.library}/{wildcards.run}/{{chunk_id}}.{assembly}.pairs.gz",
        chunk_id=chunk_ids,
    )
    return paths