

if config["map"]["mapper"] in ["bwa-mem", "bwa-mem2", "bwa-meme"]:

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

elif config["map"]["mapper"] == "bowtie2":

    rule bowtie2_index:
        input:
            reference=bowtie_index_path,
        output:
            idx=idx,
        log:
            f"logs/bowtie2_index/{assembly}.log",
        params:
            extra="",
        threads: 8
        cache: True
        wrapper:
            "v3.9.0/bio/bowtie2/build"

elif config["map"]["mapper"] == "chromap":

    rule chromap_index:
        input:
            genome=genome_path,
        output:
            idx=multiext(genome_path, ".chromap.index"),
        log:
            f"logs/chromap_index/{assembly}.log",
        conda:
            "../envs/chromap.yml"
        shell:
            r"chromap -i -r {input.genome} -o {output.idx} >{log} 2>&1"


if config["map"]["mapper"] in ["bwa-mem", "bwa-mem2", "bwa-meme"]:

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


if config["map"]["mapper"] == "bowtie2":

    # Bowtie2 rescue mapping workflow - Step 1: Global alignment for R1
    rule bowtie2_global_align_R1:
        input:
            sample=lambda wildcards: (
                [f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}_trimmed.fastq.gz"]
                if config["map"]["trim_options"]
                else [f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}.fastq.gz"]
            ),
 #           reference=bowtie_index_path,
            idx=idx,
        output:
            bam=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.global.bam"),
            unaligned=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.global.unmap.fastq"),
        params:
            extra=config["map"].get("rescue_options", {}).get("global_extra", "--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder"),
        threads: 8
        log:
            "logs/bowtie2_global/{library}.{run}.{chunk_id}_R1.log",
        benchmark:
            "benchmarks/bowtie2_global/{library}.{run}.{chunk_id}_R1.tsv"
        wrapper:
            "v3.9.0/bio/bowtie2/align"

    # Step 1b: Global alignment for R2
    rule bowtie2_global_align_R2:
        input:
            sample=lambda wildcards: (
                [f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}_trimmed.fastq.gz"]
                if config["map"]["trim_options"]
                else [f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}.fastq.gz"]
            ),
  #          reference=bowtie_index_path,
            idx=idx,
        output:
            bam=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.global.bam"),
            unaligned=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.global.unmap.fastq"),
        params:
            extra=config["map"].get("rescue_options", {}).get("global_extra", "--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder"),
        threads: 8
        log:
            "logs/bowtie2_global/{library}.{run}.{chunk_id}_R2.log",
        benchmark:
            "benchmarks/bowtie2_global/{library}.{run}.{chunk_id}_R2.tsv"
        wrapper:
            "v3.9.0/bio/bowtie2/align"

    # Step 2: Cutsite trimming for R1 unmapped reads
    rule cutsite_trim_R1:
        input:
            fastq=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.global.unmap.fastq",
        output:
            fastq=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.trimmed.fastq"),
        params:
            cutsite=config["map"].get("cutsite", "GATCGATC"),
        log:
            "logs/cutsite_trim/{library}.{run}.{chunk_id}_R1.log",
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            "workflow/scripts/cutsite_trimming --fastq {input.fastq} --cutsite {params.cutsite} --out {output.fastq} >{log} 2>&1"

    # Step 2b: Cutsite trimming for R2 unmapped reads
    rule cutsite_trim_R2:
        input:
            fastq=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.global.unmap.fastq",
        output:
            fastq=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.trimmed.fastq"),
        params:
            cutsite=config["map"].get("cutsite", "GATCGATC"),
        log:
            "logs/cutsite_trim/{library}.{run}.{chunk_id}_R2.log",
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            "workflow/scripts/cutsite_trimming --fastq {input.fastq} --cutsite {params.cutsite} --out {output.fastq} >{log} 2>&1"

    # Step 3: Local alignment for R1 trimmed reads
    rule bowtie2_local_align_R1:
        input:
            sample=[f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.trimmed.fastq"],
  #          reference=bowtie_index_path,
            idx=idx,
        output:
            bam=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.local.bam"),
        params:
            extra=config["map"].get("rescue_options", {}).get("local_extra", "--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder"),
        threads: 8
        log:
            "logs/bowtie2_local/{library}.{run}.{chunk_id}_R1.log",
        benchmark:
            "benchmarks/bowtie2_local/{library}.{run}.{chunk_id}_R1.tsv"
        wrapper:
            "v3.9.0/bio/bowtie2/align"

    # Step 3b: Local alignment for R2 trimmed reads
    rule bowtie2_local_align_R2:
        input:
            sample=[f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.trimmed.fastq"],
 #           reference=bowtie_index_path,
            idx=idx,
        output:
            bam=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.local.bam"),
        params:
            extra=config["map"].get("rescue_options", {}).get("local_extra", "--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder"),
        threads: 8
        log:
            "logs/bowtie2_local/{library}.{run}.{chunk_id}_R2.log",
        benchmark:
            "benchmarks/bowtie2_local/{library}.{run}.{chunk_id}_R2.tsv"
        wrapper:
            "v3.9.0/bio/bowtie2/align"

    # Step 4: Merge global and local alignments for R1
    rule merge_global_local_R1:
        input:
            global_bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.global.bam",
            local_bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.local.bam",
        output:
            merged=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.merged.bam"),
        threads: 4
        log:
            "logs/merge_global_local/{library}.{run}.{chunk_id}_R1.log",
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            r"""
            samtools merge -@ {threads} -n -f {output.merged} {input.global_bam} {input.local_bam} 2>{log} && \
            samtools sort -@ {threads} -m 2G -n -o {output.merged}.tmp {output.merged} 2>>{log} && \
            mv {output.merged}.tmp {output.merged}
            """

    # Step 4b: Merge global and local alignments for R2
    rule merge_global_local_R2:
        input:
            global_bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.global.bam",
            local_bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.local.bam",
        output:
            merged=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.merged.bam"),
        threads: 4
        log:
            "logs/merge_global_local/{library}.{run}.{chunk_id}_R2.log",
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            r"""
            samtools merge -@ {threads} -n -f {output.merged} {input.global_bam} {input.local_bam} 2>{log} && \
            samtools sort -@ {threads} -m 2G -n -o {output.merged}.tmp {output.merged} 2>>{log} && \
            mv {output.merged}.tmp {output.merged}
            """

    # Step 5: Pair R1 and R2 reads using mergeSAM.py
    rule pair_rescue_reads:
        input:
            r1=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R1.merged.bam",
            r2=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_R2.merged.bam",
        output:
            paired=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.bam",
            stats=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.pairstat",
        params:
            min_mapq=config["map"].get("rescue_options", {}).get("min_mapq", 10),
        threads: 4
        log:
            "logs/pair_rescue/{library}.{run}.{chunk_id}.log",
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            "python workflow/scripts/mergeSAM.py -v -t -q {params.min_mapq} -f {input.r1} -r {input.r2} -o {output.paired} >{log} 2>&1"


if config["map"]["mapper"] == "chromap":

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
            "../envs/chromap.yml"
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