

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

    # Bowtie2 rescue mapping workflow - Step 1: Global alignment
    rule bowtie2_global_align:
        input:
            sample=lambda wildcards: (
                [f"{processed_fastqs_folder}/{{library}}/{{run}}/{wildcards.side}.{{chunk_id}}_trimmed.fastq.gz"]
                if config["map"]["trim_options"]
                else [f"{processed_fastqs_folder}/{{library}}/{{run}}/{wildcards.side}.{{chunk_id}}.fastq.gz"]
            ),
            idx=idx,
        output:
            bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.global.bam",
            unaligned=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.global.unmap.fastq"),
        params:
            extra=config["map"].get("rescue_options", {}).get("global_extra", "--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder"),
        threads: 8
        log:
            "logs/bowtie2_global/{library}.{run}.{chunk_id}_{side}.log",
        benchmark:
            "benchmarks/bowtie2_global/{library}.{run}.{chunk_id}_{side}.tsv"
        wildcard_constraints:
            side="[12]"
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            r"""
            # Find bowtie2 index base (strip .1.bt2, .2.bt2, .rev.1.bt2, .rev.2.bt2, etc.)
            idx_files=({input.idx})
            index_prefix="${{idx_files[0]%%.*}}"
            (bowtie2 {params.extra} \
                -p {threads} \
                -x "${{index_prefix}}" \
                -U {input.sample[0]} \
                --un {output.unaligned} \
                2> {log}) \
            | samtools view -F 4 -@ {threads} -bS - \
            > {output.bam}
            """


    # Step 2: Cutsite trimming for unmapped reads
    rule cutsite_trim:
        input:
            fastq=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.global.unmap.fastq",
        output:
            fastq=temp(f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.trimmed.fastq"),
        params:
            cutsite=lambda wildcards: (
                SAMPLE_METADATA.get(wildcards.library, {})
                               .get(wildcards.run, {})
                               .get('ligation_site', config["map"].get("cutsite", "GATCGATC"))
            ),
        log:
            "logs/cutsite_trim/{library}.{run}.{chunk_id}_{side}.log",
        wildcard_constraints:
            side="[12]"
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            "workflow/scripts/cutsite_trimming --fastq {input.fastq} --cutsite {params.cutsite} --out {output.fastq} >{log} 2>&1"

    # Step 3: Local alignment for trimmed reads
    rule bowtie2_local_align:
        input:
            sample=[f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.trimmed.fastq"],
            idx=idx,
        output:
            bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.local.bam",
        params:
            extra=config["map"].get("rescue_options", {}).get("local_extra", "--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder"),
        threads: 8
        log:
            "logs/bowtie2_local/{library}.{run}.{chunk_id}_{side}.log",
        benchmark:
            "benchmarks/bowtie2_local/{library}.{run}.{chunk_id}_{side}.tsv"
        wildcard_constraints:
            side="[12]"
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            r"""
            # Find bowtie2 index base (strip .1.bt2, .2.bt2, .rev.1.bt2, .rev.2.bt2, etc.)
            idx_files=({input.idx})
            index_prefix="${{idx_files[0]%%.*}}"
            (bowtie2 {params.extra} \
                -p {threads} \
                -x "${{index_prefix}}" \
                -U {input.sample[0]} \
                2> {log}) \
            | samtools view -@ {threads} -bS - \
            > {output.bam}
            """
    
    def get_merge_inputs(wildcards):
        """
        Dynamically determine the input files for merging.
        If skip_ligation is True for the sample, only return the global_bam.
        Otherwise, return both global_bam and local_bam.
        """
        # Check the skip_ligation flag from the metadata
        skip_ligation = (
            SAMPLE_METADATA.get(wildcards.library, {})
                        .get(wildcards.run, {})
                        .get('skip_ligation', False)
        )

        # Start with the global_bam, which is always required
        input_files = [
            f"{mapped_parsed_sorted_chunks_folder}/{wildcards.library}/{wildcards.run}/{wildcards.chunk_id}_{wildcards.side}.global.bam"
        ]

        # If we are NOT skipping ligation, add the local_bam
        if not skip_ligation:
            input_files.append(
                f"{mapped_parsed_sorted_chunks_folder}/{wildcards.library}/{wildcards.run}/{wildcards.chunk_id}_{wildcards.side}.local.bam"
            )

        # Return the list of files
        return input_files


    # Step 4: Merge global and local alignments
    rule merge_global_local:
        input:
            get_merge_inputs 
        output:
            merged=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.merged.bam",
            mapstat=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_{{side}}.mapstat",
        threads: 4
        log:
            "logs/merge_global_local/{library}.{run}.{chunk_id}_{side}.log",
        wildcard_constraints:
            side="[12]"
        params:
            num_inputs=lambda i: len(i),
        conda:
            "../envs/bowtie2_rescue.yml"
        shell:
            r"""
            # samtools merge (one or two inputs)
            samtools merge -@ {threads} -n -f {output.merged} {input} 2>>{log}

            # sort (name sort)
            samtools sort -@ {threads} -m 2G -n -o {output.merged}.tmp {output.merged} 2>>{log} && \
            mv {output.merged}.tmp {output.merged}

            # 1. Assign the Snakemake path to a new shell variable
            merged_bam="{output.merged}"

            # 2. Use double braces {{...}} to tell Snakemake to ignore this line.
            #    The shell will now correctly perform the string trimming.
            prefix="${{merged_bam%.merged.bam}}"
            tag="{wildcards.side}"
            # Set variables for bam1 (global BAM) and bam2 (local BAM if present)
            inputs=({input})
            bam1="${{inputs[0]}}"
            if (( ${{#inputs[@]}} > 1 )); then
                bam2="${{inputs[1]}}"
            else
                bam2=""
            fi

            echo "## ${{prefix}}" > {output.mapstat}
            printf "total_${{tag}}\t" >> {output.mapstat}
            samtools view -c {output.merged} >> {output.mapstat}
            printf "mapped_${{tag}}\t" >> {output.mapstat}
            samtools view -c -F 4 {output.merged} >> {output.mapstat}
            printf "global_${{tag}}\t" >> {output.mapstat}
            samtools view -c -F 4 "$bam1" >> {output.mapstat}

            # If local bam is present, count mapped for local; else, print 0
            printf "local_${{tag}}\t"  >> {output.mapstat}
            if [[ -n "$bam2" ]]; then
                samtools view -c -F 4 "$bam2" >> {output.mapstat}
            else
                echo 0 >> {output.mapstat}
            fi
            """

    # Step 5: Pair R1 and R2 reads using mergeSAM.py
    rule pair_rescue_reads:
        input:
            r1=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_1.merged.bam",
            r2=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}_2.merged.bam",
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