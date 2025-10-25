

rule combine_dedup_stats:
    input:
        expand(
            f"{pairs_library_folder}/{{library}}.{assembly}.dedup.stats",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_dedup_stats_summary.tsv",
    log:
        "logs/combine_dedup_stats.log",
    run:
        import pandas as pd
        from pathlib import Path
        
        # Dictionary to store all data
        all_stats = {}
        
        # Read each stats file
        for stats_file in input:
            library_name = Path(stats_file).stem.replace(f".{assembly}.dedup", "")
            
            # Read the two-column file
            stats_dict = {}
            with open(stats_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:  # Skip empty lines
                        parts = line.split('\t')
                        if len(parts) == 2:
                            key, value = parts
                            stats_dict[key] = value
            
            all_stats[library_name] = stats_dict
        
        # Convert to DataFrame
        df = pd.DataFrame(all_stats)
        
        # Reset index to make library name a column
        #df.index.name = 'library'
        #df.reset_index(inplace=True)
        
        # Save to TSV
        df.to_csv(output.summary, sep='\t', index=True)
        
        with open(log[0], 'w') as log_file:
            log_file.write(f"Combined {len(all_stats)} library stats files\n")
            log_file.write(f"Output written to {output.summary}\n")

rule concat_mapstat:
    input:
        mapstats=expand(
            f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.mapstat",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_mapstats_summary.tsv",
    log:
        "logs/concat_mapstat/concat_mapstat.log",
    run:
        from pathlib import Path

        # Build library names based on the order of input files (parent folder name)
        library_names = [Path(f).parent.name for f in input.mapstats]
        input_files_str = " ".join(input.mapstats)
        names_str = " ".join(library_names)

        shell(
            "python workflow/scripts/concat_stat.py "
            "-i {input_files_str} "
            "-n {names_str} "
            "-o {output.summary} "
            "-l {log}"
        )

rule concat_pairstat:
    input:
        pairstats=expand(
            f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.pairstat",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_pairstats_summary.tsv",
    log:
        "logs/concat_pairstat/concat_pairstat.log",
    run:
        from pathlib import Path

        # Build library names based on the order of input files (parent folder name)
        library_names = [Path(f).parent.name for f in input.pairstats]
        input_files_str = " ".join(input.pairstats)
        names_str = " ".join(library_names)

        shell(
            "python workflow/scripts/concat_stat.py "
            "-i {input_files_str} "
            "-n {names_str} "
            "-o {output.summary} "
            "-l {log}"
        )


rule concat_RSstat:
    input:
        RSstats=expand(
            f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.RSstat",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_RSstats_summary.tsv",
    log:
        "logs/concat_RSstat/concat_RSstat.log",
    run:
        from pathlib import Path

        # Build library names based on the order of input files (parent folder name)
        library_names = [Path(f).parent.name for f in input.RSstats]
        input_files_str = " ".join(input.RSstats)
        names_str = " ".join(library_names)

        shell(
            "python workflow/scripts/concat_stat.py "
            "-i {input_files_str} "
            "-n {names_str} "
            "-o {output.summary} "
            "-l {log}"
        )


rule concat_hicpro_pairstat:
    input:
        hicproPairStats=expand(
            f"{pairs_library_folder}/{{library}}.allValidPairs.mergestat",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_hicpro_pairstats_summary.tsv",
    log:
        "logs/concat_RSstat/concat_hicpro_pairstats.log",
    run:
        from pathlib import Path

        # Build library names based on the order of input files (parent folder name)
        library_names = [Path(f).parent.name for f in input.hicproPairStats]
        input_files_str = " ".join(input.hicproPairStats)
        names_str = " ".join(library_names)

        shell(
            "python workflow/scripts/concat_stat.py "
            "-i {input_files_str} "
            "-n {names_str} "
            "-o {output.summary} "
            "-l {log}"
        )



def get_all_pairstats_for_library(wc):
    # Loop over all runs for this library and gather pairstats
    all_files = []
    for run in LIBRARY_RUN_FASTQS[wc.library].keys():
        # Trigger the checkpoint for this run
        checkpoints.chunk_fastq.get(library=wc.library, run=run)
        # Discover chunk_ids from the processed fastq files (which exist after chunking)
        # Use side 2 fastq to determine chunk_ids
        chunk_ids = glob_wildcards(
            f"{processed_fastqs_folder}/{wc.library}/{run}/2.{{chunk_id}}.fastq.gz"
        ).chunk_id
        # Normalize chunk_ids: strip "_trimmed" suffix if present
        chunk_ids = [cid.replace("_trimmed", "") for cid in chunk_ids]
        # Add all pairstat files for this run
        all_files.extend(
            expand(
                f"{mapped_parsed_sorted_chunks_folder}/{wc.library}/{run}/{{chunk_id}}.pairstat",
                chunk_id=chunk_ids,
            )
        )
    return sorted(all_files)

def get_all_mapstats_for_library(wc):
    # Loop over all runs for this library and gather mapstats
    all_files = []
    for run in LIBRARY_RUN_FASTQS[wc.library].keys():
        # Trigger the checkpoint for this run
        checkpoints.chunk_fastq.get(library=wc.library, run=run)
        # Discover chunk_ids from the processed fastq files (which exist after chunking)
        # Use side 2 fastq to determine chunk_ids
        chunk_ids = glob_wildcards(
            f"{processed_fastqs_folder}/{wc.library}/{run}/2.{{chunk_id}}.fastq.gz"
        ).chunk_id
        # Normalize chunk_ids: strip "_trimmed" suffix if present
        chunk_ids = [cid.replace("_trimmed", "") for cid in chunk_ids]
        # Add both sides for all chunks in this run
        all_files.extend(
            expand(
                f"{mapped_parsed_sorted_chunks_folder}/{wc.library}/{run}/{{chunk_id}}_{{side}}.mapstat",
                chunk_id=chunk_ids,
                side=["1","2"],
            )
        )
    return sorted(all_files)


def get_all_RSstats_for_library(wc):
    # Loop over all runs for this library and gather mapstats
    all_files = []
    for run in LIBRARY_RUN_FASTQS[wc.library].keys():
        # Trigger the checkpoint for this run
        checkpoints.chunk_fastq.get(library=wc.library, run=run)
        # Discover chunk_ids from the processed fastq files (which exist after chunking)
        # Use side 2 fastq to determine chunk_ids
        chunk_ids = glob_wildcards(
            f"{processed_fastqs_folder}/{wc.library}/{run}/2.{{chunk_id}}.fastq.gz"
        ).chunk_id
        # Normalize chunk_ids: strip "_trimmed" suffix if present
        chunk_ids = [cid.replace("_trimmed", "") for cid in chunk_ids]
        # Add both sides for all chunks in this run
        all_files.extend(
            expand(
                f"{mapped_parsed_sorted_chunks_folder}/{wc.library}/{run}/{{chunk_id}}.hicpro.RSstat",
                chunk_id=chunk_ids,
            )
        )
    return sorted(all_files)


rule combine_map_pair_stats:
    input:
        get_all_pairstats_for_library
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.pairstat"
    log:
        "logs/combine_map_pair_stats/{library}.log"
    shell:
        "python workflow/scripts/merge_statfiles.py -f {input} > {output} 2>{log}"

rule combine_map_mapstat:
    input:
        get_all_mapstats_for_library
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.mapstat"
    log:
        "logs/combine_map_mapstat/{library}.log"
    shell:
        "python workflow/scripts/merge_statfiles.py -f {input}> {output} 2>{log}"

rule combine_RSstat:
    input:
        get_all_RSstats_for_library
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.RSstat"
    log:
        "logs/combine_RSstats/{library}.log"
    shell:
        "python workflow/scripts/merge_statfiles.py -f {input}> {output} 2>{log}"


rule multiqc:
    input:
        fastqc,
        expand(
            f"{pairs_library_folder}/{{library}}.{assembly}.dedup.stats",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
        expand(
            f"{stats_library_group_folder}/{{library_group}}.{assembly}.stats",
            library_group=config["input"]["library_groups"].keys(),
        )
        if "library_groups" in config["input"] and len(config["input"]["library_groups"]) > 0
        else [],
    conda:
        "../envs/multiqc.yml"
    log:
        "logs/multiqc.log",
    params:
        input_dirs=lambda wildcards, input: list(set([Path(f).parent for f in input])),
        outdir=lambda wildcards, output: Path(output[0]).parent,
    output:
        report=f"{multiqc_folder}/multiqc_report.html",
        dir=directory(multiqc_folder),
    shell:
        r"""multiqc -f --outdir {output.dir} -m pairtools -m fastqc \
        {params.input_dirs} \
        >{log} 2>&1"""


