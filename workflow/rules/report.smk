

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


