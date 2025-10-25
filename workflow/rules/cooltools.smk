

rule bin_pairs_library:
    input:
        pairs=f"{pairs_library_folder}/{{library}}.{assembly}.nodups.pairs.gz",
        chromsizes=chromsizes_path,
    params:
        filter_command=lambda wildcards: (
            f'| pairtools select "{config["bin"]["filters"][wildcards.filter_name]}"'
            if config["bin"]["filters"][wildcards.filter_name]
            else ""
        ),
        assembly=assembly,
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        cool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.cool",
    log:
        "logs/bin_pairs_library/{library}.{filter_name}.{min_resolution}.log",
    benchmark:
        "benchmarks/bin_pairs_library/{library}.{filter_name}.{min_resolution}.tsv"
    shell:
        r"""
        bgzip -dc -@ {threads} {input.pairs} {params.filter_command} | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
        --assembly {params.assembly} \
        {input.chromsizes}:{wildcards.min_resolution} - {output.cool} >{log[0]} 2>&1
        """


rule zoom_library:
    input:
        cool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.cool",
    params:
        res_string=",".join([str(res) for res in config["bin"]["resolutions"]]),
        balance_args=lambda wildcards, threads: (
            f"--balance --balance-args '{config['bin'].get('balance_options', '')} --nproc {threads}'"
            if config["bin"]["balance"]
            else ""
        ),
    threads: 8
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        mcool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    log:
        "logs/zoom_library/{library}.{filter_name}.{min_resolution}.log",
    benchmark:
        "benchmarks/zoom_library/{library}.{filter_name}.{min_resolution}.tsv"
    shell:
        r"""
        cooler zoomify \
        --nproc {threads} \
        --out {output.mcool} \
        --resolutions {params.res_string} \
        {params.balance_args} \
        {input.cool} \
        >{log[0]} 2>&1
        """


rule merge_zoom_library_group_coolers:
    input:
        lambda wildcards: expand(
            f"{coolers_library_folder}/{{library}}.{assembly}.{wildcards.filter_name}.{wildcards.min_resolution}.cool",
            library=config["input"]["library_groups"][wildcards.library_group],
        ),
    params:
        res_string=",".join([str(res) for res in config["bin"]["resolutions"]]),
        balance_args=lambda wildcards, threads: (
            f"--balance --balance-args '{config['bin'].get('balance_options', '')} --nproc {threads}'"
            if config["bin"]["balance"]
            else ""
        ),
    threads: 8
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        cool=f"{coolers_library_group_folder}/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.cool",
        mcool=f"{coolers_library_group_folder}/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    log:
        "logs/merge_zoom_library_group_coolers/{library_group}.{filter_name}.{min_resolution}.log",
    benchmark:
        "benchmarks/merge_zoom_library_group_coolers/{library_group}.{filter_name}.{min_resolution}.tsv"
    shell:
        r"""
        cooler merge {output.cool} {input} >{log[0]} 2>&1

        cooler zoomify \
        --nproc {threads} \
        --out {output.mcool} \
        --resolutions {params.res_string} \
        {params.balance_args} \
        {output.cool} \
        >{log[0]} 2>&1
        """


rule mcool2hic_group:
    input:
        mcool=f"{coolers_library_group_folder}/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    output:
        hic=f"{coolers_library_group_folder}/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.hic",
    log:
        "logs/mcool2hic/{library_group}.{filter_name}.{min_resolution}.log",
    benchmark:
        "benchmarks/mcool2hic/{library_group}.{filter_name}.{min_resolution}.tsv"
    conda:
        "../envs/hictk.yml"
    threads: 8
    params:
        tmpdir=coolers_library_group_folder,
    shell:
        r"""
        hictk convert --threads {threads} --tmpdir {params.tmpdir} \
        {input.mcool} {output.hic} >{log} 2>&1
        """


use rule mcool2hic_group as mcool2hic_library with:
    input:
        mcool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    output:
        hic=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.hic",
    log:
        "logs/mcool2hic/{library}.{filter_name}.{min_resolution}.log",
    benchmark:
        r"benchmarks/mcool2hic/{library}.{filter_name}.{min_resolution}.tsv"


rule scaling_pairs_library:
    input:
        pairs=f"{pairs_library_folder}/{{library}}.{assembly}.nodups.pairs.gz",
        chromsizes=chromsizes_path,
    params:
        min_distance=config["scaling_pairs"]["min_distance"],
        max_distance=config["scaling_pairs"]["max_distance"],
        n_dist_bins_decade=config["scaling_pairs"]["n_dist_bins_decade"],
        extra=config["scaling_pairs"]["scaling_options"],
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        scaling=f"{pairs_library_folder}/{{library}}.{assembly}.nodups.scaling.tsv",
    log:
        "logs/scaling_pairs_library/{library}.log",
    benchmark:
        "benchmarks/scaling_pairs_library/{library}.tsv"
    shell:
        r"""
        pairtools scaling \
        --dist-range {params.min_distance} {params.max_distance} \
        --n-dist-bins-decade {params.n_dist_bins_decade} \
        {params.extra} \
        -o {output.scaling} \
        {input.pairs} \
        >{log[0]} 2>&1
        """