rule dist_vs_contact:
    input:
        mcool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    output:
        dist_stat = f"{downstream_dist_folder}/{{library}}.{{filter_name}}.{{min_resolution}}_loglog_fits.csv",
        pic1 = f"{downstream_dist_folder}/{{library}}.{{filter_name}}.{{min_resolution}}_hexbin_all_arms.png",
        pic2 = f"{downstream_dist_folder}/{{library}}.{{filter_name}}.{{min_resolution}}_per_arm.png"
    log:
        "logs/downstream_dist/{library}.{filter_name}.{min_resolution}.log",
    benchmark:
        "benchmarks/downstream_dist/{library}.{filter_name}.{min_resolution}.tsv"
    conda:
        "../envs/downstream_plot.yml"
    params:
        genome = config["input"]["genome"]["assembly_name"],
        out_dir = downstream_dist_folder
    threads: 8
    shell:
        r"""
        python workflow/scripts/calculate_dist_contact.py -i {input.mcool} \
        -g {params.genome} -r {wildcards.min_resolution} -t {threads} \
        -o {params.out_dir} \
        -p {wildcards.library}.{wildcards.filter_name}.{wildcards.min_resolution} >{log[0]} 2>&1
        """

rule mustache_loop_detection:
    input:
        mcool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{min_resolution}.mcool",
    output:
        loops = f"{downstream_loops_folder}/{{library}}.{{filter_name}}.{{resolution}}_loops.tsv"
    log:
        "logs/downstream_loops/{library}.{filter_name}.{resolution}.log",
    benchmark:
        "benchmarks/downstream_loops/{library}.{filter_name}.{resolution}.tsv"
    conda:
        "../envs/mustache.yml"
    params:
        threshold = lambda wildcards: config["mustache"]["thresholds"].get(wildcards.resolution, "0.1"),
        sparse = lambda wildcards: config["mustache"]["sparse_params"].get(wildcards.resolution, "0.88"),
        mustache_path = config["mustache"]["executable_path"],
        prefix = "{library}.{filter_name}.{resolution}",
        out_dir = downstream_loops_folder
    threads: 20
    shell:
        r"""
        bash workflow/scripts/mustache_loop_detection.sh \
        -i {input.mcool} \
        -r {wildcards.resolution} \
        -o {params.out_dir} \
        -p {params.prefix} \
        -t {threads} \
        -pt {params.threshold} \
        -st {params.sparse} \
        -m {params.mustache_path} >{log[0]} 2>&1
        """