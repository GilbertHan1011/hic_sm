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