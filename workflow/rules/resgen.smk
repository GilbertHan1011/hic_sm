rule resgen_upload_library_group:
    input:
        mcool=f"{coolers_library_group_folder}/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    conda:
        "../envs/resgen_python.yml"
    log:
        f"logs/resgen_upload_library_group/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.log",
    params:
        user=config["resgen"]["user"],
        project=config["resgen"]["project"],
        assembly=assembly,
    output:
        touch(
            f"{coolers_library_group_folder}/{{library_group}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool_uploaded"
        ),
    shell:
        r"""resgen sync datasets '{params.user}' '{params.project}' {input.mcool} \
        --tag filetype:cooler --tag assembly:{params.assembly} --tag datatype:matrix \
        --force-update true \
        >{log[0]} 2>&1
        """


use rule resgen_upload_library_group as resgen_upload_library with:
    input:
        mcool=f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool",
    log:
        f"logs/resgen_upload_library/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.log",
    output:
        touch(
            f"{coolers_library_folder}/{{library}}.{assembly}.{{filter_name}}.{{min_resolution}}.mcool_uploaded"
        ),