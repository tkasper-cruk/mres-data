from os import path
rule copyvae:
    input:
        h5ad_file = "results/cellranger_post/{sample_id}/filtered.h5ad"
    output:
        directory("results/copyvae/{sample_id}")
    singularity:
        path.join(workflow.basedir,"images/copyvae.sif")
    shell:
        """
        mkdir -p {output}
        copyvae -o {output} -cs {config[copyvae_cluster_size]} -eps {config[copyvae_epsilon]} {input.h5ad_file} {config[copyvae_cycle_genes]}"""