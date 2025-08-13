from os import path
rule copyKAT:
    input:
        cr_dir="results/cellranger_post/{sample_id}"
    params:
        s_id="{sample_id}",
    output:
        directory("results/copykat/{sample_id}")
    conda:
        path.join(workflow.basedir,"envs/copykat.yaml")
    script:
        path.join(workflow.basedir,"scripts/run_copyKAT.R")