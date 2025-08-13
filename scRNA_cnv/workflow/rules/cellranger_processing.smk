from os import path

rule rename_fastq:
    input:
        raw_fq_dir="results/raw_reads/{sample_id}"
    output:
        temp(directory("results/fastqs_renamed/{sample_id}"))
    params:
        s_id="{sample_id}"
    script:
        path.join(workflow.basedir,"scripts/rename_fastqs.py")

rule cellranger:
    input:
        fqdir="results/fastqs_renamed/{sample_id}"
    params:
        "--id {sample_id}",
        "--create-bam true",
        "--sample {sample_id}",
    output:
        cr_dir=directory("results/{sample_id}/cellranger/"),
    shell:
        "{config[cellranger_bin]} count  --transcriptome {config[cellranger_ref]} --fastqs {input} --localmem {config[cellranger_mem]} --localcores {threads} --output-dir {output[cr_dir]} --chemistry {config[cellranger_chemistry]} {params}"  

rule move_cellranger:
    input:
        cr_dir="results/{sample_id}/cellranger/"
    output:
        crdir=directory("results/cellranger/{sample_id}")
    shell:
        "mv {input[0]}* {output}"

rule copy_notebook:
    input: 
        notebook=path.join(workflow.basedir,"notebooks/cellranger_post.py.ipynb"),
        cbdir="results/cellbender/{sample_id}"
    output:
        temp("results/{sample_id}_notebook_copied")
    params:
        s_id="{sample_id}"
    shell:
        """
        cp {input[0]} {input[1]}/postprocessing_{params[s_id]}.ipynb
        touch {output}
        """

rule cellbender:
    input:
        crdir="results/cellranger/{sample_id}"
    output:
        cbdir=directory("results/cellbender/{sample_id}")
    params:
        checkpoint= path.join(workflow.basedir,"../ckpt.tar.gz")
    conda:
        path.join(workflow.basedir,"envs/cellbender.yaml") #only for cuda/cpu
    shell:
        """
        mkdir -p {output[0]}/..
        python -c 'import torch; print(torch.cuda.is_available())'
        cellbender remove-background --input {input[0]}/outs/raw_feature_bc_matrix.h5 --output {output[0]}/feature_bc_matrix.h5 --cuda --cpu-threads {threads} && rm -f {params}
        """
        #