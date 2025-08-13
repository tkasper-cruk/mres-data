from os import path

rule numbat_preprocessing:
    input:
        bamfile="results/cellranger/{sample_id}/outs/possorted_genome_bam.bam",
        filtered_mtx="results/cellranger/{sample_id}/outs/filtered_feature_bc_matrix"
    output:
        directory("results/numbat_pileup/{sample_id}")
    params:
        s_id="{sample_id}",
        resdir=path.join(workflow.basedir,"../results"),
        image=path.join(workflow.basedir,"images/numbat.sif"),
        homedir=workflow.basedir
    shell:
        """
        
        apptainer exec \
        -B {params.resdir}/:/data/ \
        {params.image} \
        Rscript --vanilla workflow/scripts/pileup_and_phase.R {params.s_id} {threads}
        """
        
        
        #-H 
rule numbat:
    input:
        preproc_dir="results/numbat_pileup/{sample_id}",
        filtered_mtx="results/cellranger/{sample_id}/outs/filtered_feature_bc_matrix"
    output:
        directory("results/numbat/{sample_id}")
    params:
        s_id="{sample_id}"
    conda:
        path.join(workflow.basedir,"envs/numbat.yaml")
    script:
        path.join(workflow.basedir,"scripts/run_numbat.R")
