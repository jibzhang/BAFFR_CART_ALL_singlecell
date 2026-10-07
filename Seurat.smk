rule createSeuratobj:
    input: dir_out + "/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz"
    output: dir_out + "/{sample}/Seurat/{sample}.SeuratObj.raw.rds"
    log: dir_out + "/{sample}/log/createSeuratobj.log"
    benchmark: dir_out + "/{sample}/log/createSeuratobj.bmk"
    params:
        dir_cellranger_matrix=dir_out + "/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix",
        id="{sample}",
        outdir=dir_out + "/{sample}/Seurat",
        adt_feature="adt_EGFR",
        CAR_feature="WPRE601",
        transform = "basic"
    shell:
        '''
        module load R/RStudio_R-4.4.1
        Rscript  /home/jibzhang/pipeline/Singlecell/Rscript/CreatSeuratObj_multimodal.R {params.dir_cellranger_matrix} {params.id} {params.outdir} {params.adt_feature} {params.CAR_feature} {params.transform} &> {log}
        '''

rule SeuratQC:
    input: rules.createSeuratobj.output
    output: dir_out + "/{sample}/log/SeuratQC.bmk"
    log: dir_out + "/{sample}/log/SeuratQC.log"
    benchmark: dir_out + "/{sample}/log/SeuratQC.bmk"
    params:
        rds_out="{sample}.SeuratObj.QC.rds",
        doublet_removed = "{sample}.doublet_removed.rds",
        transform = "basic"
    shell:
        '''
        module load R/RStudio_R-4.4.1
        Rscript  /home/jibzhang/pipeline/Singlecell/Rscript/SeuratQC.R {input} {params.doublet_removed} {params.rds_out} {params.transform} &> {log}
        '''