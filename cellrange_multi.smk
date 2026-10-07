rule cellrange_multi:
    input:
        dir_in + "/config_file/{sample}.csv"
    output:
        dir_out + "/{sample}/cellranger/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz"
    log:
        dir_out+"/{sample}/log/{sample}.cellranger.log"
    benchmark:
        dir_out+"/{sample}/log/{sample}.cellranger.bmk"
    params:
        id="{sample}",
	output_dir=dir_out + "/{sample}/cellranger"
    shell:
        '''
	module load cellranger
        cellranger multi --output-dir {params.output_dir} \
        --id {params.id} --csv {input} \
        --jobmode=local --localcores 8 \
        --localmem 114 &>{log}
        '''
