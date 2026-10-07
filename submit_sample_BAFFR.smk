dir_out="/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell"
intake_file = "/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell/id/intake_file.tsv"

import pandas as pd
#sample list
intake=pd.read_table(intake_file)
#sample id and bam id
samplelist=list(intake.sample_id.drop_duplicates())


rule all:
    input:
        expand(dir_out+"/smk_code/{sample}.sh",sample=samplelist)

rule write_code:
    input:
        key_smk="BAFFR_singlecell.smk"
    output:
        code=dir_out+"/smk_code/{sample}.sh"
    params:
        id="{sample}",
        cores=8,
        mem="128G",
        log=dir_out+"/smk_log/{sample}_SeuratQC_%j.log",
        key_file= dir_out+ "/{sample}/log/SeuratQC.bmk",
        dir = dir_out+"/{sample}"
    shell:
        '''
        echo "#!/bin/bash" > {output.code}
        echo "#SBATCH --job-name=seuratqc_{wildcards.sample}           # Job name" >> {output.code}
        echo "#SBATCH -n {params.cores}                          # Number of cores" >> {output.code}
        echo "#SBATCH -N 1-1                        # Min - Max Nodes" >> {output.code}
        echo "#SBATCH -p all                       # gpu queue" >> {output.code}
        echo "#SBATCH --mem={params.mem}                     # Amount of memory in GB" >> {output.code}
        echo "#SBATCH --time=24:00:00               # Time limit hrs:min:sec" >> {output.code}
        echo "#SBATCH --output={params.log}    # Standard output and error log" >> {output.code}
        echo "snakemake -s {input.key_smk} -p {params.key_file} -j{params.cores}"  >> {output.code}
        '''
