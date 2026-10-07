dir_out="/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell"
intake_file = "/net/isi-dcnl/ifs/user_data/lkwak/Data/Jibin/project/BAFFR_singlecell/id/intake_file.tsv"

import pandas as pd
intake=pd.read_table(intake_file)
samplelist=list(intake.sample_id.drop_duplicates())

include: "Seurat.smk"

Seurat = expand(dir_out + "/{sample}/log/SeuratQC.bmk", sample=samplelist)

TASKS = []
TASKS.extend(Seurat)

rule all:
    input: TASKS