"""
GVHD Single Cell RNA-Seq Pipeline
"""
import os

configfile: "config/config.yml"

out_dir = os.path.join(config['output_dir'], config['version'])

rule all:
    input:
        os.path.join(out_dir, "reports", "data_overview.html")

rule data_overview:
    input:
        os.path.join(out_dir, "seurat", "seurat_obj.rds")
    output:
        os.path.join(out_dir, "reports", "data_overview.html")
    script:
        "rmd/data_overview.Rmd"

rule build_seurat_obj:
    output:
        os.path.join(out_dir, "seurat", "seurat_obj.rds")
    script:
        "scripts/build_seurat_obj.R"
