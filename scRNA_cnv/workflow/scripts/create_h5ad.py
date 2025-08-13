import scanpy as sc

features_path = snakemake.input["cr_dir"] + "/outs/filtered_feature_bc_matrix"

adata = sc.read_10x_mtx(
    path = features_path,
    make_unique = True
)

adata.write_h5ad(snakemake.output[0])