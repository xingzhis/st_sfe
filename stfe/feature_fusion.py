import numpy as np
import pandas as pd
import scipy.sparse as sp

EPSILON = 1e-5

def transform_spatial(spatial_dict):

    std_dict = {}
    cols_list = []
    preproc_list = []
    for key, value in spatial_dict.items():
        feat_df = value.copy()
        assert not feat_df.isna().any(), f"{key} contains nan!"
        if key.endswith('_percentages'):
            # feat_df = feat_df.drop(feat_df.columns[-1], axis=1)
            feat_df = feat_df.applymap(lambda x: np.sqrt(x + 1e-7)) # take sqrt instead
        stds = np.std(feat_df, axis=0)
        std_dict[key] = stds
        feat_df = (feat_df + EPSILON) / (stds + EPSILON)
        cols_list.append({'key':key, 'cols': feat_df.columns})
        preproc_list.append(feat_df.values)
    preproc_data = np.concatenate(preproc_list, axis=1)
    # preproc_std = np.std(preproc_data)
    # preproc_data = preproc_data / preproc_std
    preproc_info = dict(
        std_dict=std_dict,
        cols_list=cols_list,
        # preproc_std=preproc_std,
    )
    return preproc_data, preproc_info

def merge_gene_spatial(adata, spatial_dict):
    gene_expr = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    preproc_data, preproc_info = transform_spatial(spatial_dict)
    gene_std = np.std(gene_expr)
    preproc_std = np.std(preproc_data)
    preproc_data = (preproc_data + EPSILON) / (preproc_std + EPSILON) * gene_std
    final_data = np.c_[gene_expr, preproc_data]
    preproc_info['gene_std'] = gene_std
    preproc_info['n_genes'] = gene_expr.shape[1]
    preproc_info['preproc_std'] = preproc_std
    return final_data, preproc_info

# def merge_gene_spatial(adata, spatial_dict):
#     gene_expr = adata.X.toarray() if sp.issparse(adata.X) else adata.X
#     # spatial_dict = adata.uns['spatial_features_raw']
#     std_dict = {}
#     cols_list = []
#     preproc_list = []
#     for key, value in spatial_dict.items():
#         feat_df = value.copy()
#         if key.endswith('_percentages'):
#             # feat_df = feat_df.drop(feat_df.columns[-1], axis=1)
#             feat_df = feat_df.applymap(lambda x: np.sqrt(x + 1e-7)) # take sqrt instead
#         stds = np.std(feat_df, axis=0)
#         std_dict[key] = stds
#         feat_df = feat_df / stds
#         cols_list.append({'key':key, 'cols': feat_df.columns})
#         preproc_list.append(feat_df.values)
#     preproc_data = np.concatenate(preproc_list, axis=1)
#     preproc_std = np.std(preproc_data)
#     gene_std = np.std(gene_expr)
#     preproc_data = preproc_data / preproc_std * gene_std
#     final_data = np.c_[gene_expr, preproc_data]
#     preproc_info = dict(
#         std_dict=std_dict,
#         cols_list=cols_list,
#         preproc_std=preproc_std,
#         gene_std=gene_std,
#     )
#     return final_data, preproc_info

def spatial_inverse_transform(preproc_data, preproc_info):
    preproc_data = preproc_data * (preproc_info['preproc_std'] + EPSILON) / preproc_info['gene_std'] - EPSILON
    rec_data = {}
    col_counter = 0
    for info in preproc_info['cols_list']:
        # rec_data['name'] = info['key']
        len_cols = len(info['cols'])
        this_data = preproc_data[:,col_counter:col_counter+len_cols]
        col_counter += len_cols
        this_data = pd.DataFrame(this_data, columns=info['cols'])
        this_data = this_data * (preproc_info['std_dict'][info['key']] + EPSILON) - EPSILON
        if info['key'].endswith('_percentages'):
            this_data = this_data**2 - 1e-7
        rec_data[info['key']] = this_data
    return rec_data

def inverse_transform(final_data, preproc_info, n_genes):
    gene_data = final_data[:, :n_genes]
    preproc_data = final_data[:, n_genes:]
    rec_data = spatial_inverse_transform(preproc_data, preproc_info)
    return gene_data, rec_data
