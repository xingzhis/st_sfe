import re
import scanpy as sc
import pandas as pd
import numpy as np
import json
import importlib.resources as pkg_resources
from stfe import reference

def find_genes(nameptrn, gene_list, start=True, end=False):
    if start:
        nameptrn = f"^{nameptrn}"
    if end:
        nameptrn = f"{nameptrn}$"
    pattern = re.compile(nameptrn, re.IGNORECASE)
    filtered_strings = [string for string in gene_list if pattern.search(string)]
    return filtered_strings

def get_markers(pattern_list, adata, max_markers=100, sort_by_dispersion=True, start=True, end=False):
    markers = []
    for g in pattern_list:
        markers += find_genes(g, adata.var_names, start=start, end=end)
    if sort_by_dispersion:
        markers = adata.var['dispersions_norm'][markers].sort_values(ascending=False).index
    markers = markers[:max_markers]
    return markers

def get_markers_from_var_names(pattern_list, var_names, start=True, end=False):
    markers = []
    for g in pattern_list:
        markers += find_genes(g, var_names, start=start, end=end)
    return markers

def get_proliferation_score_cell_cycle(adata):
    with pkg_resources.open_binary(reference, 'WOT_genesets.xlsx') as excel_file:
        proliferation_genes = pd.read_excel(excel_file)["Proliferation"].dropna()
    proliferation_genes = pd.Series([p.upper() for p in proliferation_genes])

    sc.tl.score_genes(adata, proliferation_genes, score_name="proliferation_score")

    with pkg_resources.open_text(reference, 'regev_lab_cell_cycle_genes.txt') as txt_file:
        cell_cycle_genes = [x.strip() for x in txt_file]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    scores_df = sc.get.obs_df(adata, ['phase', 'S_score', 'G2M_score', 'proliferation_score'])
    return scores_df

def get_emt_score(adata):
    EMT_genes = [
        'ZEB1',
        'LIX1L',
        'VIM',
        'AXL',
        'MMP2',
        'ANTXR2', 
        #'C3orf21', 
        'FN1',
        'NRP1',
        'TGFBI',
        'GALNT5',
        'PPARG',
        'HNMT',
        'CARD6', 
        'RBPMS',
        'TNFRSF21',
        'TMEM45B',
        'MPP7', 'SSH3', 
        #'MTAC2D1', 
        #'MUC', 
        'EPPK1', 'SHROOM3', 'EPN3', 'PRSS22', 'AP1M2', 
        'SH3YL1', 'KLC3', 'SERINC2', 
        #'EBPM', 
        'FXYD3', 'CLDN4', 
        'CRB3', 
        #'LRRC54', 
        'MAPK13', 'GALNT3', 'STAP2', 'AP1M2', 'DSP', 'ELMO3', 
        'KRTCAP3', 'MAL2', 'F11R', 
        #'GPR110', 
        #'GPR56', 
        'KRT19', 'GRHL1', 'BSPRY', 
        'C1orf116', 'S100A14', 'SPINT2', 'ANKRD22', 'ST14', 'GRHL2', 'PRR5', 
        'BSPRY', 'TJP3', 'TACSTD2', 'CDH3', 
        #'C1orf172', 
        'CDS1', 
        'MPZL2', 
        #'INADL', 
        'EPN3', 
        #'RBM35A', 
        'TMC4', 'ITGB6', 'TMEM125', 
        'EPHA1', 'CDS1', 
        'ENPP5', 'ST14', 'EPB41L5', 'ERBB3', 'RAB25', 
        'PRSS8', 'TMEM30B', 'CLDN7', 
        #'RBM35A', 
        #'TACSTD1', 
        'CDS1', 'SCNN1A',
        'CDH1',
    ]
    emt_expr = sc.get.obs_df(adata, keys=EMT_genes)
    ecad_expr = sc.get.obs_df(adata, keys=['CDH1'])
    weights = emt_expr.corr()[['CDH1']]
    emt_scores = pd.DataFrame((np.array(emt_expr) * np.array(weights).T), index=emt_expr.index, columns=emt_expr.columns).sum(axis=1)
    return emt_scores

def get_genes_from_json(json_file, adata, max_markers=1000, start=True, end=False):
    with open(json_file, 'r') as file:
        data = json.load(file)
    
    gene_list_found_dict = {}
    for gene_dict in data:
        preocess = gene_dict['process']
        gene_list = [g['gene'] for g in gene_dict['gates']]
        gene_list_found_dict[preocess] = get_markers(gene_list, adata, max_markers=max_markers, start=start, end=end)
    
    return gene_list_found_dict

def get_ligand_receptor_lists(genelist):
    with pkg_resources.open_text(reference, 'human_lr_pair.txt') as lr_file:
        lrlist = pd.read_csv(lr_file, sep='\t')
    ligand_list = []
    receptor_list = []
    for i, row in lrlist.iterrows():
        ligands = find_genes(row['ligand_gene_symbol'], genelist, True, True)
        receptors = find_genes(row['receptor_gene_symbol'], genelist, True, True)
        if len(ligands) > 0 and len(receptors) > 0:
            assert len(ligands) == 1 and len(receptors) == 1, (len(ligands), len(receptors), ligands, receptors, '\n', row)
            ligand_list.extend(ligands)
            receptor_list.extend(receptors)
    return ligand_list, receptor_list

def get_ligand_receptor_pairs(genelist):
    ligand_list, receptor_list = get_ligand_receptor_lists(genelist)
    assert len(ligand_list) == len(receptor_list), (len(ligand_list), len(receptor_list))
    n_pairs = len(ligand_list)
    lr_pairs = ligand_list + receptor_list
    return lr_pairs, n_pairs


