import os
import re
import sys
import csv
import json
import tqdm
import pyobo
import obonet
import random
import pickle
import logging
import graphviz
import datetime
import openpyxl
import networkx
import itertools
import numpy as np
import pandas as pd
import enzyme_client
from pathlib import Path
from matplotlib import rc
from bioinfokit import visuz
from graphviz import Digraph
from indra.sources import tas
import matplotlib.pyplot as plt
from indra.util import batch_iter
from collections import OrderedDict
from collections import defaultdict
import matplotlib.colors as mcolors
from indra.statements import Complex
from scipy.stats import fisher_exact
from indra.sources import indra_db_rest
import indra.tools.assemble_corpus as ac
from indra.literature import pubmed_client
from indra.assemblers.cx import hub_layout
from indra.ontology.bio import bio_ontology
from indra.databases.uniprot_client import um
from indra.assemblers.html import HtmlAssembler
from indra.statements.agent import default_ns_order
from indra.sources.omnipath import process_from_web
from indra.assemblers.cx.assembler import CxAssembler
from indra.databases import uniprot_client, hgnc_client
from indra_db.client.principal.curation import get_curations
from indra.databases.hgnc_client import get_hgnc_from_mouse, get_hgnc_name

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, os.pardir, 'input')
OUTPUT = os.path.join(HERE, os.pardir, 'output')

GO_ANNOTATIONS = os.path.join(INPUT, 'goa_human.gaf')
INDRA_DB_PKL = os.path.join(INPUT, 'db_dump_df.pkl')
DATA_SPREADSHEET = os.path.join(INPUT, 'Neuroimmune gene list .xlsx')
DRUG_BANK_PKL = os.path.join(INPUT, 'drugbank_5.1.pkl')
ION_CHANNELS = os.path.join(INPUT, 'ion_channels.txt')
SURFACE_PROTEINS_WB = os.path.join(INPUT, 'Surface Proteins.xlsx')
HUMAN_PAIN_DB = os.path.join(INPUT, 'Human_Pain_Genes_DB.tsv')
IMMUNE_CELLTYPE_LIST = ['DCs',
                        'Dermal Macs',
                        'M2a',
                        'M2b',
                        'Monocytes',
                        'Resident Mac',
                        'Mast cells']

logger = logging.getLogger('receptor_ligand_interactions')

mouse_gene_name_to_mgi = {v: um.uniprot_mgi.get(k)
                          for k, v in um.uniprot_gene_name.items()
                          if k in um.uniprot_mgi}


def _load_goa_gaf():
    """Load the gene/GO annotations as a pandas data frame."""
    # goa_ec = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP',
    #          'HGI', 'HEP', 'IBA', 'IBD'}
    goa = pd.read_csv(GO_ANNOTATIONS, sep='\t',
                      skiprows=23, dtype=str,
                      header=None,
                      names=['DB',
                             'DB_ID',
                             'DB_Symbol',
                             'Qualifier',
                             'GO_ID',
                             'DB_Reference',
                             'Evidence_Code',
                             'With_From',
                             'Aspect',
                             'DB_Object_Name',
                             'DB_Object_Synonym',
                             'DB_Object_Type',
                             'Taxon',
                             'Date',
                             'Assigned',
                             'Annotation_Extension',
                             'Gene_Product_Form_ID'])
    goa = goa.sort_values(by=['DB_ID', 'GO_ID'])
    # Filter out all "NOT" negative evidences
    goa['Qualifier'].fillna('', inplace=True)
    goa = goa[~goa['Qualifier'].str.startswith('NOT')]
    # Filter to rows with evidence code corresponding to experimental
    # evidence
    # goa = goa[goa['Evidence_Code'].isin(goa_ec)]
    return goa


goa = _load_goa_gaf()


def get_pain_mol():
    PAIN_SIGNAL_MOL = {
        "Prostaglandins": "CHEBI:26333",
        "Brandykinin": "CHEBI:3165"
    }

    CHEBI_LIST = {}
    CHEBI_NAMES = {}
    for compounds, chebi_id in PAIN_SIGNAL_MOL.items():
        CHEBI_LIST[compounds] = \
            [children[1] for children in
             bio_ontology.get_children('CHEBI',
                                       chebi_id)]

        CHEBI_NAMES[compounds] = \
            [bio_ontology.get_name('CHEBI', ids)
             for ids in CHEBI_LIST[compounds]]

    return CHEBI_NAMES


PAIN_MOL_NAMES = get_pain_mol()


def load_indra_df(fname):
    """Return an INDRA Statement data frame from a pickle file."""
    logger.info('Loading INDRA DB dataframe')
    with open(fname, 'rb') as fh:
        df = pickle.load(fh)
    logger.info('Loaded %d rows from %s' % (len(df), fname))
    return df


def get_hashes_by_gene_pair(df, ligand_genes, receptor_genes):
    hashes_by_gene_pair = defaultdict(set)
    l_genes = list(ligand_genes.values())
    l_logFC = list(ligand_genes.keys())

    for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
        if a in l_genes and b in receptor_genes:
            logFC = l_logFC[l_genes.index(a)]
            hashes_by_gene_pair[(a, b, logFC)].add(hs)
    return hashes_by_gene_pair


def download_statements(hashes):
    """Download the INDRA Statements corresponding to a set of hashes.
    """
    stmts_by_hash = {}
    for group in tqdm.tqdm(batch_iter(hashes, 200), total=int(len(hashes) / 200)):
        idbp = indra_db_rest.get_statements_by_hash(list(group),
                                                    ev_limit=10)
        for stmt in idbp.statements:
            stmts_by_hash[stmt.get_hash()] = stmt
    return stmts_by_hash


def get_genes_for_go_ids(go_ids):
    """Return genes that are annotated with a given go ID or its children."""
    all_go_ids = set()
    for go_id in go_ids:
        children_go_ids = {ch[1] for ch in bio_ontology.get_children('GO', go_id)}
        all_go_ids.add(go_id)
        all_go_ids |= children_go_ids
    df = goa[goa['GO_ID'].isin(all_go_ids)]
    up_ids = sorted(list(set(df['DB_ID'])))
    gene_names = [uniprot_client.get_gene_name(up_id) for up_id in up_ids]
    gene_names = {g for g in gene_names if g}
    return gene_names


def fix_dates(gene_names):
    replacements = {
        datetime.datetime(2020, 3, 7, 0, 0): 'March7',
        datetime.datetime(2020, 3, 2, 0, 0): 'March2',
        datetime.datetime(2020, 3, 4, 0, 0): 'March4',
        datetime.datetime(2020, 3, 5, 0, 0): 'March5',
        datetime.datetime(2020, 3, 6, 0, 0): 'March6',
        datetime.datetime(2020, 3, 9, 0, 0): 'March9',
        datetime.datetime(2020, 3, 8, 0, 0): 'March8',
        datetime.datetime(2020, 3, 11, 0, 0): 'Mar11',
        datetime.datetime(2020, 9, 1, 0, 0): 'Sept1',
        datetime.datetime(2020, 9, 2, 0, 0): 'Sept2',
        datetime.datetime(2020, 9, 3, 0, 0): 'Sept3',
        datetime.datetime(2020, 9, 4, 0, 0): 'Sept4',
        datetime.datetime(2020, 9, 5, 0, 0): 'Sept5',
        datetime.datetime(2020, 9, 6, 0, 0): 'Sept6',
        datetime.datetime(2020, 9, 7, 0, 0): 'Sept7',
        datetime.datetime(2020, 9, 8, 0, 0): 'Sept8',
        datetime.datetime(2020, 9, 9, 0, 0): 'Sept9',
        datetime.datetime(2020, 9, 10, 0, 0): 'Sept10',
        datetime.datetime(2020, 9, 11, 0, 0): 'Sept11',
        datetime.datetime(2020, 9, 15, 0, 0): 'Sept15',
    }
    fixed_gene_names = set()
    for gene_name in gene_names:
        if isinstance(gene_name, datetime.datetime):
            fixed_gene_names.add(replacements[gene_name])
        else:
            fixed_gene_names.add(gene_name)
    return fixed_gene_names


def read_workbook(workbook):
    """ This function takes Excel workbook as an input and
    returns ligand and receptor gene list respectively.
    Input: Excel workbook with single(2 columns) or two sheets
    Condition: considers first column/sheet as ligand genes and second
    column/shet as receptor genes
    """
    ligands_sheet = 'updated list of ligands '
    receptors_sheet = 'RPKM > 1.5 cfiber'
    wb = openpyxl.load_workbook(workbook)
    ligands = fix_dates(set([row[0].value for row in wb[ligands_sheet]][1:]))
    receptors = fix_dates(set([row[0].value
                               for row in wb[receptors_sheet]][1:]))
    return ligands, receptors


def _plot_de_genes(df):
    os.chdir(output_dir)
    visuz.gene_exp.volcano(df=df,
                           lfc='avg_logFC', pv='p_val',
                           plotlegend=True, legendpos='upper right',
                           legendanchor=(1.46, 1), geneid="Genes",
                           genenames="deg", gstyle=2)


def read_gene_list(infile, mode):
    gene_list = []
    try:
        with open(infile, mode) as FH:
            for eachGene in FH:
                gene_list.append(eachGene.strip("\n"))
        return gene_list

    except FileNotFoundError:
        sys.exit("Given file doesn't exist")


def filter_nuclear_receptors(receptors_go, go_term):
    # Filtering out the nuclear receptors from the receptor list
    nuclear_receptors = get_genes_for_go_ids([go_term])
    # Add any others that don't have the right annotation
    nuclear_receptors |= {'NR2C2'}
    filtered_receptors_go = receptors_go - nuclear_receptors
    return filtered_receptors_go


def filter_complex_statements(stmts, ligands, receptors):
    for stmt in stmts:
        if isinstance(stmt, Complex):
            # Statement updated by reference here
            _filter_complex(stmt, ligands, receptors)
    return stmts


def _filter_complex(stmt, lg, rg):
    """Filter out the genes from Complex statements which
    are not present in the given ligand/receptor list"""
    stmt.members = [agent for agent in stmt.members
                    if agent.name in lg or agent.name in rg]
    return stmt


def filter_op_stmts(op_stmts, lg, rg):
    """ Filter out the statements which are not ligand and receptor """
    logger.info(f'Filtering {len(op_stmts)} to ligand-receptor interactions')
    filtered_stmts = [stmt for stmt in op_stmts if
                      (any(a.name in lg for a in stmt.agent_list())
                       and any(a.name in rg for a in stmt.agent_list()))]
    logger.info(f'{len(filtered_stmts)} left after filter')
    return filtered_stmts


def html_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a HTML report"""
    html_assembler = HtmlAssembler(indra_stmts,
                                   db_rest_url='https://db.indra.bio')
    assembled_html_report = html_assembler.make_model(no_redundancy=True)
    html_assembler.save_model(fname)
    return assembled_html_report


def cx_assembler(indra_stmts, fname):
    """Assemble INDRA statements into a CX report"""
    cx_assembler = CxAssembler(indra_stmts)
    assembled_cx_report = cx_assembler.make_model()
    cx_assembler.save_model(fname)
    ndex_network_id = cx_assembler.upload_model(ndex_cred=None,
                                                private=True, style='default')
    return assembled_cx_report, ndex_network_id


def get_small_mol_report(targets_by_drug, potential_targets, fname):
    df = []
    for drug, targets in targets_by_drug.items():
        targets_in_data = targets & potential_targets
        if not targets_in_data:
            continue
        df.append(
            {
                "Drug": drug[0],
                "ID": '%s:%s' % (drug[1]),
                "Named": 0 if drug[0].startswith('CHEMBL') else 1,
                "Score": "{:.3f}".format(len(targets_in_data) / len(targets)),
                "Number of targets in data": len(targets_in_data),
                "Targets in data": ", ".join(sorted(targets_in_data)),
                "Other targets": ", ".join(sorted(targets - targets_in_data)),
            }
        )
    df = pd.DataFrame(df).sort_values(by=['Score', 'Number of targets in data',
                                          'Named'],
                                      ascending=False)
    df.to_csv(fname, sep="\t", header=True, index=False)
    return df


"""
def get_ligands_by_receptor(receptors_in_data, ligands_in_data, stmts):
    ligands_by_receptor = defaultdict(set)
    logFC = list(ligands_in_data.keys())
    lg = list(ligands_in_data.values())

    for stmt in stmts:
        agent_names = {agent.name for agent in stmt.agent_list()}
        receptors = agent_names & receptors_in_data
        ligands = agent_names & ligands_in_data
        for receptor in receptors:
            ligands_by_receptor[receptor] |= ligands
    return dict(ligands_by_receptor)
"""


def get_ligands_by_receptor(receptors_in_data, ligands_in_data, stmts):
    ligands_by_receptor = defaultdict(set)
    for stmt in stmts:
        agent_names = {agent.name for agent in stmt.agent_list()}
        receptors = agent_names & receptors_in_data
        ligands = agent_names & ligands_in_data
        for receptor in receptors:
            ligands_by_receptor[receptor] |= ligands
    return dict(ligands_by_receptor)


def filter_out_medscan(stmts):
    logger.info('Filtering out medscan evidence on %d statements' % len(stmts))
    new_stmts = []
    for stmt in stmts:
        new_evidence = [e for e in stmt.evidence if e.source_api != 'medscan']
        if not new_evidence:
            continue
        stmt.evidence = new_evidence
        if not stmt.evidence:
            continue
        new_stmts.append(stmt)
    logger.info('%d statements after filter' % len(new_stmts))
    return new_stmts


def filter_db_only(stmts):
    new_stmts = []
    for stmt in stmts:
        sources = {ev.source_api for ev in stmt.evidence}
        if sources <= {'reach', 'sparser', 'trips', 'rlimsp', 'medscan', 'eidos'}:
            continue
        new_stmts.append(stmt)
    return new_stmts


def get_cell_type_stats(stmts, ligands, receptors):
    interactome = set()
    ligand_interactions = defaultdict(set)
    for stmt in stmts:
        stmt_ligands = {a.name for a in stmt.agent_list() if
                        a.name in ligands}
        stmt_receptors = {a.name for a in stmt.agent_list() if
                          a.name in receptors}
        for ligand, receptor in itertools.product(stmt_ligands,
                                                  stmt_receptors):
            interactome.add((ligand, receptor))
            ligand_interactions[ligand].add(receptor)
    return len(interactome), ligand_interactions


def plot_interaction_potential(num_interactions_by_cell_type, fname):
    labels = {
        'DCs': 'Dendritic cells',
        'Dermal Macs': 'Dermal macrophages',
        'M2a': 'Reparative macrophages (2a)',
        'M2b': 'Reparative macrophages (2b)',
        'Monocytes': 'Monocytes',
        'Resident Mac': 'Resident macrophages',
        'Mast cells': 'Mast cells'
    }
    G = networkx.DiGraph()
    for cell_type, num_int in num_interactions_by_cell_type.items():
        G.add_node(cell_type, label=labels[cell_type])
        G.add_edge(cell_type, 'Neurons', label=num_int)
    ag = networkx.nx_agraph.to_agraph(G)
    ag.draw(fname, prog='dot')


def get_all_enzymes():
    HOME = str(Path.home())
    ec_code_path = '.obo/ec-code/ec-code.obo'
    if not os.path.exists(os.path.join(HOME, ec_code_path)):
        _ = pyobo.get_id_name_mapping('ec-code')
        obo = obonet.read_obo(os.path.join(HOME, ec_code_path))
    else:
        obo = obonet.read_obo(os.path.join(HOME, ec_code_path))
    up_nodes = set()
    for node in obo.nodes:
        if node.startswith('uniprot'):
            up_nodes.add(node[8:])
    human_ups = {u for u in up_nodes if uniprot_client.is_human(u)}
    enzymes = {uniprot_client.get_gene_name(u) for u in human_ups}
    enzymes = {g for g in enzymes if not hgnc_client.is_kinase(g)}
    enzymes = {g for g in enzymes if not hgnc_client.is_phosphatase(g)}
    logger.info(f'Filtered {len(enzymes)} enzymes in total')
    return enzymes


def process_seurat_csv(infile, fc):
    """ Process Seurat dataframe and only filter in
    genes with the given Fold change """
    l_df = pd.read_csv(infile, header=0, sep=",")
    l_df.columns = l_df.columns.str.replace('Unnamed: 0', 'Genes')
    # filtered_markers = df[(df.avg_logFC > fc) &
    #                      (df.p_val_adj <= pval)]['Genes']
    # filtered_markers = df[(df.avg_logFC > fc)]['Genes']
    filtered_df = l_df[l_df['avg_logFC'] > 0.25][['Genes', 'avg_logFC']]
    filtered_df = filtered_df.sort_values(by='avg_logFC', ascending=False)
    filtered_dict = {}
    for r, c in filtered_df.iterrows():
        filtered_dict[c[1]] = c[0]
    # Volcano plot of DE genes
    _plot_de_genes(l_df)
    # return set(filtered_markers)
    return filtered_dict


def get_de_product_list(de_enzyme_product_list,
                        de_enzyme_stmts):
    if len(de_enzyme_product_list) > 1:
        de_enzyme_product_list = pd.merge(de_enzyme_stmts, de_enzyme_product_list,
                                          on=['Enzyme', 'Interaction', 'product', 'logFC'],
                                          how="outer").fillna('')
        return de_enzyme_product_list.sort_values(by='logFC', ascending=False)

    elif len(de_enzyme_product_list) < 1:
        de_enzyme_product_list = de_enzyme_stmts
        return de_enzyme_product_list


def get_enzyme_product_interactions(df, de_en_df, receptors_in_data):
    hashes_by_gene_pair = defaultdict(set)
    seen_product = set()
    product_and_fc = defaultdict(set)
    for r, c in de_en_df.iterrows():
        if c[2] not in seen_product:
            product_and_fc[c[2]].add((c[0], c[3]))
        seen_product.add(c[2])

    for a, b, hs in zip(df.agA_name, df.agB_name, df.stmt_hash):
        if a in product_and_fc and b in receptors_in_data:
            enzyme_logFC = [e for v in product_and_fc[a]
                            for e in v]
            enzyme, logFC = enzyme_logFC[0], enzyme_logFC[1]
            hashes_by_gene_pair[(a, b, enzyme, logFC)].add(hs)
    return hashes_by_gene_pair


def get_pain_phenotype(lg, pain_db):
    for r, c in pain_db.iterrows():
        if isinstance(pain_db.iloc[r]['gene_symbols'], str):
            l = set(pain_db.iloc[r]['gene_symbols'].split(","))
            pheno = pain_db.iloc[r]['phenotype_description']
            for g in l:
                if g in lg:
                    r_phenotype[(g)].add(pheno)
    return r_phenotype


def make_rnk(infile):
    df = pd.read_csv(infile, header=0, sep=",")
    df.columns = df.columns.str.replace('Unnamed: 0', 'Genes')
    df = df.loc[0:, ['Genes', 'p_val']]
    return df


def make_pheno_file(l_phenotype):
    pheno_df = []
    for keys, values in l_phenotype.items():
        pheno_df.append(
            {
                "Receptor": keys,
                "Phenotype_description": ", ".join(values)
            }
        )
    return pd.DataFrame(pheno_df)


def ligand_mgi_to_hgnc_name(seurat_ligand_genes):
    filtered_mgi = defaultdict(set)
    for logfc, gene in seurat_ligand_genes.items():
        if gene in mouse_gene_name_to_mgi:
            filtered_mgi[(gene, logfc)].add(mouse_gene_name_to_mgi[gene])

    hgnc_gene_dict = defaultdict(set)
    seen_genes = set()
    for key, value in filtered_mgi.items():
        mgi_id = next(iter(value))
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_symbol = get_hgnc_name(hgnc_id)
        if hgnc_symbol not in seen_genes:
            hgnc_gene_dict[(key[1])].add(hgnc_symbol)
        else:
            pass
        seen_genes.add(hgnc_symbol)
    return hgnc_gene_dict


def mgi_to_hgnc_name(gene_list):
    """Convert given mouse gene symbols to HGNC equivalent symbols"""
    filtered_mgi = {mouse_gene_name_to_mgi[gene] for gene in gene_list
                    if gene in mouse_gene_name_to_mgi}
    hgnc_gene_set = set()
    for mgi_id in filtered_mgi:
        hgnc_id = get_hgnc_from_mouse(mgi_id)
        hgnc_gene_set.add(get_hgnc_name(hgnc_id))
    return hgnc_gene_set


def make_interaction_df(interaction_dict):
    interaction_list = [
        {
            'Agent_A': [stmt[0].agent_list()][0][0].name,
            'Agent_B': [stmt[0].agent_list()][0][1].name,
            'Interaction type': re.match("\w+", str(stmt[0])).group(),
            'Enzyme': stmt[1],
            'logFC': fc
        }
        for fc, stmts in interaction_dict.items()
        for stmt in stmts
        if len(stmt[0].agent_list()) > 1

    ]
    df = pd.DataFrame(interaction_list)
    df = df.sort_values(by=['logFC'],
                        ascending=False)
    return df


def create_interaction_digraph(ligand_receptors,
                               sorted_enzyme_FC,
                               fname):
    '''
    This function takes two dictionaries as input,
    ligand receptors and enzyme fold change and creates
    a interaction Digraph of ligands, enzymes and receptors.

    Parameters
    ----------
    celtype_stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.
    network_name : Optional[str]
        The name of the network to be assembled. Default: indra_assembled

    Attributes
    ----------
    ligands_dict : dict
        Dict of foldchange and ligands as keys and receptors as values
    enzyme dict : dict
        Dict of foldchange as keys and enzymes as values
    fname : str
        output file name
    '''

    ligand_receptors = dict(sorted(ligand_receptors.items(),
                                   reverse=True))
    G = networkx.DiGraph()

    top_lg_rc = dict(sorted(itertools.islice(ligand_receptors.items(), 10)))
    top_en = dict(itertools.islice(sorted_enzyme_FC.items(), 10))

    for FC_lg, rcs in top_lg_rc.items():
        for rc in rcs:
            G.add_node(FC_lg[1], color='green')
            G.add_edge(FC_lg[1], rc, label="{:.2f}".format(FC_lg[0]))
    for en_FC, en in top_en.items():
        for chem in enzyme_product_dict[en]:
            for rcs in products_receptors[chem]:
                G.add_node(en, color='red')
                G.add_edge(en, chem, label="{:.2f}".format(en_FC))
                G.add_edge(chem, rcs)

    G.graph.setdefault('graph', {})['rankdir'] = 'LR'
    ag = networkx.nx_agraph.to_agraph(G)
    fname = os.path.join(OUTPUT, fname + "interactions_digraph.pdf")
    ag.draw(fname, prog='dot')


if __name__ == '__main__':
    # Read and extract cell surface proteins from CSPA DB
    wb = openpyxl.load_workbook(SURFACE_PROTEINS_WB)
    surface_protein_set = set(row[4].value for row in wb['Sheet 1']
                              if row[6].value == 'yes')
    logger.info('Got %d surface proteins from spreadsheet' %
                len(surface_protein_set))
    ligand_terms = ['cytokine activity', 'hormone activity',
                    'growth factor activity']
    receptor_terms = ['signaling receptor activity']

    ligand_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                     for term in ligand_terms]
    receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                       for term in receptor_terms]

    ligand_genes_go = get_genes_for_go_ids(ligand_go_ids)
    receptor_genes_go = get_genes_for_go_ids(receptor_go_ids)

    manual_ligands = {'THBS1'}

    # remove all the receptors from the surface_protein_set
    full_ligand_set = \
        (surface_protein_set - receptor_genes_go) | ligand_genes_go | \
        manual_ligands

    # Filtering out the nuclear receptors from the receptor list
    receptor_genes_go = filter_nuclear_receptors(receptor_genes_go,
                                                 'GO:0004879')

    # Add ION channels to the receptor list
    ion_channels = set()
    with open(ION_CHANNELS, 'r') as fh:
        for line in fh:
            ion_channels.add(line.strip())
    receptor_genes_go |= ion_channels

    # Load the INDRA DB DF
    df = load_indra_df(INDRA_DB_PKL)

    # Fetch omnipath database biomolecular interactions and
    # process them into INDRA statements
    op = process_from_web()

    ### Small molecule search
    if not os.path.exists(os.path.join(INPUT,
                                       'stmts_inhibition.pkl')):
        # Process TAS statements
        tp = tas.process_from_web()

        # Read drugbank database pickle
        with open(DRUG_BANK_PKL, "rb") as fh:
            dp = pickle.load(fh)

        # Run preassembly on a list of statements
        stmts = ac.run_preassembly(tp.statements + dp, return_toplevel=False,
                                   run_refinement=False)

        # Filter the statements to a given statement type
        stmts_inhibition = ac.filter_by_type(stmts, 'Inhibition')
    else:
        with open(os.path.join(INPUT, 'stmts_inhibition.pkl'), 'rb') as fh:
            stmts_inhibition = pickle.load(fh)

    targets_by_drug = defaultdict(set)

    # Create a dictionary of Drugs and targets
    for stmt in stmts_inhibition:
        drug_grounding = stmt.subj.get_grounding(
            ns_order=default_ns_order + ['CHEMBL', 'PUBCHEM', 'DRUGBANK',
                                         'HMS-LINCS'])
        targets_by_drug[(stmt.subj.name, drug_grounding)].add(stmt.obj.name)

    # Collect lists of receptors based on GO annotations and
    # by reading the data
    _, raw_receptor_genes = read_workbook(DATA_SPREADSHEET)
    receptor_genes = mgi_to_hgnc_name(raw_receptor_genes)
    receptors_in_data = receptor_genes & receptor_genes_go
    with open(os.path.join(OUTPUT, "receptors.csv"), 'w') as fh:
        fh.write('\n'.join(sorted(receptors_in_data)))

    # all_enzymes = get_controller_enzymes(['CHEBI:26333', 'CHEBI:3165'])
    all_enzymes = get_all_enzymes()

    stmts_by_cell_type = {}
    stmts_db_by_cell_type = {}
    num_interactions_by_cell_type = {}
    ligand_interactions_by_cell_type = {}
    possible_drug_targets = set()
    possible_db_drug_targets = set()
    de_enzyme_list = set()
    de_enzyme_product_list = set()
    r_phenotype = defaultdict(set)
    ligands_df = pd.DataFrame(columns=['Genes', 'p_val'])
    all_ranked_lg_df = pd.DataFrame(columns=['Interaction statement',
                                             'logFC'])
    enzyme_product_dict = defaultdict(set)
    possible_rc_drug_targets = set()
    possible_en_drug_targets = defaultdict(set)
    ligands_FC = {}
    enzymes_FC = {}
    cell_type_markers = {}


    # Looping over each file (cell type) and perform anylysis
    # for each cell type
    seurat_ligand_genes = {}
    for cell_type in IMMUNE_CELLTYPE_LIST:
        output_dir = os.path.join(OUTPUT, cell_type)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # read the input (immune cell type) ligand file
        cell_type_full = 'TwoGroups_DEG1_%s_AJ' % cell_type
        LIGANDS_INFILE = os.path.join(INPUT, '%s.csv' % cell_type_full)

        # Extract markers from seurat dataframe with logFC >= 0.25 and
        # pval <= 0.05
        seurat_ligand_genes = process_seurat_csv(LIGANDS_INFILE,
                                                 fc=0.25)

        # Pool all ligands along with its respective logFC and create a rank file
        # for GSEA
        ligands_df = pd.concat([ligands_df, make_rnk(LIGANDS_INFILE)],
                               ignore_index=True)

        if len(seurat_ligand_genes) == 0:
            logger.info('Skipping %s' % cell_type)
            continue

        # Get logFC as key and ligands as values
        ligand_genes = ligand_mgi_to_hgnc_name(seurat_ligand_genes)

        # Retain only ligands
        ligands_in_data = {k: next(iter(v)) for k, v in ligand_genes.items()
                           if next(iter(v)) in full_ligand_set}

        # Keep all ligands with FC from all cell types
        ligands_FC.update(ligands_in_data)

        # Retain only enzymes
        enzymes_in_data = {k: next(iter(v)) for k, v in ligand_genes.items()
                           if next(iter(v)) in all_enzymes}

        cell_type_markers[cell_type] = ligands_in_data
        cell_type_markers[cell_type].update(enzymes_in_data)

        # Keep all enzymes with FC from all cell typesa
        enzymes_FC.update(enzymes_in_data)

        # Get enzyme products by taking pathway commons DB as reference
        de_enzyme_stmts = enzyme_client.get_enzyme_products(enzymes_in_data)

        # Store enzymes as keys and its respective product
        # as values
        for r, c in de_enzyme_stmts.iterrows():
            enzyme_product_dict[(c[0])].add(c[2])

        # Keep merging enzyme products from all the celltypes
        de_enzyme_product_list = get_de_product_list(de_enzyme_product_list,
                                                     de_enzyme_stmts)

        # Get enzyme interactions with pain molecules
        pain_interactions = enzyme_client.get_pain_interactions(de_enzyme_stmts,
                                                                PAIN_MOL_NAMES)

        # write de enzyme products to a file
        de_enzyme_stmts.to_csv(os.path.join(output_dir, cell_type + "_de_enzymes_stmts.tsv"),
                               sep="\t", header=True, index=False)

        # write de enzyme pain interactions to a file
        pain_interactions.to_csv(os.path.join(output_dir, cell_type + "_enzyme_pain_interactions.tsv"),
                                 sep="\t", header=True, index=False)

        # Get the union of all the enzymes in the data
        possible_drug_targets |= set(enzymes_in_data.values())

        # Get the union of all the enzymes in the data
        de_enzyme_list |= set(enzymes_in_data.values())

        ligands_fc_df = {'Ligands': [*ligands_in_data.values()],
                         'logFC': [*ligands_in_data.keys()]}
        ligands_fc_df = pd.DataFrame(ligands_fc_df)

        ligands_fc_df.to_csv(os.path.join(output_dir, cell_type + "_ligands_fc.csv"),
                             header=True, index=False)

        logger.info(f'Loaded {len(ligands_in_data)} ligand genes from data')
        logger.info(f'Loaded {len(receptors_in_data)} receptor genes from data')

        # Now get INDRA DB Statements for the receptor-ligand pairs
        hashes_by_gene_pair = get_hashes_by_gene_pair(df, ligands_in_data,
                                                      receptors_in_data)

        all_hashes = set.union(*hashes_by_gene_pair.values())

        stmts_by_hash = download_statements(all_hashes)

        indra_db_stmts = list(stmts_by_hash.values())

        # Filtering out the indirect INDRA statements
        indra_db_stmts = ac.filter_direct(indra_db_stmts)

        # Filter statements which are not ligands/receptors
        op_filtered = filter_op_stmts(op.statements, ligands_in_data.values(),
                                      receptors_in_data)

        # Merge omnipath/INDRA statements and run assembly
        indra_op_stmts = ac.run_preassembly(indra_db_stmts + op_filtered,
                                            run_refinement=False)

        # Filter incorrect curations
        db_curations = get_curations()
        indra_op_filtered = ac.filter_by_curation(indra_op_stmts,
                                                  curations=db_curations)
        indra_op_filtered = filter_complex_statements(indra_op_filtered,
                                                      ligands_in_data,
                                                      receptors_in_data)

        # We do this again because when removing complex members, we
        # end up with more duplicates
        indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                               run_refinement=False)

        stmts_public = filter_out_medscan(indra_op_filtered)

        with open(os.path.join(
                OUTPUT, cell_type,
                'indra_ligand_receptor_statements.pkl'), 'wb') as fh:
            pickle.dump(indra_op_filtered, fh)

        stmts_by_cell_type[cell_type] = indra_op_filtered
        stmts_db_by_cell_type[cell_type] = filter_db_only(indra_op_filtered)
        num_interactions_by_cell_type[cell_type], \
        ligand_interactions_by_cell_type[cell_type] = \
            get_cell_type_stats(stmts_db_by_cell_type[cell_type],
                                ligands_in_data.values(),
                                receptors_in_data)

        # Creating a dict of logFC as key and
        # ligand as its value
        lg_logFC = {lg_fc[2]: lg_fc[0]
                    for lg_fc in hashes_by_gene_pair.keys()}

        logFC_stmts = defaultdict(set)
        lg_list = list(lg_logFC.values())
        fc_list = list(lg_logFC.keys())

        # Create a dict of logFC as key and statements as values
        for stmt in stmts_public:
            for ag in stmt.agent_list():
                if ag.name in lg_list:
                    fc = fc_list[lg_list.index(ag.name)]
                    logFC_stmts[(fc)].add(stmt)

        # Sorting the values in descending order
        sorted_lg_stmts = dict(sorted(logFC_stmts.items(),
                                      reverse=True))

        # create a dataframe of ranked statements
        sorted_stmts_df = [
            {
                'Interaction statement': stmt,
                'logFC': fc
            }
            for fc, stmts in sorted_lg_stmts.items()
            for stmt in stmts

        ]
        ranked_df = pd.DataFrame(sorted_stmts_df)
        all_ranked_lg_df = pd.concat([all_ranked_lg_df, ranked_df],
                                     ignore_index=False)

        # write the dataframe to a TSV file
        ranked_df.to_csv(os.path.join(output_dir, cell_type + "_ligand_receptor_ranked_stmts.tsv"),
                         sep="\t", header=True, index=False)

        # unpack the set of ranked statemetns into a list
        # for assembling into a html file
        sorted_lg_stmts_list = [stmt for stmts in sorted_lg_stmts.values()
                                for stmt in stmts]

        # Assemble the statements into HTML formatted report and save into a file
        indra_op_html_report = \
            html_assembler(
                sorted_lg_stmts_list,
                fname=os.path.join(output_dir,
                                   'indra_ligand_receptor_report.html'))

        # Assemble the statements into Cytoscape networks and save the file
        # into the disk
        # Optional: Please configure the indra config file in
        # ~/.config/indra/config.ini with NDEx credentials to upload the
        # networks into the server

        indra_op_cx_report, ndex_network_id = \
            cx_assembler(
                stmts_public,
                fname=os.path.join(output_dir,
                                   'indra_ligand_receptor_report.cx'))

        ligands_by_receptor = get_ligands_by_receptor(receptors_in_data,
                                                      set(ligands_in_data.values()),
                                                      indra_op_filtered)

        ligands_by_receptor_db = get_ligands_by_receptor(receptors_in_data,
                                                         set(ligands_in_data.values()),
                                                         stmts_db_by_cell_type[cell_type])

        for fc, en in enzymes_in_data.items():
            possible_en_drug_targets[(fc)].add(en)

        possible_rc_drug_targets |= set(ligands_by_receptor.keys())
        possible_drug_targets |= set(ligands_by_receptor.keys())
        possible_db_drug_targets |= set(ligands_by_receptor_db.keys())

    # Save all the ligand genes into a ranked files
    ligands_df.to_csv(os.path.join(OUTPUT, "ligands_pval.rnk"),
                      index=False, sep="\t")

    get_small_mol_report(targets_by_drug, possible_drug_targets,
                         os.path.join(OUTPUT, 'drug_targets.tsv'))

    get_small_mol_report(targets_by_drug, possible_db_drug_targets,
                         os.path.join(OUTPUT, 'drug_targets_db.tsv'))

    plot_interaction_potential(num_interactions_by_cell_type,
                               os.path.join(OUTPUT,
                                            'interaction_potential.pdf'))

    # Save the DE enzyme list to a file
    with open(os.path.join(OUTPUT,
                           'human_de_enzyme_list.txt'), 'w') as fh:
        fh.writelines("%s\n" % enzyme for enzyme in de_enzyme_list)

    # Save the DE enzyme product list to a csv file
    de_enzyme_product_list.to_csv(os.path.join(OUTPUT, "de_enzyme_products.csv"),
                                  index=False)

    # Get interactions for enzyme products expressed
    # in neurons
    de_enzyme_product_hash = get_enzyme_product_interactions(df, de_enzyme_product_list,
                                                             receptors_in_data)
    all_hashes = set.union(*de_enzyme_product_hash.values())

    # prdct_logFC stores enzyme product as a key and
    # its respective enzyme and logFC as value
    prdct_logFC = defaultdict(set)
    for prdct in de_enzyme_product_hash.keys():
        prdct_logFC[(prdct[0])].add((prdct[2], prdct[3]))

    stmts_by_hash = download_statements(all_hashes)
    indra_db_stmts = list(stmts_by_hash.values())

    # Filtering out the indirect INDRA statements
    indra_db_stmts = ac.filter_direct(indra_db_stmts)

    # Filter incorrect curations
    db_curations = get_curations()
    indra_op_filtered = ac.filter_by_curation(indra_db_stmts,
                                              curations=db_curations)
    indra_op_filtered = filter_complex_statements(indra_op_filtered,
                                                  set(de_enzyme_product_list['product']),
                                                  receptors_in_data)

    # We do this again because when removing complex members, we
    # end up with more duplicates
    indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                           run_refinement=False)

    # Creating a dictionary of logFC and
    # its respective statement and enzyme
    logFC_stmts = defaultdict(set)
    products_receptors = defaultdict(set)
    for stmt in indra_op_filtered:
        for ag in stmt.agent_list():
            if ag.name in prdct_logFC:
                for k in prdct_logFC[ag.name]:
                    en, fc = k[0], k[1]
                logFC_stmts[(fc)].add((stmt, en))
                chem, rc = stmt.agent_list()[0].name, stmt.agent_list()[1].name
                products_receptors[(chem)].add(rc)

    sorted_stmts = dict(sorted(logFC_stmts.items(), reverse=True))
    sorted_stmts_list = [stmt[0] for stmts in sorted_stmts.values()
                         for stmt in stmts]

    # Assemble the statements into HTML formatted report and save into a file
    indra_op_html_report = \
        html_assembler(
            sorted_stmts_list,
            fname=os.path.join(OUTPUT,
                               'logfc_indra_enzyme_product_neuron_report.html'))

    # Creating a dataframe of ranked enzyme products
    # and its interactions on the neuron side
    en_prdct_interaction_df = \
        make_interaction_df(sorted_stmts)

    # Write out the dataframe to a csv
    pd.DataFrame(sorted_stmts_df).to_csv(os.path.join(OUTPUT,
                                                      "ranked_de_enzyme_products_ineractions.csv"),
                                         index=True)

    # Make a dataframe of ligand
    # receptor interactions
    ranked_lg_dict = defaultdict(set)
    for r, c in all_ranked_lg_df.iterrows():
        ranked_lg_dict[(c[1])].add((c[0], 'NA'))

    lg_interaction_df = \
        make_interaction_df(ranked_lg_dict)

    # Concat enzyme product dataframe and
    # ligand receptor interaction dataframe
    full_interaction_df = pd.concat([en_prdct_interaction_df,
                                     lg_interaction_df]).sort_values(by=['logFC'],
                                                                     ascending=False)
    full_interaction_df.to_csv(os.path.join(OUTPUT, 'full_ranked_interaction.csv'),
                               index=True, sep=",", header=True)

    human_pain_df = pd.read_csv(HUMAN_PAIN_DB, sep="\t", header=0)
    r_phenotype_dict = get_pain_phenotype(receptors_in_data,
                                          human_pain_df)

    # Make the phenotype file and save it
    r_phenotype_csv = make_pheno_file(r_phenotype_dict)
    pd.DataFrame(r_phenotype_csv).to_csv(os.path.join(OUTPUT, "receptor_phenotypes.csv"),
                                         index=False, header=True, sep=",")

    # full drug interactions
    possible_en_drug_targets = OrderedDict(sorted(possible_en_drug_targets.items(),
                                                  reverse=True))
    fc_en_drug_targets = list(possible_en_drug_targets.keys())
    en_drug_targets = list(possible_en_drug_targets.values())

    en_targets_in_data = defaultdict(set)
    rc_targets_in_data = defaultdict(set)
    fc_in_data = defaultdict(set)

    for drug, targets in targets_by_drug.items():
        for t in targets:
            for e in en_drug_targets:
                if t in e:
                    fc = fc_en_drug_targets[en_drug_targets.index(e)]
                    fc_in_data[(drug[0])].add(fc)
                    en_targets_in_data[(drug[0])].add(t)
                if t in ligands_by_receptor.keys():
                    rc_targets_in_data[(drug[0])].add(t)

    drug_interaction_list = []
    for drug, targets in targets_by_drug.items():
        intermediates = set()
        l = []
        if drug[0] in en_targets_in_data.keys():
            en_target = list(en_targets_in_data[drug[0]])
            for enzymes in en_target:
                intermediates.update(enzyme_product_dict[enzymes])
            intermediates = list(intermediates)

            for val in intermediates:
                if val != None:
                    l.append(val)
            avgFC = sum(fc_in_data[drug[0]]) / len(fc_in_data[drug[0]])
        else:
            en_target = ''
            avgFC = 0

        if drug[0] in rc_targets_in_data.keys():
            rc_target = rc_targets_in_data[drug[0]]
        else:
            rc_target = ""

        other_targets = targets - set(en_target) - set(rc_target)

        drug_interaction_list.append(
            {
                "Drug": drug[0],
                "Named": 0 if drug[0].startswith('CHEMBL') else 1,
                "Enzyme targets": ", ".join(en_target),
                "Enzyme products": ", ".join(l),
                "Receptor targets": ", ".join(rc_target),
                "Score": "{:.3f}".format((len(en_target) + len(rc_target)) / len(targets)),
                "Other targets": ", ".join(sorted(other_targets)),
                "Other target hits": "{:.3f}".format(len(other_targets)),
                "Total hits": "{:.3f}".format((len(en_target) + len(rc_target) + len(other_targets))),
                "avgFC": float(avgFC)
            }
        )
    pd.DataFrame(drug_interaction_list).sort_values(by=['avgFC'], ascending=False).to_csv(os.path.join(
        OUTPUT, "ranked_enzyme_drug_target.csv"))

    # Di-graph of interactions between top 10 DE ligands, enzyme products
    # and receptors
    ligands_logFC = defaultdict(set)
    ligandsFC_by_receptor = defaultdict(set)
    sorted_ligands_FC = dict(sorted(ligands_FC.items(), reverse=True))
    sorted_enzyme_FC = dict(sorted(enzymes_FC.items(), reverse=True))
    lg_fc = list(sorted_ligands_FC.keys())
    lg = list(sorted_ligands_FC.values())
    filtered_stmts_by_cell_type = []
    cell_type_stmts = defaultdict(set)

    for cell_type in stmts_db_by_cell_type.keys():
        stmts = stmts_db_by_cell_type[cell_type]
        for stmt in stmts:
            agent_names = {agent.name for agent in stmt.agent_list()}
            receptors = agent_names & receptors_in_data
            ligands = "".join(agent_names & set(lg))
            if len(ligands) > 0 and ligands not in receptors:
                filtered_stmts_by_cell_type.append(stmt)
                en_logFC = lg_fc[lg.index(ligands)]
                for receptor in receptors:
                    # Storing the interactions in a dictionary
                    # for each cell type and plot di-graphs
                    cell_type_stmts[(en_logFC, ligands)].add(receptor)

                    # Storing interactions for all the cell-types
                    ligandsFC_by_receptor[(en_logFC, ligands)].add(receptor)

        # Plot interactions for each cell type
        create_interaction_digraph(cell_type_stmts, sorted_enzyme_FC,
                                   fname=os.path.join(cell_type, cell_type + "_"))
        # Reset the cell type statements dict
        cell_type_stmts = defaultdict(set)

    filtered_stmts_by_cell_type = ac.run_preassembly(filtered_stmts_by_cell_type,
                                                     run_refinement=False)
    # Assemble the statements into HTML formatted report and save into a file
    indra_op_html_report = \
        html_assembler(
            filtered_stmts_by_cell_type,
            fname=os.path.join(OUTPUT,
                               'filtered_ligand_enzyme_interactions.html'))

    # Plotting interaction Di-graph for all the cell-types
    create_interaction_digraph(ligandsFC_by_receptor, sorted_enzyme_FC, fname='')

    ### Downstream analysis

    # Human pain genetics db
    gene_symbol = [str(g).split(",") for g in human_pain_df['gene_symbols']]
    pain_phenotype = [str(p).split(",") for p in human_pain_df['phenotype_description']]
    pheno_direction = human_pain_df['phenotype_direction']
    pain_pheno_dict = defaultdict(set)

    # Import ctd database
    # Opening JSON file and loading the data
    # into the variable data
    with open(os.path.join(INPUT, 'ctd_db.json')) as json_file:
        data = json.load(json_file)

    ctd_db_json = data['associations']

    # now we will open a file for writing
    data_file = open(os.path.join(INPUT, 'ctd_db.csv'), 'w')

    # create the csv writer object
    csv_writer = csv.writer(data_file)

    # Counter variable used for writing
    # headers to the CSV file
    count = 0

    for emp in ctd_db_json:
        if count == 0:
            # Writing headers of CSV file
            header = emp.keys()
            csv_writer.writerow(header)
            count += 1

        # Writing data of CSV file
        csv_writer.writerow(emp.values())
    data_file.close()

    # Read the csv back and clean it
    ctd_csv = pd.read_csv(os.path.join(INPUT, 'ctd_db.csv'), header=0)
    count = 0
    # Clean the gene column
    for r, c in ctd_csv.iterrows():
        gene = re.search("symbol': '([A-Z0-9]+)", c[0]).group(1)
        ctd_csv['gene'][count] = gene
        count += 1

    # Creating a dictionary of ctd genes
    # as keys and its respective Standarized values
    # as values
    ctd_geneset = dict()
    for r, c in ctd_csv.iterrows():
        ctd_geneset[c[0]] = c[2]

    with open(INDRA_DB_PKL, 'rb') as indra_sif:
        indra_sif = pickle.load(indra_sif)

    downstream_hits = defaultdict(set)
    downstream_evidence = dict()
    downstream_stmts = defaultdict(set)
    rc_in_data_dwnstrm = set()
    upstream = defaultdict(set)
    not_upstream = defaultdict(set)

    for a, b, ev, hs in zip(indra_sif.agA_name,
                            indra_sif.agB_name,
                            indra_sif.evidence_count,
                            indra_sif.stmt_hash):

        upstream[(b)].add(a)
        if a in receptor_genes_go:
            downstream_hits[(b)].add(a)
            if b not in downstream_evidence:
                downstream_evidence[b] = ev

            else:
                downstream_evidence[b] = downstream_evidence[b] + ev

            downstream_stmts[(b)].add(hs)
        if a in receptors_in_data:
            rc_in_data_dwnstrm.add(b)

    targets_keys = list(targets_by_drug.keys())
    targets_values = list(targets_by_drug.values())

    ds_in_data = defaultdict(set)
    drug_names = defaultdict(set)
    drug_no_names = defaultdict(set)

    for drug, targets in targets_by_drug.items():
        for t in targets:
            if t in downstream_hits.keys():
                if drug[0].startswith('CHEMBL'):
                    drug_no_names[(t)].add(drug[0])
                else:
                    drug_names[(t)].add(drug[0])

    # Below code is to aggregate all the
    # downstream hits and its pain pheno direction
    ds_direction_dict = dict()
    pheno_direction = human_pain_df['phenotype_direction']
    for k in downstream_hits.keys():
        for g, p, pdir in zip(gene_symbol,
                              pain_phenotype,
                              pheno_direction):
            if k in g:
                p = "".join(p)
                if k + ":" + p + ":" + pdir in ds_direction_dict:
                    ds_direction_dict[k + ":" + p + ":" + pdir] = ds_direction_dict[k + ":" + p + ":" + pdir] + 1
                else:
                    ds_direction_dict[k + ":" + p + ":" + pdir] = 1
    d = defaultdict(set)
    for t in ds_direction_dict:
        d[re.match("([a-zA-Z0-9-]+)(.*)", t).group(1)].add(
            re.match("([a-zA-Z0-9-]+:)(.*)", t).group(2) + ":" + str(ds_direction_dict[t]))

    # Creating a dictionary of downstream targets as keys
    # its values as pain phenotype description
    for k in downstream_hits.keys():
        for g, p in zip(gene_symbol, pain_phenotype):
            if k in g:
                pain_pheno_dict[(k)].add("".join(p))

    downstream_df = []
    for ds, us in downstream_hits.items():
        if ds in receptor_genes:
            total_interactome = len(receptors_in_data)
            total_non_interactome = len(receptor_genes_go - receptors_in_data)
            interactome_genes = len(receptors_in_data & downstream_hits[ds])
            non_interactome_genes = len(downstream_hits[ds] - receptors_in_data)

            pval = fisher_exact([[interactome_genes, non_interactome_genes],
                                 [total_interactome, total_non_interactome]],
                                alternative='greater')[1]

            if ds in drug_names:
                drugs_names = ", ".join(drug_names[ds])
            else:
                drugs_names = 'NA'

            if ds in drug_no_names:
                drugs_no_names = ", ".join(drug_no_names[ds])
            else:
                drugs_no_names = 'NA'

            if ds in ctd_geneset:
                ctd_gene = 1
                std_val = ctd_geneset[ds]
            else:
                ctd_gene = 0
                std_val = 0

            if ds in d:
                pain_phenotype_des = ", ".join(d[ds])
            else:
                pain_phenotype_des = 'NA'

            downstream_df.append(
                {
                    "Upstream genes": ", ".join(us),
                    "Upstream in interacome": ", ".join(receptors_in_data & downstream_hits[ds]),
                    "Downstream target": ds,
                    "CTD_DB": ctd_gene,
                    "CTD_Std_Val": std_val,
                    "Pain description": pain_phenotype_des,
                    "Named drugs": drugs_names,
                    "Unnamed drugs": drugs_no_names,
                    "evidence count": downstream_evidence[ds],
                    "Total upstream": len(us),
                    "Total upstream in interactome": interactome_genes,
                    "pval": pval,
                    "Statement hash": ', '.join(str(x) for x in downstream_stmts[ds])
                }
            )
    downstream_df = pd.DataFrame(downstream_df).sort_values(by=['pval'],
                                                            ascending=True)

    downstream_df.to_csv(os.path.join(OUTPUT, "downstream_hits_filtered_pval.csv"),
                         sep=",", header=True, index=False)

    # Creating layout aspect

    ligandsFC_by_receptor = dict(sorted(ligandsFC_by_receptor.items(),
                                        reverse=True))

    G = Digraph(engine='dot')

    top_lg_rc = dict(sorted(itertools.islice(ligandsFC_by_receptor.items(), 10)))
    top_en = dict(itertools.islice(sorted_enzyme_FC.items(), 10))

    for FC_lg, rcs in top_lg_rc.items():
        for rc in rcs:
            G.node(FC_lg[1], color='green')
            G.edge(FC_lg[1], rc, label="{:.2f}".format(FC_lg[0]))
    for en_FC, en in top_en.items():
        for chem in enzyme_product_dict[en]:
            for rcs in products_receptors[chem]:
                G.node(en, color='red')
                G.edge(en, chem, label="{:.2f}".format(en_FC))
                G.edge(chem, rcs)

    G.graph_attr['rankdir'] = 'LR'

    # decode the Digraph to JSON format
    graphviz_json_string = G.pipe('json').decode()

    # parse the resulting json_string
    graphviz_json_dict = json.loads(graphviz_json_string)

    # Assembling the top interactions in the interactome into
    # INDRA statements
    top_lg_rc = dict(sorted(itertools.islice(ligandsFC_by_receptor.items(), 10)))
    top_en = dict(itertools.islice(sorted_enzyme_FC.items(), 10))

    ## Get enzyme product statements
    en_product = defaultdict(set)
    en_fc = dict()
    en_fc_stmts = defaultdict(set)

    for k, v in top_en.items():
        for chem in enzyme_product_dict[v]:
            if len(products_receptors[chem]) >= 1:
                en_product[(v)].add(chem)
                en_fc[v] = k

    # get statements for the enzyme products
    for a, b, hs in zip(df.agA_name,
                        df.agB_name,
                        df.stmt_hash):
        if a in en_product and b in en_product[a]:
            en_fc_stmts[(en_fc[a], a, b)].add(hs)

    ## Get ligand receptor statements
    lg_fc = dict()
    lg_rc = defaultdict(set)
    lg_rc_stmts = defaultdict(set)
    for k, v in top_lg_rc.items():
        lg_fc[k[1]] = k[0]
        lg_rc[(k[1])].update(v)

    for a, b, hs in zip(df.agA_name,
                        df.agB_name,
                        df.stmt_hash):
        if a in lg_rc and b in lg_rc[a]:
            lg_rc_stmts[(lg_fc[a], a, b)].add(hs)

    ## Get enzyme product receptor statements
    chem_rc = defaultdict(set)
    chem_stmts = defaultdict(set)

    for en_FC, en in top_en.items():
        for chem in enzyme_product_dict[en]:
            for rcs in products_receptors[chem]:
                chem_rc[(chem)].add(rcs)

    # get statements for the enzyme products and receptors
    for a, b, hs in zip(df.agA_name,
                        df.agB_name,
                        df.stmt_hash):
        if a in chem_rc and b in chem_rc[a]:
            chem_stmts[(a, b)].add(hs)

    # Concatenating all the statements into a list
    agents = set(list(lg_rc.keys()) + list(chem_rc.keys()) + list(en_product.keys()) \
                 + [r for v in lg_rc.values() for r in v] + [r for v in chem_rc.values() for r in v] + \
                 [r for v in en_product.values() for r in v])

    en_hs = [h for s in en_fc_stmts.values() for h in s]
    lg_hs = [h for s in lg_rc_stmts.values() for h in s]
    chem_rc_hs = [h for s in chem_stmts.values() for h in s]
    merged_hs = en_hs + lg_hs + chem_rc_hs

    stmts_by_hash = download_statements(merged_hs)
    indra_db_stmts = list(stmts_by_hash.values())

    # Filtering out the indirect INDRA statements
    indra_db_stmts = ac.filter_direct(indra_db_stmts)
    indra_filtered = ac.run_preassembly(indra_db_stmts,
                                        run_refinement=False)

    # Remove agents from complex statements
    # which are not in the interactome
    for stmt in indra_filtered:
        if isinstance(stmt, Complex):
            # Statement updated by reference here
            stmt.members = [agent for agent in stmt.members
                            if agent.name in agents or agent.name in agents]

    # bens snippet
    cxa = CxAssembler(indra_filtered,
                      network_name='neuro_interactome')
    cxa.make_model()
    cx_node_name_to_id = {node['n']: node['@id']
                          for node in cxa.cx['nodes']}

    # Creating layout aspect
    layout_aspect = []
    for obj in graphviz_json_dict['objects']:
        cord = obj['pos'].split(',')
        node_name = obj['name']
        cx_node_id = cx_node_name_to_id.get(node_name)
        if not cx_node_id:
            logger.info('The node %s is not in the CX network' % node_name)
            continue
        layout_aspect.append(
            {
                'node': cx_node_id,
                'x': float(cord[0]),
                'y': float(cord[1])
            }
        )

    cxa.cx['cartesianLayout'] = layout_aspect
    ndex_network_id = cxa.upload_model(ndex_cred=None,
                                       private=True, style='default')