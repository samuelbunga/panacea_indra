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
db_curations = get_curations()


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

# Load the INDRA DB DF
df = load_indra_df(INDRA_DB_PKL)


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
    r_phenotype = defaultdict(set)
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


def filter_incorrect_curations(stmts):
    # Filter incorrect curations
    indra_op_filtered = ac.filter_by_curation(stmts,
                                              curations=db_curations)
    return indra_op_filtered


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
                               fname,
                               enzyme_product_dict,
                               products_receptors):
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

    top_lg_rc = dict(sorted(ligand_receptors.items()))
    top_en = dict(sorted_enzyme_FC.items())

    for FC_lg, rcs in top_lg_rc.items():
        for rc in rcs:
            G.add_node(FC_lg[1], color='green')
            G.add_edge(FC_lg[1], rc, label="{:.2f}".format(FC_lg[0]))
    for en_FC, en in top_en.items():
        for chem in enzyme_product_dict[en]:
            for rcs in products_receptors[chem]:
                print(rcs)
                G.add_node(en, color='red')
                G.add_edge(en, chem, label="{:.2f}".format(en_FC))
                G.add_edge(chem, rcs)

    G.graph.setdefault('graph', {})['rankdir'] = 'LR'
    ag = networkx.nx_agraph.to_agraph(G)
    fname = os.path.join(OUTPUT, fname + "interactions_digraph.pdf")
    ag.draw(fname, prog='dot')