from functions import *

# Read and extract cell surface proteins from CSPA DB
wb = openpyxl.load_workbook(SURFACE_PROTEINS_WB)
surface_protein_set = set(row[4].value for row in wb['Sheet 1']
                          if row[6].value == 'yes')
logger.info('Got %d surface proteins from spreadsheet' %
            len(surface_protein_set))
ligand_terms = ['cytokine activity', 'hormone activity',
                'growth factor activity']
receptor_terms = ['signaling receptor activity']

# Getting GO id's for ligands and receptors by using
# GO terms
ligand_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                 for term in ligand_terms]
receptor_go_ids = [bio_ontology.get_id_from_name('GO', term)[1]
                   for term in receptor_terms]

# Converting GO id's to gene symbols
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

    # cache the statemnts to pickle for later use
    with open(os.path.join(INPUT, 'stmts_inhibition.pkl'), 'wb') as fh:
        pickle.dump(stmts_inhibition, fh)

else:
    # Read the tas statements which are cached
    with open(os.path.join(INPUT, 'stmts_inhibition.pkl'), 'rb') as fh:
        stmts_inhibition = pickle.load(fh)

targets_by_drug = defaultdict(set)

# Create a dictionary of Drugs and targets
for stmt in stmts_inhibition:
    drug_grounding = stmt.subj.get_grounding(
        ns_order=default_ns_order + ['CHEMBL', 'PUBCHEM', 'DRUGBANK',
                                     'HMS-LINCS'])
    targets_by_drug[(stmt.subj.name, drug_grounding)].add(stmt.obj.name)

with open(os.path.join(OUTPUT, 'targets_by_drug.pkl'), 'wb') as fh:
    pickle.dump(targets_by_drug, fh)

# Collect lists of receptors based on GO annotations and
# by reading the data
# Read list of neuro immune genes from the spread sheet
_, raw_receptor_genes = read_workbook(DATA_SPREADSHEET)
receptor_genes = mgi_to_hgnc_name(raw_receptor_genes)
receptors_in_data = receptor_genes & receptor_genes_go
with open(os.path.join(OUTPUT, "receptors.csv"), 'w') as fh:
    fh.write('\n'.join(sorted(receptors_in_data)))

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
    logger.info('Processing %s' % cell_type)
    output_dir = os.path.join(OUTPUT, cell_type)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Checking for indra ligand statements pkl in the 
    # current directory
    if os.path.isfile(os.path.join(output_dir, 'indra_ligand_receptor_statements.pkl')):
        with open(os.path.join(output_dir, 'indra_ligand_receptor_statements.pkl'), 'rb') as fh:
            indra_op_filtered = pickle.load(fh)
    else:
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

        # Save enzyme_in_data object
        with open(os.path.join(output_dir, 'enzymes_in_data.pkl'), 'wb') as fh:
            pickle.dump(enzymes_in_data, fh)

        cell_type_markers[cell_type] = ligands_in_data
        cell_type_markers[cell_type].update(enzymes_in_data)

        # Keep all enzymes with FC from all cell typesa
        enzymes_FC.update(enzymes_in_data)

        # Get enzyme products by taking pathway commons DB as reference
        de_enzyme_product_interaction = enzyme_client.get_enzyme_products(enzymes_in_data)

        # Store enzymes as keys and its respective product
        # as values
        for r, c in de_enzyme_product_interaction.iterrows():
            enzyme_product_dict[(c[0])].add(c[2])

        # Keep merging enzyme products interactions from all the celltypes
        de_enzyme_product_list = get_de_product_list(de_enzyme_product_list,
                                                     de_enzyme_product_interaction)

        # Get enzyme interactions with pain molecules for the cell type
        pain_interactions = enzyme_client.get_pain_interactions(de_enzyme_product_interaction,
                                                                PAIN_MOL_NAMES)

        # write de enzyme products to a file
        de_enzyme_product_interaction.to_csv(os.path.join(output_dir, cell_type + "_de_enzymes_product_interaction.tsv"),
                               sep="\t", header=True, index=False)

        # write de enzyme pain interactions to a file
        pain_interactions.to_csv(os.path.join(output_dir, cell_type + "_enzyme_pain_interactions.tsv"),
                                 sep="\t", header=True, index=False)

        # Get the union of all the enzymes in the data
        # and include them in drug target set
        possible_drug_targets |= set(enzymes_in_data.values())

        # Get the union of all the enzymes in the data
        # from all the cell types
        de_enzyme_list |= set(enzymes_in_data.values())

        # Make a data frame of ligands and logFC
        # for the respective cell type
        ligands_fc_df = {'Ligands': [*ligands_in_data.values()],
                         'logFC': [*ligands_in_data.keys()]}
        ligands_fc_df = pd.DataFrame(ligands_fc_df)

        # Save the dataframe into a csv
        ligands_fc_df.to_csv(os.path.join(output_dir, cell_type + "_ligands_fc.csv"),
                             header=True, index=False)
        logger.info(f'Loaded {len(ligands_in_data)} ligand genes from data')
        logger.info(f'Loaded {len(receptors_in_data)} receptor genes from data')

        # Now get INDRA DB Statements for the receptor-ligand pairs
        hashes_by_gene_pair = get_hashes_by_gene_pair(df, ligands_in_data,
                                                      receptors_in_data)

        with open(os.path.join(output_dir, 'hashes_by_gene_pair.pkl'), 'wb') as fh:
            pickle.dump(hashes_by_gene_pair, fh)

        # get the union of all the statement hashes 
        all_hashes = set.union(*hashes_by_gene_pair.values())
        # Download the statements by hashes
        stmts_by_hash = download_statements(all_hashes)
        # get only the list of all the available statemtns
        indra_db_stmts = list(stmts_by_hash.values())
        # Filtering out the indirect INDRA statements
        indra_db_stmts = ac.filter_direct(indra_db_stmts)
        # Filter statements which are not ligands/receptors from 
        # OmniPath database
        op_filtered = filter_op_stmts(op.statements, ligands_in_data.values(),
                                      receptors_in_data)
        # Merge omnipath/INDRA statements and run assembly
        indra_op_stmts = ac.run_preassembly(indra_db_stmts + op_filtered,
                                            run_refinement=False)

        # Filter incorrect curations        
        indra_op_filtered = filter_incorrect_curations(indra_op_stmts)

        # Filter complex statements
        indra_op_filtered = filter_complex_statements(indra_op_filtered,
                                                      ligands_in_data,
                                                      receptors_in_data)

        # We do this again because when removing complex members, we
        # end up with more duplicates
        indra_op_filtered = ac.run_preassembly(indra_op_filtered,
                                               run_refinement=False)
        # Filter out medscan statements
        stmts_public = filter_out_medscan(indra_op_filtered)

        # Save indra ligand receptor statements for the
        # respective cell type
        with open(os.path.join(output_dir, 
                               'indra_ligand_receptor_statements.pkl'), 'wb') as fh:
            pickle.dump(stmts_public, fh)

    # Save enzyme possible drug targtes
    if not os.path.isfile(os.path.join(OUTPUT, 'possible_en_drug_targets.pkl')):
        with open(os.path.join(OUTPUT, 'possible_en_drug_targets.pkl'), 'wb') as fh:
            pickle.dump(possible_en_drug_targets, fh)

    # Load enzymes in data
    with open(os.path.join(output_dir, 'enzymes_in_data.pkl'), 'rb') as fh:
        enzymes_in_data = pickle.load(fh)

    # load hashes by gene pair
    with open(os.path.join(output_dir, 'hashes_by_gene_pair.pkl'), 'rb') as fh:
        hashes_by_gene_pair = pickle.load(fh)
    ligands_fc_df = pd.read_csv(os.path.join(output_dir, cell_type + "_ligands_fc.csv"), 
                                header=0)
    # Get the list of ligands from the cell type
    ligands_in_data = list(ligands_fc_df['Ligands'])

    stmts_public = indra_op_filtered

    # Store the cell type specific indra statements in a dictionary
    stmts_by_cell_type[cell_type] = indra_op_filtered

    # Filter statements to database only and store them
    # in a separate dictionary
    stmts_db_by_cell_type[cell_type] = filter_db_only(indra_op_filtered)

    # Getting cell type stats
    num_interactions_by_cell_type[cell_type], \
    ligand_interactions_by_cell_type[cell_type] = \
        get_cell_type_stats(stmts_db_by_cell_type[cell_type],
                            ligands_in_data,
                            receptors_in_data)

    # Creating a dict of logFC as key and
    # ligand as its value for the respective cell type
    # from the hashes of ligands and receptors
    lg_logFC = {lg_fc[2]: lg_fc[0]
                for lg_fc in hashes_by_gene_pair.keys()}

    logFC_stmts = defaultdict(set)
    fc_list = list(lg_logFC.keys())
    lg_list = list(lg_logFC.values())

    # Create a dict of logFC as key and statements as values
    for stmt in indra_op_filtered:
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
    # Concat all the ligands receptors interaction statements 
    # from all the cell types
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
            stmts_by_cell_type[cell_type],
            fname=os.path.join(output_dir,
                               'indra_ligand_receptor_report.cx'))

    # create a dictionary of receptors as keys and its repective
    # ligands as values
    ligands_by_receptor = get_ligands_by_receptor(receptors_in_data,
                                                  set(ligands_in_data),
                                                  stmts_by_cell_type[cell_type])
    with open(os.path.join(OUTPUT, 'ligands_by_receptor.pkl'), 'wb') as fh:
        pickle.dump(ligands_by_receptor, fh)

    # create a dictionary of receptors as keys and its repective
    # ligands as values from ligands and receptors with filtered database 
    # statements
    ligands_by_receptor_db = get_ligands_by_receptor(receptors_in_data,
                                                     set(ligands_in_data),
                                                     stmts_db_by_cell_type[cell_type])
    # Create a dictionary of enzymes and its foldchange
    for fc, en in enzymes_in_data.items():
        possible_en_drug_targets[(fc)].add(en)

    # Take a union of receptors from ligands by receptor dictionary
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

# Save ligands_FC
if not os.path.isfile(os.path.join(OUTPUT, "ligands_FC.pkl")):
    with open(os.path.join(OUTPUT, 'ligands_FC.pkl'), 'wb') as fh:
        pickle.dump(ligands_FC, fh)

# Save Enzymes_FC
if not os.path.isfile(os.path.join(OUTPUT, "enzymes_FC.pkl")):
    with open(os.path.join(OUTPUT, 'enzymes_FC.pkl'), 'wb') as fh:
        pickle.dump(enzymes_FC, fh)

# Save statements by DB
if not os.path.isfile(os.path.join(OUTPUT, "stmts_db_by_cell_type.pkl")):
    with open(os.path.join(OUTPUT, 'stmts_db_by_cell_type.pkl'), 'wb') as fh:
        pickle.dump(stmts_db_by_cell_type, fh)

# Save statements by cell type
if not os.path.isfile(os.path.join(OUTPUT, "stmts_by_cell_type.pkl")):
    with open(os.path.join(OUTPUT, 'stmts_by_cell_type.pkl'), 'wb') as fh:
        pickle.dump(stmts_by_cell_type, fh)

# Save enzyme product list
if not os.path.isfile(os.path.join(OUTPUT, "enzyme_product_dict.pkl")):
    with open(os.path.join(OUTPUT, 'enzyme_product_dict.pkl'), 'wb') as fh:
         pickle.dump(enzyme_product_dict, fh)

# Save the DE enzyme list to a file
with open(os.path.join(OUTPUT,
                       'human_de_enzyme_list.txt'), 'w') as fh:
    fh.writelines("%s\n" % enzyme for enzyme in de_enzyme_list)

# Save the DE enzyme product list to a csv file
if not os.path.isfile(os.path.join(OUTPUT, "de_enzyme_products.pkl")):
    with open(os.path.join(OUTPUT, 'de_enzyme_products.pkl'), 'wb') as fh:
        pickle.dump(de_enzyme_product_list, fh)

# Save ligand receptor interactions from all
# the cell types
if not os.path.isfile(os.path.join(OUTPUT, "all_ligand_receptor_statements.pkl")):
    with open(os.path.join(OUTPUT, 'all_ligand_receptor_statements.pkl'), 'wb') as fh:
        pickle.dump(all_ranked_lg_df, fh)
    all_ranked_lg_df.to_csv(os.path.join(OUTPUT, 'all_ligand_receptor_statements.csv'),
                            header=True,
                            index=False)
