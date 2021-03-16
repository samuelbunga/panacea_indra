from functions import *

### ENZYMES
with open(os.path.join(OUTPUT, 'de_enzyme_products.pkl'), 'rb') as fh:
    de_enzyme_product_list = pickle.load(fh)

with open(os.path.join(OUTPUT, "receptors.csv"), 'r') as fh:
    receptors_in_data = {l.strip() for l in fh}


# Get interactions for enzyme products expressed
# in neurons
de_enzyme_product_hash = get_enzyme_product_interactions(df, de_enzyme_product_list,
                                                         receptors_in_data)

# Get union of all the hashes
all_hashes = set.union(*de_enzyme_product_hash.values())

# prdct_logFC stores enzyme product as a key and
# its respective enzyme and logFC as value
prdct_logFC = defaultdict(set)
for prdct in de_enzyme_product_hash.keys():
    prdct_logFC[(prdct[0])].add((prdct[2], prdct[3]))

# Download the statements from hashes
stmts_by_hash = download_statements(all_hashes)
indra_db_stmts = list(stmts_by_hash.values())

# Filtering out the indirect INDRA statements
indra_db_stmts = ac.filter_direct(indra_db_stmts)
# Filtering out incorrect statemetns
indra_op_filtered = filter_incorrect_curations(indra_db_stmts)
# Filtering out the complex statements
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

# Sort the statements 
sorted_stmts = dict(sorted(logFC_stmts.items(), reverse=True))

# Create a list of sorted statements
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
en_prdct_interaction_df.to_csv(os.path.join(OUTPUT,
                                            "ranked_de_enzyme_products_ineractions.csv"),
                               index=True)

# Save Products receptors object
if not os.path.isfile(os.path.join(OUTPUT, 'products_receptors.pkl')):
    with open(os.path.join(OUTPUT, 'products_receptors.pkl'), 'wb') as fh:
        pickle.dump(products_receptors, fh)

### Ligand receptor interactions
# Make a dataframe of ligand
# receptor interactions
if os.path.isfile(os.path.join(OUTPUT, 'all_ligand_receptor_statements.pkl')):
    with open(os.path.join(OUTPUT, 'all_ligand_receptor_statements.pkl'), 'rb') as fh:
        all_ranked_lg_df = pickle.load(fh)
    ranked_lg_dict = defaultdict(set)
    for r, c in all_ranked_lg_df.iterrows():
        ranked_lg_dict[(c[1])].add((c[0], 'NA'))

    lg_rc_interaction_df = \
        make_interaction_df(ranked_lg_dict)

    # Concat enzyme product dataframe and
    # ligand receptor interaction dataframe
    full_interaction_df = pd.concat([en_prdct_interaction_df,
                                     lg_rc_interaction_df]).sort_values(by=['logFC'],
                                                                    ascending=False)
    #Write the full ligand, enzyme -> receptor to csv
    full_interaction_df.to_csv(os.path.join(OUTPUT, 'full_ranked_interaction.csv'),
                               index=True, sep=",", header=True)
    

### Enzyme product drug interaction
with open(os.path.join(OUTPUT, 'possible_en_drug_targets.pkl'), 'rb') as fh:
            possible_en_drug_targets = pickle.load(fh)
        
with open(os.path.join(OUTPUT, 'ligands_by_receptor.pkl'), 'rb') as fh:
    ligands_by_receptor = pickle.load(fh)
        
# full drug interactions
possible_en_drug_targets = OrderedDict(sorted(possible_en_drug_targets.items(),
                                              reverse=True))
# Load logFC of enzyme drug targets into a list
fc_en_drug_targets = list(possible_en_drug_targets.keys())

# Load enzymes into list
en_drug_targets = list(possible_en_drug_targets.values())

en_targets_in_data = defaultdict(set)
rc_targets_in_data = defaultdict(set)
fc_in_data = defaultdict(set)

with open(os.path.join(OUTPUT, 'targets_by_drug.pkl'), 'rb') as fh:
        targets_by_drug = pickle.load(fh)
        
for drug, targets in targets_by_drug.items():
    for t in targets:
        for e in en_drug_targets:
            if t in e:
                fc = fc_en_drug_targets[en_drug_targets.index(e)]
                fc_in_data[(drug[0])].add(fc)
                en_targets_in_data[(drug[0])].add(t)
            if t in ligands_by_receptor.keys():
                rc_targets_in_data[(drug[0])].add(t)


with open(os.path.join(OUTPUT, 'enzyme_product_dict.pkl'), 'rb') as fh:
    enzyme_product_dict = pickle.load(fh)
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