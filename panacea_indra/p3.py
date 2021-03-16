from functions import *


# Di-graph of interactions between top 10 DE ligands, enzyme products
# and receptors

with open(os.path.join(OUTPUT, 'enzyme_product_dict.pkl'), 'rb') as fh:
    enzyme_product_dict = pickle.load(fh)

with open(os.path.join(OUTPUT, 'products_receptors.pkl'), 'rb') as fh:
    products_receptors = pickle.load(fh)
    
with open(os.path.join(OUTPUT, 'ligands_FC.pkl'), 'rb') as fh:
    ligands_FC = pickle.load(fh)

with open(os.path.join(OUTPUT, 'enzymes_FC.pkl'), 'rb') as fh:
    enzymes_FC = pickle.load(fh)

with open(os.path.join(OUTPUT, 'stmts_db_by_cell_type.pkl'), 'rb') as fh:
    stmts_db_by_cell_type = pickle.load(fh)

with open(os.path.join(OUTPUT, "receptors.csv"), 'r') as fh:
    receptors_in_data = {l.strip() for l in fh}


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
                               os.path.join(cell_type,cell_type + "_"),
                               enzyme_product_dict,
                               products_receptors)
    
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
create_interaction_digraph(ligandsFC_by_receptor, 
                           sorted_enzyme_FC,
                           '', 
                           enzyme_product_dict,
                           products_receptors)