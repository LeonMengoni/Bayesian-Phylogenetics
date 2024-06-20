from TreeClass import Phylogenetic_Tree

def initialize_random_tree(taxa_names, taxa_sequences, outgroup, distribution_params):
    tree_init = Phylogenetic_Tree().generate_random_topology(taxa_names=taxa_names) # Topology initialization
    tree_init.custom_root_with_outgroup(outgroup) # root with specified outgroup
    tree_init.generate_random_branch_lengths(distribution_params=distribution_params) # Branch lengths initialization
    
    # Clamp sequences to leaves
    for i, clade in enumerate(tree_init.get_terminals()):
        clade.sequence = taxa_sequences[i]

    return tree_init