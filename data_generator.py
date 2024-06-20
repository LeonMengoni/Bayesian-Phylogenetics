from TreeClass import Phylogenetic_Tree

def simulate_artificial_data(n_taxa, seq_len, outgroup, distribution_params):
    true_tree = Phylogenetic_Tree().generate_random_topology(n_taxa=n_taxa)
    true_tree.custom_root_with_outgroup(outgroup)
    true_tree.generate_random_branch_lengths(distribution_params=distribution_params)
    true_tree.generate_sequences(seq_len)
    data = {clade.name:clade.sequence for clade in true_tree.get_terminals()}
    return true_tree, data
