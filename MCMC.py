from TreeClass import Phylogenetic_Tree
import copy
from Bio import Phylo


def initialize_random_tree(taxa_names, taxa_sequences, outgroup, distribution_params):
    tree_init = Phylogenetic_Tree().generate_random_topology(taxa_names=taxa_names) # Topology initialization
    tree_init.custom_root_with_outgroup(outgroup) # root with specified outgroup
    tree_init.generate_random_branch_lengths(distribution_params=distribution_params) # Branch lengths initialization
    
    # Clamp sequences to leaves
    for i, clade in enumerate(tree_init.get_terminals()):
        clade.sequence = taxa_sequences[i]

    return tree_init

class Sampler():
    def __init__(self, taxa_names, taxa_sequences, outgroup, distribution_params, nwalkers=1):
        self.initial_tree_list = [initialize_random_tree(taxa_names, taxa_sequences, outgroup, distribution_params) for i in range(nwalkers)]
        self.nwalkers = nwalkers
        self.distribution_params = distribution_params

    def run_mcmc(self, nsteps=1000, burn_in=100, save_frequency=100): # save = True
        self.nsteps = nsteps
        self.burn_in = burn_in
        self.save_frequency = save_frequency

        self.saved_trees = []
        self.accept_ratio = []
        
        for tree_init in self.initial_tree_list:
            saved_trees_per_walker = []
            accept = 0
            tree_init_copy = copy.deepcopy(tree_init)
            for n in range(self.nsteps):
                tree_init_copy = tree_init_copy.NNI_step(self.distribution_params)
                if tree_init_copy.accepted: accept +=1
                if n > self.burn_in and n % self.save_frequency == 0:
                    saved_trees_per_walker.append(tree_init_copy)
            self.saved_trees.append(saved_trees_per_walker)
            self.accept_ratio.append(accept / self.nsteps)

    def save_to_file(self, file_name):
        flattened_saved_trees = [tree for saved_trees_per_walker in self.saved_trees for tree in saved_trees_per_walker]
        Phylo.write(flattened_saved_trees, file_name, "nexus")
