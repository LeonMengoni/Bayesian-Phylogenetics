from settings import *
import numpy as np
from keras.utils import to_categorical
import utils as ut
import warnings
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from matplotlib import pyplot as plt
import copy
import random

# Root: MRCA (Most Recent Common Ancestor) of all known given present known taxa ("terminal" taxa of "leaves")
# Parent: "parent" node, given an edge
# Child/children: "children" nodes or "child" node, given an edge 
# Ancestor: any node along the chain of parents, starting from a node 
# Descendant: any node belonging to a clade, starting from a node
# Terminal: "leaf" nodes in the tree

class Phylogenetic_Tree(Phylo.BaseTree.Tree):
    def __init__(self, root=None):  
        super().__init__(root=root)
        self.root.name = "Root"
        self.parents = self.dict_parents()
        self.n_taxa = len(self.get_terminals())

    @classmethod
    def read_from_file(cls, file_path, file_format):
        tree = Phylo.read(file_path, file_format)
        return cls(root=tree.clade)

    @classmethod
    def generate_random_topology(cls, taxa_names=None, n_taxa=None): # 
        if taxa_names is None: # generate tree from a random number of randomly named taxa
            if n_taxa is None: n_taxa = np.random.randint(3,100)
            taxa_names = [f"S{i}" for i in range(n_taxa)]
        
        else: # generate tree from predefined terminal taxa
            if n_taxa: warnings.warn(f"Parameter n_taxa = {n_taxa} is redundant")
            n_taxa = len(taxa_names)

        tree = Tree().randomized(taxa_names, branch_length=None)
        instance = cls(root=tree.clade)
        return instance
    
    def dict_parents(self): # dictionary with child-parent relationships 
        parents = {child:parent for parent in self.find_clades(order="level") for child in parent.clades}
        return parents
    
    def get_sibling(self, node):
        sibling = [clade for clade in self.parents[node].clades if clade != node][0]
        return sibling

    def sample_branch_length(self, params):
        if params["distribution"] == "exponential":
            return np.random.exponential(scale=params["scale"])
        elif params["distribution"] == "uniform":
            return np.random.uniform(low=params["low"], high=params["high"])
        else:
            raise ValueError("Unsupported distribution type")
        
    def calculate_total_tree_length(self):
        L = np.sum([clade.branch_length for clade in self.find_clades() if clade != self.root])
        return L
    
    def generate_random_branch_lengths(self, distribution_params, type=None):
        for clade in self.find_clades():
            if clade != self.root:
                clade.branch_length = self.sample_branch_length(distribution_params)
                # if distribution == "exponential": # Exponential distribution (scale is mean)
                #     clade.branch_length = np.random.exponential(scale=0.1)
                # elif distribution == "uniform": # Uniform distribution
                #     clade.branch_length = np.random.random()
        
        if type == "ultrametric": # Ultrametric tree: distance to root is the same for all leaves 
            self.make_tree_ultrametric()

    def make_tree_ultrametric(self):
        distances_to_root = [(clade.name, self.distance(clade)) for clade in self.get_terminals()]
        present = max(distances_to_root, key=lambda x: x[1])[1]

        for terminal_clade in self.get_terminals():
            if terminal_clade == self.root: continue
            parent =  self.parents[terminal_clade]
            terminal_clade.branch_length = present - self.distance(parent)

    def generate_sequences(self, seq_len): # generate root sequence and descendant sequences according to Jukes-Cantor substitution model
        self.seq_len = seq_len
        for parent in self.find_clades(): 
            if parent == self.root:
                parent.sequence = ut.generate_sequence(self.seq_len)
            for child in parent.clades:
                child.sequence = ut.generate_child_from_parent(parent.sequence, child.branch_length)
                child.mutations = ut.count_mutations(parent.sequence, child.sequence)

    def custom_root_with_outgroup(self, clade_name):
        self.root_with_outgroup(clade_name)
        self.parents = self.dict_parents() # Update parents dictionary

    # Felsenstein pruning algorithm for calculating likelihood
    def calculate_likelihood(self, p_root=np.ones(n_bases)/n_bases, log=False):
        for parent in self.get_nonterminals(order="postorder"):
            for child in parent.clades:
                if child.is_terminal():
                    child.C = to_categorical(ut.sequence_to_index(child.sequence), num_classes=n_bases)
                child.W = ut.Markov_matrix(child.branch_length) @ child.C.T
            parent.C = (np.multiply(parent.clades[0].W, parent.clades[1].W)).T

        joint_prob = p_root @ self.root.C.T
        likelihood = np.prod(joint_prob)
        if log:
            log_likelihood = np.sum(np.log(joint_prob))
            return likelihood, log_likelihood
        else: 
            return likelihood
        
    # Nearest Neighbor Interchange (NNI) methods
    def NNI_eligible_nodes(self):
        i = 0
        NNI_nodes = []
        for node in self.get_nonterminals(order="postorder"):
            if node == self.root: continue
            elif self.parents[node] == self.root:
                i += 1
                if i == 2: # If root node has two children within the nonterminals, then I skip the first one and get only the second one. 
                    NNI_nodes.append(node)
            else:
                NNI_nodes.append(node)

        assert self.n_taxa-3 == len(NNI_nodes) 
        return NNI_nodes
    
    def NNI_permutable_nodes(tree, chosen_node): # 3 nodes that can be permuted in a rooted tree NNI step
        permutable = []
        permutable += chosen_node.clades # 2 children
        sibling = tree.get_sibling(chosen_node)
        if tree.parents[chosen_node] == tree.root: # if parent is root, get one nibling
            niblings = sibling.clades
            permutable.append(niblings[0])
        else: # if parent is not root, get sibling
            permutable.append(sibling)
        assert len(permutable) == 3 
        return permutable

    def NNI_generate_topology(tree, chosen_node, permutable):
        NNI_tree = copy.deepcopy(tree)
        for NNI_clade in NNI_tree.get_nonterminals():
            if NNI_clade.name == chosen_node.name:
                NNI_clade.clades = permutable[:2]
                NNI_parent = NNI_tree.parents[NNI_clade]
                if NNI_parent == NNI_tree.root: # if parent is root, get the other nibling (the one that wasn't placed in permutable)
                    NNI_sibling = NNI_tree.get_sibling(NNI_clade)
                    NNI_niblings = NNI_sibling.clades
                    NNI_sibling.clades = [permutable[2], NNI_niblings[1]]
                else: # if parent is not root
                    NNI_parent.clades = [NNI_clade, permutable[2]]
        NNI_tree.parents = NNI_tree.dict_parents() # Update parents dictionary (only structure whose elements have undergone transformation when copying tree)
        return NNI_tree
    
    def NNI_step(self, distribution_params):
        NNI_nodes = self.NNI_eligible_nodes()
        chosen_node = random.choice(NNI_nodes)
        permutable = self.NNI_permutable_nodes(chosen_node)
        _, self.loglikelihood = self.calculate_likelihood(log=True)

        for i in range(2):
            new_nodes = permutable[-(i+1):] + permutable[:-(i+1)]
            NNI_tree = self.NNI_generate_topology(chosen_node, new_nodes)
            NNI_tree.generate_random_branch_lengths(distribution_params=distribution_params)
            _, NNI_tree.loglikelihood = NNI_tree.calculate_likelihood(log=True)

            delta_log_likelihood = NNI_tree.loglikelihood - self.loglikelihood
            if distribution_params["distribution"] == "exponential": 
                if not hasattr(self, 'L'):
                    self.L = self.calculate_total_tree_length() # total tree length initial tree
                NNI_tree.L = NNI_tree.calculate_total_tree_length() # total tree length proposal tree
                delta_log_likelihood -= (NNI_tree.L - self.L) / distribution_params["scale"] # lambda of exponential distirbution, interval of uniform distribution
            elif distribution_params["distribution"] == "uniform":
                pass
            rho = np.min([1, np.exp(delta_log_likelihood)])

            if np.random.random() < rho:
                NNI_tree.accepted = True
                return NNI_tree
        
        self.accepted = False
        return self
    
    # Visualization
    def draw(self, figsize=(4,4)):
        _, ax = plt.subplots(figsize=figsize)
        Phylo.draw(self, axes=ax)