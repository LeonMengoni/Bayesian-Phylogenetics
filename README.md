<h1 align="center">Information theory and Inference<br> University of Padua <br>2023/2024</h1>

<h3 align="center"><b>Group members:</b> Leon Mengoni, Alessio Saccomani</h3>

# Bayesian inference of phylogenetic trees

<div align="center">
  <img src=Images/Dinosaur_meme.jpg/ width=400 height=372>
</div>

**Phylogenetics** is the study of evolutionary relationships among biological entities – often species, individuals, or genes. These relationships are represented by a **phylogenetic tree** (or _phylogeny_), a diagram that represents the evolutionary history and common ancestry of the entities being studied. The branching pattern of the tree indicates how these different entities are related, with each branch point, or node, representing a common ancestor. The species whose phylogeny we want to reconstruct are placed at the **leaves**, the terminal nodes of the tree.

The species that we want to relate can be characterized by a series of observable heritable traits. We will be using **aligned DNA sequences**, but also morphological traits or protein amino acid sequences are commonly used. Each DNA sequence in our dataset corresponds to one species. For example, our data can look like these primate DNA sequences: 

```
Saimiri_sciureus               ATGACCTCAC...
Callicebus_donacophilus        ATGACCATTA...
Cebus_albifrons                ATGACCTCTC...
Aotus_trivirgatus              ATGACTTCTC...
Hylobates_lar                  ATGACCCCCC...
Pan_paniscus                   ATGACCCCAA...
```

We assume that these species are descendants of a common ancestor, and that their DNA sequences have evolved from their ancestors' DNA by nucleotide substitution. **Nucleotide substitution** refers to the process by which one nucleotide (A, G, C or T) in a DNA sequence is replaced by another over time. These substitutions occur due to various evolutionary forces such as mutations, natural selection, genetic drift, and recombination. To accurately infer evolutionary relationships, scientists use models of nucleotide substitution. These models describe the rates at which one nucleotide is expected to be replaced by another over evolutionary time. Several models exist, ranging from simple to complex (Jukes-Cantor Model, Kimura Two-Parameter Model, General Time Reversible (GTR) Model): we will use the **Jukes-Cantor** model, that assumes that all nucleotide substitutions occur at the same rate.

An example of a tree relating the above primates could be the following:

<div align="center">
  <img src=Images/Example_tree.jpg/ width=800 height=445>
</div>

A common approach for inferring this sort of tree is **Bayesian phylogenetics**. This approach applies Bayesian statistics to infer the most probable phylogenetic tree given the observed data and prior information. Bayesian methods use _Markov Chain Monte Carlo_ (MCMC) algorithms to sample from the posterior distribution, generating a set of trees and parameters that are consistent with the observed data and the prior information. The process of Bayesian phylogenetic analysis typically involves several key steps:

1. **Model selection**: choosing the appropriate model of nucleotide substitution.
2. **Data preparation**: aligning sequences and preparing the data matrix.
3. **Prior specification**: defining the prior distributions for the parameters.
4. **Likelihood calculation**: calculating the likelihood based on the model. 
5. **MCMC simulation**: running the MCMC algorithm to sample from the posterior distribution.
6. **Convergence diagnosis**: checking that the MCMC simulation has adequately explored the parameter space.
7. **Tree summarization**: summarizing the sampled trees into a consensus tree, which represents the most probable phylogenetic relationships.

We will now explore these steps.

## 1. Model selection: Jukes-Cantor base substitution model

In order to infer a phylogeny, the first thing to do is to specify a model of base substitution. The Jukes-Cantor model (1969) is the simplest one: it assumes a uniform distribution of bases at the root node, indexed by $\boldsymbol{\rho}$, and defines a per-site mutation rate $\alpha=\frac{\textit{number of base substitutions}}{\textit{site}\cdot\textit{unit time}}$ that is equal for all base mutations. The alphabet that we will be using is $\lbrace A, G, C, T \rbrace$.

$$ p_\boldsymbol{\rho} = \left(p_A=\frac{1}{4}, p_G=\frac{1}{4}, p_C=\frac{1}{4}, p_T=\frac{1}{4}\right) $$

The transition matrix $Q$ is therefore defined so that the rows sum up to 1. The element $Q(1,2) = q_{AG} = \frac{\alpha}{3}$ is, for example, the rate at which the substitution $A\rightarrow G$ occurs. 

$$
Q = \begin{pmatrix}
-\alpha & \alpha/3 & \alpha/3 & \alpha/3 \\
\alpha/3 & -\alpha & \alpha/3 & \alpha/3 \\
\alpha/3 & \alpha/3 & -\alpha & \alpha/3 \\
\alpha/3 & \alpha/3 & \alpha/3 & -\alpha \\
\end{pmatrix}
$$

In continuous time, the differential equation that we want to solve is:

$$
\frac{d}{dt}p_t = p_t Q
$$

which yields the solution 

$$
M(t) = e^{Qt} = 
\begin{pmatrix}
1-a & a/3 & a/3 & a/3 \\
a/3 & 1-a & a/3 & a/3 \\
a/3 & a/3 & 1-a & a/3 \\
a/3 & a/3 & a/3 & 1-a \\
\end{pmatrix}
\qquad a = a(t) = \frac{3}{4}\left(1-e^{-\frac{4}{3}\alpha t}\right)
$$

The parameter $a$ is interpreted as the probability of a base mutation over time $t$. The matrix $M$ is a **Markov matrix**, in that its rows are normalized to 1, so that $p(t) = p_0 M(t)$.

At this point, we assume a **molecular clock**, i.e. a constant mutation rate $\alpha$ over evolutionary times. This allows us to define the _branch length_ $d = \alpha t$, which is the average number of base substitutions per site over elapsed time $t$. This number can be greater than 1 since it admits the possibility of multiple hidden substitutions at a site. Our model parameters at this stage are $\boldsymbol{\theta} = (T, \lbrace t_e\rbrace_{e \in E(T)}, \alpha)$, i.e. the tree topology $T$, the branch lengths $\lbrace t_e\rbrace_{e \in E(T)}$ and the mutation rate $\alpha$. In our code, we will assume $\alpha = 1$, so that $d = t$, and we reduce parameter space by one dimension. 

## 2. Data preparation: artificial or real data? (TODO)

We will be conducting a double analysis, one with "real" data, provided by RevBayes 


## 3. Prior specification

## 4. Likelihood calculation

## 5. MCMC simulation: traversing tree space

## 6. Convergence diagnosis

## 7. Tree summarization


## Bayesian inference TODO

The parameters of our Jukes-Cantor model are: 
$$\boldsymbol{\theta} = (T, \lbrace t_e\rbrace_{e \in E(T)}, \alpha)$$

Our data is composed of a list of genetic sequences:
$$D = (AGCTGA, GGCTGA, \dots , AGCTCC)$$

The posterior can be written as:
$$P(\boldsymbol{\theta}|D) = P(T, \lbrace t_e\rbrace, \alpha|D) = P(T|\lbrace t_e\rbrace, \alpha, D) \cdot P(\lbrace t_e\rbrace|D) \cdot P(\alpha|D)$$

The first term can be expanded as:
$$P(T|\lbrace t_e\rbrace, \alpha, D) = \frac{P(D|T, \lbrace t_e\rbrace, \alpha) \cdot P(T|\lbrace t_e\rbrace, \alpha)}{P(D|\lbrace t_e\rbrace, \alpha)} = \frac{P(D|T) \cdot P(T)}{P(D)} \propto \mathcal{L}(T)$$

where $P(D|\lbrace t_e\rbrace, \alpha) = P(D)$ and $P(T|\lbrace t_e\rbrace, \alpha) = P(T)$ hold for the marginal likelihood and the prior probability distribution for the tree, respectively, and where $P(D|T)=\mathcal{L}(T)$ is the likelihood function for the selected tree $T$.

The second and third terms can be simplified as:
$$P(\lbrace t_e\rbrace|D) = \frac{P(D|\lbrace t_e\rbrace) \cdot P(\lbrace t_e\rbrace)}{P(D)} = P(\lbrace t_e\rbrace)$$
$$P(\alpha|D) = \frac{P(D|\alpha) \cdot P(\alpha)}{P(D)} = P(\alpha)$$
since the data doesn't solely depend on the values of the edge lengths and of $\alpha$, but only when these are combined with a tree topology. In other words, the data doesn't update our beliefs on what should be the probability distribution for these parameters. Certainly, if the data is very variable, we could say that it is a consequence of the fact that the substitution rate $\alpha$ is big, or that the edges are especially long. Therefore, they are simply represented by their prior probability distributions, which we may consider uniform on some range of values. 

In summary,

$$P(\boldsymbol{\theta}|D) =\frac{\mathcal{L} \cdot P(T)}{P(D)} \cdot P(\lbrace t_e\rbrace) \cdot P(\alpha) \propto \mathcal{L}$$

In the framework of our MCMC sampling algorithm, the likelihood is the only function that is of interest for computing the Metropolis acceptance ratio.

## Felsenstein's pruning algorithm TODO

**TODO**
Add tree files that are used in the code 

**NOTE**
Branch lengths sampled by the RevBayes algorithm:
- terminal branches are typically longer
- internal branches are shorter

## References

<a id="1">[1]</a> 
Felsenstein, J. (2004). _Inferring phylogenies_. Sinauer Associates, Inc. 

<a id="2">[2]</a> 
Höhna, Landis, Heath, Boussau, Lartillot, Moore, Huelsenbeck, Ronquist. 2016. "RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language." Systematic Biology, 65:726-736.

