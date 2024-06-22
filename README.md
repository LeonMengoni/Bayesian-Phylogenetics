<h1 align="center">Information Theory and Inference<br> University of Padua <br>2023/2024</h1>

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

We will now explore these steps and present our code.

## 1. Model selection: Jukes-Cantor base substitution model

In order to infer a phylogeny, the first thing to do is to specify a model of base substitution. The Jukes-Cantor model (1969) is the simplest one: it assumes a uniform distribution of bases at the root node, indexed by $\boldsymbol{\rho}$, and defines a per-site mutation rate $\alpha=\frac{\textit{number of base substitutions}}{\textit{site}\cdot\textit{unit time}}$ that is equal for all base mutations. The alphabet that we will be using is $\lbrace A, G, C, T \rbrace$. The transition rate matrix $Q$ is therefore defined so that the rows sum up to 1. The element $Q(1,2) = q_{AG} = \frac{\alpha}{3}$ is, for example, the rate at which the substitution $A\rightarrow G$ occurs. 

$$
Q = \begin{pmatrix}
-\alpha & \alpha/3 & \alpha/3 & \alpha/3 \\
\alpha/3 & -\alpha & \alpha/3 & \alpha/3 \\
\alpha/3 & \alpha/3 & -\alpha & \alpha/3 \\
\alpha/3 & \alpha/3 & \alpha/3 & -\alpha \\
\end{pmatrix}
$$

In continuous time, the differential equation $\frac{d}{dt}p_t = p_t Q$ yields the solution $p(t) = p_0 M(t)$ where $M$ is a **Markov matrix**, in that its rows are normalized to 1, and $a$ is the probability of a base mutation over time $t$.

$$
M(t) = e^{Qt} = 
\begin{pmatrix}
1-a & a/3 & a/3 & a/3 \\
a/3 & 1-a & a/3 & a/3 \\
a/3 & a/3 & 1-a & a/3 \\
a/3 & a/3 & a/3 & 1-a \\
\end{pmatrix}
\qquad 
a = a(t) = \frac{3}{4}\left(1-e^{-\frac{4}{3}\alpha t}\right)
$$

At this point, we assume a **molecular clock**, i.e. a constant mutation rate $\alpha$ over evolutionary times. This allows us to define the _branch length_ $d = \alpha t$, which is the average number of base substitutions per site over elapsed time $t$. This number can be greater than 1 since it admits the possibility of multiple hidden substitutions at a site. Our model parameters at this stage are the tree topology $T$, the branch lengths $\lbrace t_e\rbrace_{e \in E(T)}$ and the mutation rate $\alpha$. In our code, we will assume $\alpha = 1$, so that $d = t$, and we reduce parameter space by one dimension: $\boldsymbol{\theta} = (T, \lbrace t\rbrace)$.

## 2. Data preparation: artificial or real data?

Generally, DNA sequences have to be aligned so that differing regions can be compared using the previously chosen model of nucleotide substitution. For simplicity, we will skip this step, working directly with aligned data. 

In order to evaluate the validity of our model and justify the choices made in setting up the MCMC algorithm, we will generate artificial data, starting from a tree with known topology and branch lengths. The code for generating an artificial tree is given by the function `simulate_artificial_data` in `data_generator.py`. 

Then, we will also use DNA sequences provided by RevBayes, to test the RevBayes scripts. 

## 3. Prior specification

This is one of the key steps in setting up a Bayesian inference. The prior distribution specifies our previous knowledge on the model parameters, and different distributions can lead to quite different posteriors. The tree topology and the branch lengths are independent, so we can write the prior as:

$$
P(\boldsymbol{\theta}) = P(T, \lbrace t\rbrace) = P(T)\cdot P(\lbrace t\rbrace)
$$

In literature [[1]](#1), the topology prior distribution is taken to be uniform, while for branch lengths, an exponential or uniform distribution seems to do the job. We will try both of them. If we choose an exponential distribution, we get:

$$
P(\lbrace t\rbrace) = \prod_{i=1}^E P(t_i) = \prod_{i=1}^E \lambda e^{-\lambda t_i} \propto \exp\left({-\lambda\sum_{i=1}^E t_i}\right) \propto e^{-\lambda L}
$$

where $L$ is the total length of the tree, defined as the sum of all the branch lengths (and $E$ is of course the total number of branches in the tree). If the individual branch lengths are distributed exponentially, the distribution of their sum will be a gamma distribution. Therefore, a different approach could be to sample the total tree length $L$ from a gamma distribution and then obtain the individual branch lengths proportionally with respect to a Dirichlet distribution. (**TODO**)

## 4. Likelihood calculation: Felsenstein's pruning algorithm

We can write the likelihood function as:

$$
\mathcal{L}(\boldsymbol{\theta}) = P(D|\boldsymbol{\theta}) = \prod_{s=1}^{S} P(D_s|\boldsymbol{\theta}) = \prod_{s=1}^{S} P(\lbrace B_i^s\rbrace|\boldsymbol{\theta})
$$

where $D_s = \lbrace B_i^s\rbrace$ indicates the sequence of bases $B_i$ (for all species considered) present at site $s\in\lbrace 1, 2,\cdots S\rbrace$, where $S$ is the sequence length. Taking as an example part of the data shown previously

```
Saimiri_sciureus               ATGACCTCAC...
Callicebus_donacophilus        ATGACCATTA...
Cebus_albifrons                ATGACCTCTC...
```

we would write $P(D_1|\boldsymbol{\theta}) = P(AAA|\boldsymbol{\theta})$ or $P(D_7|\boldsymbol{\theta}) = P(TAT|\boldsymbol{\theta})$. The individual terms are calculated using Felsenstein's pruning algorithm [[1]](#1).

## 5. MCMC simulation: traversing tree space

As we said previously, having set $\alpha = 1$, the parameters of our model are the tree topology and the branch lengths: $\boldsymbol{\theta} = (T, \lbrace t_e\rbrace)$. Our data is simply a list of genetic sequences: $D = (AGCTGA, GGCTGA, \dots , AGCTCC)$. The posterior is:

$$
P(\boldsymbol{\theta}|D) = \frac{P(D|\boldsymbol{\theta})\cdot P(\boldsymbol{\theta})}{P(D)} \propto \mathcal{L}(\boldsymbol{\theta})\cdot P(\boldsymbol{\theta})
$$

where $\mathcal{L}(\boldsymbol{\theta})$ is the Felsenstein likelihood and $P(\boldsymbol{\theta})$ is the prior.

To sample the posterior we use the **Metropolis-Hastings algorithm**; therefore, given $\pi(\boldsymbol{\theta}) = P(\boldsymbol{\theta}|D)$ as the stationary distribution of the Markov Chain, the Metropolis-Hastings acceptance ratio is:

$$
r = \min\left\lbrace 1, \frac{\pi(\boldsymbol{\theta'})}{\pi(\boldsymbol{\theta})} \frac{Q(\boldsymbol{\theta}|\boldsymbol{\theta'})}{Q(\boldsymbol{\theta'}|\boldsymbol{\theta})}\right\rbrace 
= \min\left\lbrace 1, \frac{\mathcal{L}(\boldsymbol{\theta'})}{\mathcal{L}(\boldsymbol{\theta})}\frac{P(\boldsymbol{\theta'})}{P(\boldsymbol{\theta})}\frac{Q(\boldsymbol{\theta}|\boldsymbol{\theta'})}{Q(\boldsymbol{\theta'}|\boldsymbol{\theta})}\right\rbrace
= \min\left\lbrace 1, \frac{\mathcal{L}(\boldsymbol{\theta'})}{\mathcal{L}(\boldsymbol{\theta})}\frac{P(\lbrace t_e'\rbrace)}{P(\lbrace t_e\rbrace)}\frac{Q(\boldsymbol{\theta}|\boldsymbol{\theta'})}{Q(\boldsymbol{\theta'}|\boldsymbol{\theta})}\right\rbrace
$$

where we have simplified $P(T)$ out of the equation since it is uniform. 

- If the prior distribution on the branch lengths is **uniform**, the acceptance ratio simplifies to:

$$
r = \min\left\lbrace 1, \frac{\mathcal{L}(\boldsymbol{\theta'})}{\mathcal{L}(\boldsymbol{\theta})}\frac{Q(\boldsymbol{\theta}|\boldsymbol{\theta'})}{Q(\boldsymbol{\theta'}|\boldsymbol{\theta})}\right\rbrace
$$

- If we assume the branch lengths to follow an **exponential** distribution, we can introduce the total tree length $L$ again, and calculate the acceptance ratio as:

$$
r = \min\left\lbrace 1, \frac{\mathcal{L}(\boldsymbol{\theta'})}{\mathcal{L}(\boldsymbol{\theta})}e^{-\lambda (L'-L)}\frac{Q(\boldsymbol{\theta}|\boldsymbol{\theta'})}{Q(\boldsymbol{\theta'}|\boldsymbol{\theta})}\right\rbrace
$$

Every step in the Markov chain is determined by a double proposal: first, a new tree topology is proposed according to a Nearest Neighbor Interchange (NNI) or a Subtree Prune and Regraft (described in [[1]](#1) and [[2]](#2)); second, a new set of branches is proposed according to $Q(\lbrace t'\rbrace|\lbrace t\rbrace)$, which is described by a so-called _Scaling move_. In practice, a scaling factor is generated uniformly from $\mathcal{U}\left[\frac{1}{\lambda},\lambda\right]$, which multiplies the current branch length to obtain the new branch length, as shown in the figure below [[4]](#4)).

<div align="center">
  <img src=Images/Scaling_move.jpg/ width=250 height=160>
</div>

## 6. Convergence diagnosis

Once the Markov Chain is set up and running, we have to run it long enough until convergence, to obtain samples from the posterior distribution. This is one of the main issues with MCMC, in general. In our case, the sheer size of the parameter space makes it practically impossible to be confident of having reached convergence. One strategy is to run multiple chains from different starting points and compare the results, since these should converge to the same distribution [[4]](#4).

## References

<a id="1">[1]</a>
Joseph Felsenstein. _Inferring Phylogenies_. Sinauer Associates, Sunderland, MA, 2004.

<a id="2">[2]</a>
Elizabeth S. Allman and John A. Rhodes. _Mathematical Models in
Biology: An Introduction_. Cambridge University Press, Cambridge,
2004.

<a id="3">[3]</a> 
Sebastian Höhna, Michael J. Landis, Tracy A. Heath, Bastien Boussau, Nicolas Lartillot, Brian R. Moore, John P. Huelsenbeck, Fredrik Ronquist. "RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language." Systematic Biology, 65:726-736, 2016.

<a id="4">[4]</a>
Barido-Sottani J, Schwery O, Warnock RCM et al. Practical guidelines for Bayesian phylogenetic inference using Markov Chain Monte Carlo (MCMC) [version 1; peer review: 3 approved, 1 approved with reservations]. Open Res Europe 2023, 3:204 (https://doi.org/10.12688/openreseurope.16679.1)

