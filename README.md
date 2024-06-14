# Bayesian inference of phylogenetic trees

**Phylogenetic trees**, or phylogenies, are structures that describe the ancestral and evolutionary relationships between species, starting from a series of observable heritable traits. The type of observable traits that we will be using are **aligned data sequences**, but also morphological traits or protein amino acid sequences are commonly used.

Each DNA sequence in our dataset corresponds to one species. For example, our data can look like these primate DNA sequences: 

```
Saimiri_sciureus               ATGACCTCAC...
Callicebus_donacophilus        ATGACCATTA...
Cebus_albifrons                ATGACCTCTC...
Aotus_trivirgatus              ATGACTTCTC...
Hylobates_lar                  ATGACCCCCC...
Pan_paniscus                   ATGACCCCAA...
```

From this data, we want to infer a tree that is able to explain the observed base mutations. To do this, we will be adopting Bayesian methods of phylogenetic inference, based on _Markov Chain Monte Carlo_ (MCMC). 

<div align="center">
  <img src=Images/Example_tree_with_photos.jpg/>
</div>

## Jukes-Cantor base substitution model

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

The parameter $a$ is interpreted as the probability of a base mutation over time $t$. The matrix $M$ is a **Markov matrix**, in that its rows are normalized to 1, and it holds that

$$
p(t) = p_0 M(t)
$$

At this point, we assume a **molecular clock**, i.e. a constant mutation rate $\alpha$ over evolutionary times. This allows us to define the _branch length_ $d = \alpha t$, which is the average number of base substitutions per site over elapsed time $t$. This number can be greater than 1 since it admits the possibility of multiple hidden substitutions at a site.

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


**NOTE**
Branch lengths sampled by the RevBayes algorithm:
- terminal branches are typically longer
- internal branches are shorter
