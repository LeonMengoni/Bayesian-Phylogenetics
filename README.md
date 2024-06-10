# Bayesian inference on phylogenetic trees

## Jukes-Cantor base substitution model

$$
p_\boldsymbol{\rho} = \left(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}\right)
$$

Define $\alpha$ as the rate of base substitution per unit time: $\alpha = \frac{\textit{number of base substitutions}}{\textit{unit time}}$

$$
Q = 
\begin{pmatrix}
-\alpha & \alpha/3 & \alpha/3 & \alpha/3 \\
\alpha/3 & -\alpha & \alpha/3 & \alpha/3 \\
\alpha/3 & \alpha/3 & -\alpha & \alpha/3 \\
\alpha/3 & \alpha/3 & \alpha/3 & -\alpha \\
\end{pmatrix}
$$

Markov matrix:

$$
M = e^{Qt} = 
\begin{pmatrix}
1-a & a/3 & a/3 & a/3 \\
a/3 & 1-a & a/3 & a/3 \\
a/3 & a/3 & 1-a & a/3 \\
a/3 & a/3 & a/3 & 1-a \\
\end{pmatrix}
$$

where $a = a(t) = \frac{3}{4}\left(1-e^{-\frac{4}{3}\alpha t}\right)$ is the probability that a base mutates to a different base over time $t$.

## Sampling parameter space

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


**NOTE**
Branch lengths sampled by the RevBayes algorithm:
- terminal branches are typically longer
- internal branches are shorter
