# Bayesian inference on phylogenetic trees

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

In the framework of our MCMC sampling algorithm, the likelihood is the only function that is of interest for computing the Metropolis acceptance ratio.
