# Bayesian inference on phylogenetic trees

## Sampling parameter space

The parameters of our Jukes-Cantor model are: 
$$\boldsymbol{\theta} = (T, \lbrace t_e\rbrace_{e \in E(T)}, \alpha)$$

Our data is composed of a list of genetic sequences:
$$D = (AGCTGA, GGCTGA, \dots , AGCTCC)$$

The posterior can be written as:
$$P(\boldsymbol{\theta}|D) = P(T, \lbrace t_e\rbrace, \alpha|D) = P(T|\lbrace t_e\rbrace, \alpha, D) \cdot P(\lbrace t_e\rbrace|D) \cdot P(\alpha|D)$$
