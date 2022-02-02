# 14-Reward modulated - Engel 2015

## Model presentation

The network is composed of the coding and decoding layer as usual. The goal is to compare the classic reward modulated learning and our confidence modulated learning. 

With the reward modulated learning, after each trial the synaptic weights modification is defined as:
$$
\Delta w_{ij} = \lambda R_{\theta} r_j
$$
with $R_{\theta}$ the reward modulation at stimulus $\theta$.
$$
R_{\theta} = R_{\theta} + (R-R_{\theta})/\tau
$$
with $R$ the reward at this trial.

In order to be able to compare with our learning algorithm we impose a normalization of the weights. Even if this reward modulation leads to a non-divergent learning algorithm.

## Without confidence

The first question one can ask is whether the reward modulation $R_{\theta}$ has an impact on the learning process. The following figures investigate the performances of the network and the ones where confidence does not play a role.

##### Uniform TC and $\alpha  = 0.2$



![fig1](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig1.svg)

##### Optimized TC and $\alpha  = 0.2$

![fig2](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig2.svg)

##### Uniform TC and $\alpha  = 0.4$

![fig5](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig5.svg)

##### Optimized TC and $\alpha  = 0.4$



![fig6](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig6.svg)

One can see that having a reward modulation of the learning rate leads to better performances of the network (or can't really conclude as some approximations have been made in the computations).





## With confidence modulation

##### Uniform TC and $\alpha  = 0.2$

![fig7](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig7.svg)

##### Optimized TC and $\alpha  = 0.2$



![fig8](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig8.svg)

##### Uniform TC and $\alpha  = 0.4$

![fig11](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig11.svg)

##### Optimized TC and $\alpha  = 0.4$

![fig12](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook14\fig12.svg)

When confidence modulates the learning, the performances are increased. Even when comparing to the fully reward modulated learning, the confidence modulations leads to better performances.