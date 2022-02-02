# 15-Fcode through learning

## Model

The goal of this analysis is to study in more details the behavior of the network during learning. More specifically the impact of the neural code and/or the learning algorithm on the learning process.

I will make use of the quantity $F_{code}$ and I recall that:
$$
\begin{equation}
	I(\mu,x) - I(\mu,r) = \frac12 \int dx \, p(x) \frac{F_{\text{cat}}(x)}{F_{\text{code}}(x)}
\end{equation}
$$


This quantity is computed by using bell-shaped tuning curves for the coding neurons with specific gains. However the gains need to be positive numbers. After the learning half of the weights are negative. I will assume that the gain of the tuning curve correspond to the absolute value of the synaptic weight after learning. This assumption allows to compute $F_{code}$ and depends on the decaying rate of the tuning curves.

## Numerical simulations

As already observed, the network is able to rapidly learn  the task. I will first detail what happens during the first 1000 learning steps.

### $\alpha = 0.25$

![First-steps-N15-Alpha25](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\First-steps-N15-Alpha25.png)

The curves in red correspond to a coding layer that is uniform and in blue an optimized one. The dashed lines stand for a learning without confidence and the other lines for a learning with confidence.

Each color stands for a specific number of trials during learning.

One can see that confidence has little impact when the coding is uniform. The Fisher information ism ore or less constant during learning. However, when the coding is optimized, confidence allows a strong increase of Fisher information at the boundary between categories.

![First-steps-N35-Alpha25](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\First-steps-N35-Alpha25.png)

If the number of neurons $N$ of the coding layer increases, the behavior is similar. Yet is more noisy as the tuning curves are sharper.

### $\alpha = 0.4$

![First-steps-N35-Alpha4](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\First-steps-N35-Alpha4.png)![First-steps-N15-Alpha4](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\First-steps-N15-Alpha4.png)

When the overlap between categories increases one can note that the difference between the two types of coding disappear.

### Longer learning time

The following figures are the same as above but for longer learning time. 

![Long-steps-N15-Alpha25](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\Long-steps-N15-Alpha25.png)



![Long-steps-N35-Alpha25](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\Long-steps-N35-Alpha25.png)





![Long-steps-N15-Alpha4](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\Long-steps-N15-Alpha4.png)

![Long-steps-N35-Alpha4](C:\Documents\1 - PhD Projects\4-project_coding_categories\plots\notebooks\notebook15\Long-steps-N35-Alpha4.png)