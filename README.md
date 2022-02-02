This is the github folder for the paper published in Neural Computation **Confidence-Controlled Hebbian Learning Efficiently Extracts Category Membership From Stimuli Encoded in View of a Categorization Task** in 2021 (https://direct.mit.edu/neco/article-abstract/34/1/45/107909/Confidence-Controlled-Hebbian-Learning-Efficiently)
For open-access: https://www.biorxiv.org/content/10.1101/2020.08.06.239533v2.article-info

The main part of the code lies within the *src folder*. 

The other folders are a mix of exploratory analysis and generation of the figures. Thus they are not necessarily polished and highly commented. 

Don't hesitate to contact me if you have any questions about this project at kevin.berlemont@gmail.com. 

**src folder**

This folder has been cleaned up to help for comprehension, which is not necessarily the case for all the other files of the project.

- *dynamical_system.jl*: this file contains all the Julia function to generate the decision dynamics of the attractor neural network with the coding and decision layer
- *information_auto.jl*: This file contains all the necessary funcitons to compute the information theory related quantities using automatic differentiation
- *information_theory.jl*: same as the previous file but with all the quantities written explicitily within the functions
- *load_data.jl*: functions that are used to load specific data for the paper figures
- *optimal_tuning_curves.jl*: contains the functions to compute the optimal tuning curves of a neural coding layer to perform a categorization task (it follows the formulas in the paper)
- *structures.jl*: contains the structures that are used to perform the numerical simulations
- *training_network.jl*: contains all the learning algorithms for the two-layers neural network

**scripts folder**

This folder contains various scripts that will produce the data necessary for the paper and the project. Please note that some part will have to be uncommented to run through all the parameters. This happens when I was running new parameters and didn't wanted to run an extensive simulations through parameters already studied.

**notebooks folder**

This fodler contains the various notebooks I used through the project. Some of them are just exploratory notebooks, other contains the code to reproduce the figures of the paper. For instance **notebook 12 and 15** will reproduce the main figures of the paper. I left the collection of all my notebooks for informational purpose.
