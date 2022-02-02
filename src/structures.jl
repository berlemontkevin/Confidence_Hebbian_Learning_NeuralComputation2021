
## Loading the different structures
using Parameters

##########
# Structures for the attractor network dynamics
###########
@with_kw struct MyAttractorNetwork
  # Definition du type "My Attractor Network"
  # chaque champ est un des paramètres du modèle
  noiseamp::Float64  =0.02# noise amp of the OU process
  i0::Float64 = 0.3255 #resting state of the OU process
  jn11::Float64 =0.2609# Synaptic strength unit 1 to 1
  jn21::Float64 = 0.0497# Synaptic strength unit 2 to 1
  jn12::Float64 = 0.0497# Synaptic strength unit 1 to 2
  jn22::Float64 = 0.2609# Synaptic strength unit 2 to 2
  tnmda::Float64 =100.0# time constant of NMDA receptors
  γ::Float64 = 0.000641
  tampa::Float64 = 2.0# time constant of AMPA receptors
  threshold::Float64# threhsold of the decision
  end


struct MyEulerParameters
  # Contient tous les paramètres nécessaires pour la simulation
  #Méthode d'intégration Euler
  dt::Float64
  timewind::Float64 # time windows for computing the mean of the variables
  slidewind::Float64 # sliding winows for computing the mean of variables
end

struct MyTransferFunction
  # Contient les détails de la fonction de transfer input-firing rate of the neurons
  # cette fonction vient de ... (TODO inclure ref.)
  a::Float64
  b::Float64
  d::Float64
end

mutable struct MyInput
  # define the parameters of the input sent into both units
  i1::Float64
  i2::Float64
  μ0::Float64
  jaext::Float64
end

mutable struct MyUnit
  # structure that defines the variables of a neural population
  rate::Float64
  s::Float64 # synaptic gatign variable
  isyn_coming::Float64
  inoise::Float64 # current noise going into the unit
  meanrate::Array{Float64} # current sliding rate for the mean
end


###############
# Coding layer structures
################


mutable struct MyCodingLayer
      centers::Array{Float64}
      widths::Array{Float64}
      weightsto1::Array{Float64}
      weightsto2::Array{Float64}
      gain::Array{Float64}
end



######################
# Structures of the different categories
##################

abstract type MyCat
end

@with_kw struct MyCatGauss <: MyCat
    centers::Array{Float64} = [0.0,1.0]
    width::Array{Float64}
    proba::Array{Float64} = [0.5,0.5]
    name::String = "Gaussian_Categories"
    end

@with_kw struct MyCatWeibull <: MyCat
    centers::Array{Float64} = [0.0,1.0]
    alpha::Float64
    beta::Float64
    width::Float64 # this width corresponds to the Gaussian equivalent
    name::String = "Weibull_Categories"

end


#################
# Structure for the simulations parameters
############

@with_kw struct ListParameters_θconfidenceSimulation
  ## Definition List Parameters
  # used in script 01-...
  # TODO customize save for ListParameters
  # TODO check output of the different functions
  #TODO check the type of folder we want to save
N::Int
T_simu::Int
ΔTsave::Int
Seuil_confidence::Float64
seuil_decision::Float64
type_TC::String
α::Float64
cats::MyCat
box_conv::Float64=1.0 # default value of box convolution
learning_rate::Float64 = 0.005 # default value
cst_normalisation::Float64 = 0.25 # this value is a test #0.5 #default value # Dans tous les cas cette valeur est obtenue en modiafiant jaext et en supposant un courant externe constant
#corrupted::Float64
# this value of 0.25 is because 100% acc at this level of perf
end
