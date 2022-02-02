using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using Distributed
using BSON
using DataFrames

addprocs(9)
@everywhere using DrWatson
@everywhere quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
@everywhere using Distributed,BSON,DataFrames
@everywhere include(srcdir("training_network.jl"))



# In this file: various parameters are tested
# We save various data TODO

### Important
# We do not change the value of α for the Gaussian ctaegorie
# It stays constant through the learning
# We save the data for Uniform, Hebbian and Optimized tuning curves
# We use Poisson Input for the codign layer andn ot jsut the mean anymore
#########
############
#%% Default values of savenames
#########

@everywhere DrWatson.default_prefix(e::ListParameters_θconfidenceSimulation) = ""
# dossier de type
# sims/scriptname/results/Modele/TypeCat/Seuil
@everywhere DrWatson.allaccess(::ListParameters_θconfidenceSimulation) = (:N,
:T_simu,
:ΔTsave,
:Seuil_confidence,
:seuil_decision,
:type_TC,
:α,
:box_conv)

@everywhere function foldername_create(dummy::ListParameters_θconfidenceSimulation,box::Float64)
    seuil = dummy.seuil_decision
    Modele = "Decision_Seuil=$seuil"
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    f = "$Modele\\$tcatname\\"
    return f
end
#<<<<<<<<<<<<<<<<<<<<<<<<<<<

######
# Function to run the simulations
############
@everywhere function create_simu_from_dict(diclist)
	@unpack N,T_simu,ΔTsave,Seuil_confidence,seuil_decision,type_TC,α,box_conv = diclist

    mycatGauss = MyCatGauss(width=[α,α])
    simulations = ListParameters_θconfidenceSimulation(N=N,
    T_simu=T_simu,
    ΔTsave=ΔTsave,
    Seuil_confidence=Seuil_confidence,
    seuil_decision=seuil_decision,
    type_TC=type_TC,
    α = α,
    cats=mycatGauss,
	box_conv = box_conv)

    return simulations
end

@everywhere function mean_learning_Gaussian_Cat(diclist,box=1.0,rep=10)
    list_weights_to_unit1 = [] # chaque elt correspond à une simu
    list_weights_to_unit2 = [] # chaque elt correspond à une simu
	box = diclist["box_conv"]
    simulations = create_simu_from_dict(diclist)

    for i=1:rep
        Weightsto1,Weightsto2 = learning_poisson_input(simulations,box)
        push!(list_weights_to_unit1,Weightsto1)
        push!(list_weights_to_unit2,Weightsto2)
        Weightsto1 = nothing
        Weightsto2 = nothing
    end



    data = @dict list_weights_to_unit1 list_weights_to_unit2 simulations

    #save(datadir(fname,sname),data)
    return data
end


@everywhere function save_meanlearning_Gauss(diclist)
    box = diclist["box_conv"]

    dummy = mean_learning_Gaussian_Cat(diclist,box)
    @unpack list_weights_to_unit1,list_weights_to_unit2,simulations = dummy
    sname = savename(simulations,"bson")
    fnametemp = foldername_create(simulations,box)
    fname = "sims\\03-script\\results\\$fnametemp"
    fnameparam = "sims\\03-script\\param\\"

    mkpath(datadir(fname))
    mkpath(datadir(fnameparam))

    save(datadir(fname,sname),dummy)
    save(datadir(fnameparam,sname),diclist)
    return true
end


## We have a defualt learning rate of 0.005
# a default normlisation cst of 0.5

# initialize lsit of number of neurones
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

αlist = [i for i=10:40]
αlist = αlist./100


thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]

general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",#["Uniform","Optimized"],
    "α" => αlist, # largeur des catégories gaussiennes
    "box_conv" => 0.1#[1.0,0.1]
)

dicts = dict_list(general_args)
@time save_meanlearning_Gauss(dicts[1])
#@time pmap(save_meanlearning_Gauss,dicts)

# way to run back the simulation when failed
# for simulationstemp in dicts
#    simulations = create_simu_from_dict(simulationstemp)
# sname = savename(simulations,"bson")
#  fnametemp = foldername_create(simulations)
#  fname = "sims\\01-script\\results\\$fnametemp"
#
# file = produce_or_load(
#     datadir("$fname"), # path
#     simulationstemp, # container
#     mean_learning_Gaussian_Cat, # function
#     prefix = "" # prefix for savename
#  )
#  fnameparam = "sims\\01-script\\param\\"
#
#  mkpath(datadir(fnameparam))
#
#
#  save(datadir(fnameparam,sname),simulationstemp)
#
# end



thlist_bis=[5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]

general_args_bis = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist_bis, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 15.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => αlist, # largeur des catégories gaussiennes,
    "box_conv" => [1.0,0.1]
)

dicts_bis = dict_list(general_args_bis)
@time pmap(save_meanlearning_Gauss,dicts_bis)

########
# Again simulation crashed, following to run it back
#for simulationstemp in dicts
#    simulations = create_simu_from_dict(simulationstemp)
 #fnametemp = foldername_create(simulations)
 #fname = "sims\\01-script\\results\\$fnametemp"
#
# file = produce_or_load(
#     datadir("$fname"), # path
#     simulationstemp, # container
#     mean_learning_Gaussian_Cat, # function
#     prefix = "" # prefix for savename
 #)
 #fnameparam = "sims\\01-script\\param\\"

 #mkpath(datadir(fnameparam))


# save(datadir(fnameparam,sname),simulationstemp)

#end

#################################################################


# How to check all the results TODO
# res = collect_results!(datadir("sims\\01-script\\results"),subfolders=true)
#
# res = collect_results!(datadir("sims\\01-script\\param"),subfolders=true)
#
#
# res[:1]
#
# res
# How to read a file from this script (TODO)
# simulations=create_simu_from_dict(general_args)
#
# sname = savename(simulations,"bson")
# fnametemp = foldername_create(simulations)
# fname = "sims\\01-script\\results\\$fnametemp"
#
# file = produce_or_load(
#     datadir("$fname"), # path
#     general_args, # container
#     mean_learning_Gaussian_Cat, # function
#     prefix = "" # prefix for savename
# )
#
#
# file
