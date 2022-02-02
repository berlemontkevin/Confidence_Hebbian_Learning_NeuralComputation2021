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
############
# Default values of savenames
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
:α)

@everywhere function foldername_create(dummy::ListParameters_θconfidenceSimulation)
    seuil = dummy.seuil_decision
    Modele = "Decision_Seuil=$seuil"
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    f = "$Modele\\$tcatname\\"
    return f
end


######
# Function to run the simulations
############
@everywhere function create_simu_from_dict(diclist)
    @unpack N,T_simu,ΔTsave,Seuil_confidence,seuil_decision,type_TC,α = diclist

    mycatGauss = MyCatGauss(width=[α,α])
    simulations = ListParameters_θconfidenceSimulation(N=N,
    T_simu=T_simu,
    ΔTsave=ΔTsave,
    Seuil_confidence=Seuil_confidence,
    seuil_decision=seuil_decision,
    type_TC=type_TC,
    α = α,
    cats=mycatGauss)

    return simulations
end

@everywhere function mean_learning_Gaussian_Cat(diclist,rep=10)
    list_weights_to_unit1 = [] # chaque elt correspond à une simu
    list_weights_to_unit2 = [] # chaque elt correspond à une simu

    simulations = create_simu_from_dict(diclist)

    for i=1:rep
        Weightsto1,Weightsto2 = learning(simulations)
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

    dummy = mean_learning_Gaussian_Cat(diclist)
    @unpack list_weights_to_unit1,list_weights_to_unit2,simulations = dummy
    sname = savename(simulations,"bson")
    fnametemp = foldername_create(simulations)
    fname = "sims\\01-script\\results\\$fnametemp"
    fnameparam = "sims\\01-script\\param\\"

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

αlist = zeros(31)
for i=1:31
    αlist[i] = (9+i)*0.01
end

thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]

general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => αlist # largeur des catégories gaussiennes
)

dicts = dict_list(general_args)

@time pmap(save_meanlearning_Gauss,dicts)

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
    "α" => αlist # largeur des catégories gaussiennes
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
