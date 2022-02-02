### Currently Windows version

#%% Premabule
using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using Distributed
using BSON
using DataFrames
#<<<<<<<

#%% Multi processeurs
@everywhere using DrWatson
@everywhere quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
@everywhere using Distributed,BSON,DataFrames
@everywhere include(srcdir("training_network.jl"))
#<<<<<<<<

# In this file:
#TODO

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
    alpha = dummy.α#tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    fname = "$Modele\\$tcatname\\"
    return fname
end

@everywhere function foldername_create_accuracy(dummy::ListParameters_θconfidenceSimulation)
    seuil = dummy.seuil_decision
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = dummy.α#tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    fname = "$tcatname\\"
    return fname
end

#<<<<<<<<<<<<<<<<<<<<<<<<<<<


#%% Functions for Gaussian Categories
@everywhere function create_simu_from_dict(diclist)
	# create a structure simulaiton from a dictionnary of parameters
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
	# currently the box_conv timescale needs to be given outside the function
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


@everywhere function save_meanlearning_Gauss_script03(diclist)
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


@everywhere function save_meanlearning_Gauss_script04(diclist)
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

#<<<<<<<<<<


#%% Functions for Weibull Categories

@everywhere function create_simu_Weibull_from_dict(diclist)
	@unpack N,T_simu,ΔTsave,Seuil_confidence,seuil_decision,type_TC,α,box_conv = diclist


    α_weibull = 0.1
    β_weibull = 1.61
    gammag = 0.82
    alphag = sqrt((1-2.0.*(α_weibull+0.5))/(2*log((1-gammag)/gammag)))

    mycatWeib = MyCatWeibull(alpha=α_weibull,beta=β_weibull, width = alphag)
    simulations = ListParameters_θconfidenceSimulation(N=N,
    T_simu=T_simu,
    ΔTsave=ΔTsave,
    Seuil_confidence=Seuil_confidence,
    seuil_decision=seuil_decision,
    type_TC=type_TC,
    α = α,#alphag, Approximation
    cats=mycatWeib,
	box_conv=box_conv)

    return simulations
end


@everywhere function mean_learning_Weibull_Cat(diclist,box=1.0,rep=10)
    list_weights_to_unit1 = [] # chaque elt correspond à une simu
    list_weights_to_unit2 = [] # chaque elt correspond à une simu
	box = diclist["box_conv"]

    simulations = create_simu_Weibull_from_dict(diclist)

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


@everywhere function save_meanlearning_Weib(diclist)

    box = diclist["box_conv"]
    dummy = mean_learning_Weibull_Cat(diclist,box)
    @unpack list_weights_to_unit1,list_weights_to_unit2,simulations = dummy
    sname = savename(simulations,"bson")
    fnametemp = foldername_create(simulations,box)
    fname = "sims\\04-script\\results\\$fnametemp"
    fnameparam = "sims\\04-script\\param\\"

    mkpath(datadir(fname))
    mkpath(datadir(fnameparam))

    save(datadir(fname,sname),dummy)
    save(datadir(fnameparam,sname),diclist)
    return true
end


#<<<<<<<<<<<<


#%% Load the weights after learning and compute accuracy

@everywhere function generate_mean_weights_script04(simulation)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL,box)
	fname = "sims\\04-script\\results\\$fnametemp"
	file = produce_or_load(
	datadir("$fname"), # path
	simulation, # container
	mean_learning_Gaussian_Cat, # function
	prefix = "" # prefix for savename
	)


	weights1_temp = file[1][:list_weights_to_unit1]
	weights2_temp = file[1][:list_weights_to_unit2]

	mw1 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)
	mw2 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)

	for i=1:10

		for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)
			for k=1:simulationL.N
				mw1[j,k] += weights1_temp[i][j,k] .- weights1_temp[i][j,end+1-k]
				mw2[j,k] += weights2_temp[i][j,k] .- weights2_temp[i][j,end+1-k]
			end
		end

	end
	mw1 = mw1 .- mw2
	mw2 = -mw1
	mw1 = mw1./40
	mw2 = mw2./40

	for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)
		temp_w1 = sum(mw1[j,:].*mw1[j,:])
		temp_w2 = sum(mw2[j,:].*mw2[j,:])
		for k=1:simulationL.N
			mw1[j,k] = norm_w*mw1[j,k]/sqrt(temp_w1)
			mw2[j,k] = norm_w*mw2[j,k]/sqrt(temp_w2)

		end
	end

	dummy1 = zeros(length(mw1[1,:]))
	dummy2 = zeros(length(mw1[1,:]))
	for i=1:10
		dummy1 = dummy1 .+mw1[end+1-i,:]./10
		dummy2 = dummy2 .+mw2[end+1-i,:]./10

	end


	return mw1,mw2,dummy1,dummy2
end

@everywhere function generate_mean_weights_script04_Weibull(simulation)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL,box)
	fname = "sims\\04-script\\results\\$fnametemp"
	file = produce_or_load(
	datadir("$fname"), # path
	simulation, # container
	mean_learning_Weibull_Cat, # function
	prefix = "" # prefix for savename
	)


	weights1_temp = file[1][:list_weights_to_unit1]
	weights2_temp = file[1][:list_weights_to_unit2]

	mw1 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)
	mw2 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)

	for i=1:10

		for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)
			for k=1:simulationL.N
				mw1[j,k] += weights1_temp[i][j,k] .- weights1_temp[i][j,end+1-k]
				mw2[j,k] += weights2_temp[i][j,k] .- weights2_temp[i][j,end+1-k]
			end
		end

	end
	mw1 = mw1 .- mw2
	mw2 = -mw1
	mw1 = mw1./40
	mw2 = mw2./40

	for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)
		temp_w1 = sum(mw1[j,:].*mw1[j,:])
		temp_w2 = sum(mw2[j,:].*mw2[j,:])
		for k=1:simulationL.N
			mw1[j,k] = norm_w*mw1[j,k]/sqrt(temp_w1)
			mw2[j,k] = norm_w*mw2[j,k]/sqrt(temp_w2)

		end
	end

	dummy1 = zeros(length(mw1[1,:]))
	dummy2 = zeros(length(mw1[1,:]))
	for i=1:10
		dummy1 = dummy1 .+mw1[end+1-i,:]./10
		dummy2 = dummy2 .+mw2[end+1-i,:]./10

	end


	return mw1,mw2,dummy1,dummy2
end

@everywhere function generate_mean_weights_script03(simulation,indice)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL,box)
	fname = "sims\\03-script\\results\\$fnametemp"
	file = produce_or_load(
	datadir("$fname"), # path
	simulation, # container
	mean_learning_Gaussian_Cat, # function
	prefix = "" # prefix for savename
	)


	weights1_temp = file[1][:list_weights_to_unit1]
	weights2_temp = file[1][:list_weights_to_unit2]

	mw1 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)
	mw2 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)

	for i=1:10

		for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)
			for k=1:simulationL.N
				mw1[j,k] += weights1_temp[i][j,k] .- weights1_temp[i][j,end+1-k]
				mw2[j,k] += weights2_temp[i][j,k] .- weights2_temp[i][j,end+1-k]
			end
		end

	end
	mw1 = mw1 .- mw2
	mw2 = -mw1
	mw1 = mw1./40
	mw2 = mw2./40

	for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)
		temp_w1 = sum(mw1[j,:].*mw1[j,:])
		temp_w2 = sum(mw2[j,:].*mw2[j,:])
		for k=1:simulationL.N
			mw1[j,k] = norm_w*mw1[j,k]/sqrt(temp_w1)
			mw2[j,k] = norm_w*mw2[j,k]/sqrt(temp_w2)

		end
	end

	dummy1 = mw1[indice,:]
	dummy2 = mw2[indice,:]

	return mw1,mw2,dummy1,dummy2
end


@everywhere function true_performance(simulationstemp,
				stimlist::Array{Float64,1},wto1::Array{Float64,1},
					wto2::Array{Float64,1},tempnumber = 1000)

		simulations = create_simu_from_dict(simulationstemp)
		box = simulationstemp["box_conv"]

		network = MyAttractorNetwork(threshold=simulations.seuil_decision)
	    timeparameters = MyEulerParameters(0.5,4.0,4.0)
	    tf = MyTransferFunction(270.0,108.0,0.1540)

	    if simulations.type_TC == "Uniform"
	        codinglayer =init_layer(simulations.N)
	    elseif simulations.type_TC == "Hebbian"
	        codinglayer =init_layer_Hebbian(simulations.N,simulations.cats)
	    else
	        codinglayer =init_layer(simulations.N,simulations.cats)
	    end
		codinglayer.weightsto1 = wto1
		codinglayer.weightsto2 = wto2

	    inputToUnits = MyInput(0.0,0.0,30.0,0.00052)

		perf_list = zeros(length(stimlist))
	for k=1:length(stimlist)
			stim=stimlist[k]
			dummy=0.0
		@simd	for i=1:tempnumber
				unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
				unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
	        r1,r2 = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
				if mean(r1.-r2) >0
					dummy+=1
				end
			end
			perf_list[k] = dummy/tempnumber
		end
			return perf_list
end

#<<<<<<<<<<

#%% Save accuracy from a dictionnary
# explain the structure of data
# Pas besoin de sauvegarder tous les choix car au final la décision correspond à une
# loi binomiale. On peut donc retirer les choix et bootstrapper.
#
# Solution de facilité : Un fichier par paramètres
# enregistre sous .bson un dict avec les stim et les perfs
@everywhere function save_accuracy_from_script03_batch(dicttemp,
				stimlist=[i*0.005 for i=0:100],number=2000)

				indice = dicttemp["batch"]
				dict =Dict(
				    "N" => dicttemp["N"], # nombre de neurones dans le codage
				    "T_simu" => 10000,
				    "ΔTsave" => 100, # pas de temps de sauvegarde
				    "Seuil_confidence" => dicttemp["Seuil_confidence"], # seuil sur la confiancep our l'appr.
				    "seuil_decision" => 20.0, # seuil de la décision
				    "type_TC" => dicttemp["type_TC"],
				    "α" => dicttemp["α"],#αlist, # largeur des catégories gaussiennes
				    "box_conv" => 0.1
				)
		w1,w2,mw1,mw2=	generate_mean_weights_script03(dict,indice)
		accuracy=true_performance(dict,stimlist,mw1,mw2,number)


	    box = dict["box_conv"]
		simulations = create_simu_from_dict(dict)
 	    sname = savename(simulations,"bson")
	    fnametemp = foldername_create_accuracy(simulations)
	    fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	    #fnameparam = "sims\\05-script\\param\\"

	    mkpath(datadir(fname))
	    #mkpath(datadir(fnameparam))

		data = @dict stimlist accuracy dict number
	    save(datadir(fname,sname),data)
	    #save(datadir(fnameparam,sname),diclist)
end
#<<<<<<



# Currently the performances are computed at a normalisation of 0.25
#%% Accuracy for script 03 ---- 20Hz
#αlist = [i for i=10:40]
#αlist = αlist./100
#Nlist = zeros(31)
#for i=1:31#
#    Nlist[i] = 9+i
#end

#thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
#general_args = Dict(
#    "N" => Nlist, # nombre de neurones dans le codage
#    "T_simu" => 10000,
#    "ΔTsave" => 100, # pas de temps de sauvegarde
#    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
#    "seuil_decision" => 20.0, # seuil de la décision
#    "type_TC" => ["Uniform","Optimized"],
#    "α" => [0.20,0.25,0.30,0.40],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
#    "box_conv" => [0.1],
#	"batch" =>[10,20,30,40]
#)
#dicts = dict_list(general_args)

#@time save_accuracy_from_script03_batch(dicts[1])

#@time pmap(save_accuracy_from_script03_batch,dicts)

#<<<<<<<<<<
#
# #%% Full data accuracy at N=20 and the other (for supp. fig), alpha = 0.25
# αlist = [i for i=10:40]
# αlist = αlist./100
# Nlist = zeros(31)
# for i=1:31
#     Nlist[i] = 9+i
# end
#
# thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
# general_args = Dict(
#     "N" => [25],#]Nlist, # nombre de neurones dans le codage
#     "T_simu" => 10000,
#     "ΔTsave" => 100, # pas de temps de sauvegarde
#     "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
#     "seuil_decision" => 20.0, # seuil de la décision
#     "type_TC" => ["Uniform","Optimized"],
#     "α" => [0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
#     "box_conv" => [0.1],
# 	"batch" =>[5,10,20,30,40,50,60,70,80,90,100]
# )
# dicts = dict_list(general_args)
#
# @time save_accuracy_from_script03_batch(dicts[1])
#
# @time pmap(save_accuracy_from_script03_batch,dicts)
#
# general_args = Dict(
#     "N" => [10,15,20,30,35],#]Nlist, # nombre de neurones dans le codage
#     "T_simu" => 10000,
#     "ΔTsave" => 100, # pas de temps de sauvegarde
#     "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
#     "seuil_decision" => 20.0, # seuil de la décision
#     "type_TC" => ["Uniform","Optimized"],
#     "α" => [0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
#     "box_conv" => [0.1],
# 	"batch" =>[5,10,20,30,40,50,60,70,80,90,100]
# )
# dicts = dict_list(general_args)
#
# @time save_accuracy_from_script03_batch(dicts[1])
#
# @time pmap(save_accuracy_from_script03_batch,dicts)
# #<<<<<<<<<<<
#
# #%% New bacth save
# αlist = [i for i=10:40]
# αlist = αlist./100
# Nlist = zeros(31)
# for i=1:31
#     Nlist[i] = 9+i
# end
#
# thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
# general_args = Dict(
#     "N" => [10,15,20,30,35],#]Nlist, # nombre de neurones dans le codage
#     "T_simu" => 10000,
#     "ΔTsave" => 100, # pas de temps de sauvegarde
#     "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
#     "seuil_decision" => 20.0, # seuil de la décision
#     "type_TC" => ["Uniform","Optimized"],
#     "α" => [0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
#     "box_conv" => [0.1],
# 	"batch" =>[1,2,3,4,5]
# )
# dicts = dict_list(general_args)
#
# @time save_accuracy_from_script03_batch(dicts[1])
#
# @time pmap(save_accuracy_from_script03_batch,dicts)
# #<<<<<<<<<

#%% New bacth save
αlist = [i for i=10:40]
αlist = αlist./100
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
general_args = Dict(
    "N" => Nlist,#]Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.3,0.4],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[20]
)
dicts = dict_list(general_args)

@time save_accuracy_from_script03_batch(dicts[1])

@time pmap(save_accuracy_from_script03_batch,dicts)
#<<<<<<<<<
