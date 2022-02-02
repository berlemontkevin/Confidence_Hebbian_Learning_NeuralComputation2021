## In this file I analyze the performance obtained after learning of the task
## This file will save the results in order to be able to observe the fig. with R
#%% Premabule
using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using Distributed
using BSON
using DataFrames
@everywhere include(srcdir("training_network.jl"))

#<<<<<<<



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


#%% Load the weights after learning and compute accuracy


@everywhere function generate_mean_weights_script03(simulation)
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

	dummy1 = zeros(length(mw1[1,:]))
	dummy2 = zeros(length(mw1[1,:]))
	for i=1:10
		dummy1 = dummy1 .+mw1[end+1-i,:]./10
		dummy2 = dummy2 .+mw2[end+1-i,:]./10

	end


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
@everywhere function save_accuracy_from_script03(dict,
				stimlist=[i*0.005 for i=0:100],number=2000)


		w1,w2,mw1,mw2=	generate_mean_weights_script03(dict)
		accuracy=true_performance(dict,stimlist,mw1,mw2,number)


	    box = dict["box_conv"]
		simulations = create_simu_from_dict(dict)
 	    sname = savename(simulations,"bson")
	    fnametemp = foldername_create_accuracy(simulations)
	    fname = "sims\\05-script\\from_script03\\$fnametemp"
	    #fnameparam = "sims\\05-script\\param\\"

	    mkpath(datadir(fname))
	    #mkpath(datadir(fnameparam))

		data = @dict stimlist accuracy dict number
	    save(datadir(fname,sname),data)
	    #save(datadir(fnameparam,sname),diclist)
end

@everywhere function save_accuracy_from_script04(dict,
		stimlist=[i*0.005 for i=0:100],number=2000)

		w1,w2,mw1,mw2=	generate_mean_weights_script04(dict)
		accuracy=true_performance(dict,stimlist,mw1,mw2,number)


	    box = dict["box_conv"]
		simulations = create_simu_from_dict(dict)
 	    sname = savename(simulations,"bson")
	    fnametemp = foldername_create_accuracy(simulations)
	    fname = "sims\\05-script\\from_script04\\$fnametemp"
	    #fnameparam = "sims\\05-script\\param\\"

	    mkpath(datadir(fname))
	    #mkpath(datadir(fnameparam))

		data = @dict stimlist accuracy dict number
	    save(datadir(fname,sname),data)
	    #save(datadir(fnameparam,sname),diclist)
end
@everywhere function save_accuracy_from_script04_Weibull(dict,
		stimlist=[i*0.005 for i=0:100],number=2000)

		w1,w2,mw1,mw2=	generate_mean_weights_script04_Weibull(dict)
		accuracy=true_performance(dict,stimlist,mw1,mw2,number)


	    box = dict["box_conv"]
		simulations = create_simu_from_dict(dict)
 	    sname = savename(simulations,"bson")
	    fnametemp = foldername_create_accuracy(simulations)
	    fname = "sims\\05-script\\from_script04\\$fnametemp"
	    #fnameparam = "sims\\05-script\\param\\"

	    mkpath(datadir(fname))
	    #mkpath(datadir(fnameparam))

		data = @dict stimlist accuracy dict number
	    save(datadir(fname,sname),data)
	    #save(datadir(fnameparam,sname),diclist)
end
#<<<<<<

#%% Threshold to percentage conversion

function threshold_to_percentage(θlist,seuil)

	network = MyAttractorNetwork(0.02,0.3255,0.2609,0.0497,0.0497,0.2609,100.0,0.000641,2.0,seuil)

	timeparameters = MyEulerParameters(0.5,4.0,4.0)
	percentagelist = Float64[]
	tf = MyTransferFunction(270.0,108.0,0.1540)
	deltaR = Float64[]
	for coherence in 0.0:0.01:0.20
	for i=1:1000
	unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
	unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
	inputToUnits = MyInput(coherence,-coherence,30.0,0.00052)
	r1,r2 = simu(unit1,unit2,inputToUnits,network,timeparameters,tf)
	push!(deltaR,abs.(mean(r1[end-4:end]-r2[end-4:end])))
	end
	end

	for θ in θlist
		delta = copy(deltaR)
		delta[deltaR.<θ]=1.0*ones(length(delta[deltaR.<θ]))
		delta[deltaR.>θ]=0.0*ones(length(delta[deltaR.>θ]))

		push!(percentagelist,sum(delta)./length(delta))
	end


	return percentagelist
end




#<<<<



#%% Functions pour les performances après optimisation non linéaire

## function that load the weights
## function that compute the performance (approximation)
function load_optimality(dict)
	N = Int(dict["N"])
	alpha = dict["α"]
	if dict["type_TC"] == "Optimized"
		sname = "C:\\Documents\\1 - PhD Projects\\3-project_coding-categories\\data\\02_optimization\\CatWidth$alpha\\N$N\\OptimalWeights-OptimizedTC.csv"
	else
		sname = "C:\\Documents\\1 - PhD Projects\\3-project_coding-categories\\data\\02_optimization\\CatWidth$alpha\\N$N\\OptimalWeights-UniformTC.csv"

	end
	dfG = CSV.read(sname,header=:false,transpose=true)


	return dfG
end

function performance_approx(x::Float64,codinglayer::MyCodingLayer)
	# parameters used for the performance function are the one at 20Hz

	dummy = 0
	for j=1:length(codinglayer.centers)
		dummy = dummy .+codinglayer.weightsto1[j].*fi(x,codinglayer.centers[j],codinglayer.widths[j])
	end

	beta = 1.35#1.25
	alpha_rescaled = 0.049#0.072
dummy = dummy[1]
		if sign(dummy)<0
			perf =0.5*exp(-exp(beta.*log.(-dummy/alpha_rescaled)))
		else
			perf =(1.0 - 0.5*exp(-exp(beta.*log.(dummy/alpha_rescaled))))
		end


	return perf
end
#<<<<<<<<<<



#%% Création DataFrame pour poids optimisés
using CSV
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end
αlist = [9+i for i=1:31]
αlist = αlist./100

thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
args_nonlin_optimisation = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => 20.0,#thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Uniform",
    "box_conv" => 1.0,
	"α"=>αlist)

dicts = dict_list(args_nonlin_optimisation)

df = DataFrame(N=Int[],alpha=Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts


	temp = load_optimality(dict)
	codinglayer = init_layer(Int(dict["N"]))
	codinglayer.weightsto1 = temp

	for x in [i*0.005 for i=0:100]
		dummy = performance_approx(x,codinglayer)
		push!(df,[dict["N"],dict["α"],dummy,x])

	end

end

CSV.write(datadir("notebooks","06-Script03-NonLinear-Optimisation-TC-Uniform.csv"),df)



thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
args_nonlin_optimisation = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => 20.0,#thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => 1.0,
	"α"=>αlist)

dicts = dict_list(args_nonlin_optimisation)

df = DataFrame(N=Int[],alpha=Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts


	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	temp = load_optimality(dict)
	codinglayer = init_layer(Int(dict["N"]),cats)
	codinglayer.weightsto1 = temp

	for x in [i*0.005 for i=0:100]
		dummy = performance_approx(x,codinglayer)
		push!(df,[dict["N"],dict["α"],dummy,x])
		push!(df,[dict["N"],dict["α"],1.0-dummy,1.0-x])


	end

end

CSV.write(datadir("notebooks","06-Script03-NonLinear-Optimisation-TC-Optimized.csv"),df)



#<<<<<<<<<<<<<



## TODO: amélioration des données avec du bootstrap pour les courbes


###################################################################
#%% Analyse des données du script 06
# Voici les différentes figures voulues:
# Une base de données avec toutes les données
# La différence Opti - Uni
# La différence Best - Opti
# La différence Best - Uni
# Juste les perfs
using CSV
# lecture desp erformances après opti non linéaire
opti_w_uniform = CSV.read(datadir("notebooks","06-Script03-NonLinear-Optimisation-TC-Uniform.csv"))
opti_w_optimized = CSV.read(datadir("notebooks","06-Script03-NonLinear-Optimisation-TC-Optimized.csv"))


# Opti - Uniform
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

αlist = [i for i=10:40]
αlist = αlist./100
thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
percentagelist20 = threshold_to_percentage(thlist,20.0)

args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform"],
    "α" => [0.20],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[10,20,30,40]
)

dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]

	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"

	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => "Uniform",
		"box_conv" => dict["box_conv"],
		"α"=>dict["α"])
		simulations = create_simu_from_dict(dict_temp_Uniform )

		fnametemp = foldername_create_accuracy(simulations)
		sname = savename(simulations,"bson")
		fname = "sims\\05-script\\from_script03\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy] .- dataUniform[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],indice])
push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i],indice])
	end
end
CSV.write(datadir("notebooks","11-Script06-Delta-Optimized-Uniform-20Hz.csv"),df)




# Best - Opti
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Optimized"],
    "α" => [0.20],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[10,20,30,40]
)

df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	temp = load_optimality(dict)
#	dummy = temp[:Column1]./(sqrt(sum(temp[:Column1].*temp[:Column1]))).*0.25
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)
	codinglayer.weightsto1 = temp[:Column1]
dummy_opti = zeros(101)
j=0
	for x in [i*0.005 for i=0:100]
		j=j+1
		dummy_opti[j] = performance_approx(x,codinglayer)
	end
		dummy = dummy_opti .- data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

		for i=1:length(data[:stimlist])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
							dataDict["box_conv"],dummy[i],data[:stimlist][i],indice])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
												dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i],indice])
		end
end
CSV.write(datadir("notebooks","11-Script06-Delta-Best-Optimized-20Hz.csv"),df)

# best - Uni
using Bootstrap
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Uniform",
    "box_conv" => 0.1,
	"α"=>αlist)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],bacth=Float64[])
for dict in dicts

	box = dict["box_conv"]


	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\05-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	temp = load_optimality(dict)
	codinglayer = init_layer(Int(dict["N"]))
	codinglayer.weightsto1 = temp
	dummy_opti = zeros(101)
	j=0
		for x in [i*0.005 for i=0:100]
			j=j+1
			dummy_opti[j] = performance_approx(x,codinglayer)
		end


		dummy = dummy_opti .- data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]
		for i=1:length(data[:stimlist])
		#	templist = Float64[]
		#	for k=1:100
		#		d1 = rand(Binomial(1,dummy_opti[i]))
		#		d2 = rand(Binomial(1,temp[i]))
		#		push!(templist,d1-d2)
		#	end
		#	bs1 = bootstrap(mean,templist,BasicSampling(100))
		#	bci1 = confint(bs1, BasicConfInt(0.95))
		#	if sign(bci1[1][2]) == sign(bci1[1][3])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
							dataDict["box_conv"],dummy[i],data[:stimlist][i],indice])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
												dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i],indice])
		#									else
		#										push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		#															dataDict["box_conv"],0.0,data[:stimlist][i]])
		#										push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		#																				dataDict["box_conv"],0.0,1.0-data[:stimlist][i]])
		#									end
		end
end
CSV.write(datadir("notebooks","11-Script06-Delta-Best-Uniform-20Hz.csv"),df)



## Juste les perfs
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform"],
    "α" => [0.20,0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[10,20,30,40]
)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => "Uniform",
		"box_conv" => dict["box_conv"],
		"α"=>dict["α"])
		indice = dict["batch"]

		simulations = create_simu_from_dict(dict_temp_Uniform)
		sname = savename(simulations,"bson")

		fnametemp = foldername_create_accuracy(simulations)

		fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"],indice])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"],indice])

	end
end
CSV.write(datadir("notebooks","11-Script06-Perf-20Hz.csv"),df)



# Simple perf with deduction from no confidence
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform"],
    "α" => [0.20,0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[10,20,30,40]
)
dicts = dict_list(args_temp)

df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => 20.0, # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"box_conv" => dict["box_conv"],
		"α"=>dict["α"])
		simulations = create_simu_from_dict(dict_temp_Uniform )

		fnametemp = foldername_create_accuracy(simulations)
		sname = savename(simulations,"bson")
		fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy] - dataUniform[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"],indice])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"],indice])

	end
end
CSV.write(datadir("notebooks","11-Script06-Delta-Perf-20Hz.csv"),df)


#<<<<<<<<<<<<<
#######################################################################


#%% Data for figure influence of batch on learning
αlist = [i for i=10:40]
αlist = αlist./100
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
general_args = Dict(
    "N" => [25],#]Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[5,10,20,30,40,50,60,70,80,90,100]
)
dicts = dict_list(general_args)


df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]


		dummy = data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"],indice])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"],indice])

	end
end
CSV.write(datadir("notebooks","11-Script06-Perf-Full-Batch-20Hz.csv"),df)

general_args = Dict(
    "N" => [10,15,20],#]Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform"],
    "α" => [0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[5,10,20,30,40,50,60,70,80,90,100]
)
dicts = dict_list(general_args)


df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]


		dummy = data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"],indice])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"],indice])

	end
end
CSV.write(datadir("notebooks","11-Script06-Perf-Full-Batchv2-20Hz.csv"),df)

# Simple perf with deduction from no confidence
args_temp = Dict(
    "N" => [10,15,20,25,30,35], # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.25],#0.10 lancé plus tot#αlist, # largeur des catégories gaussiennes
    "box_conv" => [0.1],
	"batch" =>[5,10,20,30,40,50,60,70,80,90,100]
)
dicts = dict_list(args_temp)

df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[],batch=Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")

	fnametemp = foldername_create_accuracy(simulations)

	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => 20.0, # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"box_conv" => dict["box_conv"],
		"α"=>dict["α"])
		simulations = create_simu_from_dict(dict_temp_Uniform )

		fnametemp = foldername_create_accuracy(simulations)
		sname = savename(simulations,"bson")
		fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy] - dataUniform[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"],indice])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"],indice])

	end
end
CSV.write(datadir("notebooks","11-Script06-Delta-Perfv2-20Hz.csv"),df)

#<<<<<<<<<<<<<<<<<<<


#%% Analyze confidence depending on batch during learning
# batch 5 pour commencer
using CSV

dict = Dict(
    "N" => 20, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => 14.0, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => 0.1,
	"α"=>0.25)

mw1,mw2,dummy1,dummy2 = generate_mean_weights_script03(dict)
df = DataFrame(Indice = Float64[], batch = Int[], weight = Float64[])

fname = "notebooks"



function generate_confidence(c::MyCodingLayer,indice,mw1,mw2)
templist = Int[]
templist2 = []
c = init_layer(Int(dict["N"]),MyCatGauss(width=[0.25,0.25]))
c.weightsto1 = mw1[indice,:]
c.weightsto2 = mw2[indice,:]
for stim in 0.2:0.01:0.5
	temp = 0
	temp2 = []
 for k=1:1000
	 network = MyAttractorNetwork(0.02,0.3255,0.2609,0.0497,0.0497,0.2609,100.0,0.000641,2.0,20.0)
	 tf = MyTransferFunction(270.0,108.0,0.1540)

	 timeparameters = MyEulerParameters(0.5,4.0,4.0)
	response = responsetostim(stim,c)
	i1,i2=inputCodingtoNetwork(c,response)
	unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
	unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
	inputToUnits = MyInput(i1,i2,30.0,0.00052)
	r1,r2 = simu(unit1,unit2,inputToUnits,network,timeparameters,tf)
		if abs(mean(r1[end-3:end])-mean(r2[end-3:end])) >14
			temp = temp+1
		end
		push!(temp2,abs(mean(r1[end-3:end])-mean(r2[end-3:end])))
 end
	push!(templist,temp)
	push!(templist2,mean(temp2))

	end
	return templist./1000,templist2
end


a,b = generate_confidence(c,1,mw1,mw2)
plot(0.2:0.01:0.5,b)

a,b = generate_confidence(c,2,mw1,mw2)
plot!(0.2:0.01:0.5,b)

a,b = generate_confidence(c,3,mw1,mw2)
plot!(0.2:0.01:0.5,b)

a,b = generate_confidence(c,4,mw1,mw2)
plot!(0.2:0.01:0.5,b)

a,b = generate_confidence(c,5,mw1,mw2)
plot!(0.2:0.01:0.5,b)


a,b = generate_confidence(c,20,mw1,mw2)
plot!(0.2:0.01:0.5,b)

a,b = generate_confidence(c,10,mw1,mw2)
plot!(0.2:0.01:0.5,b)

a,b = generate_confidence(c,100,mw1,mw2)
plot!(0.2:0.01:0.5,b)

a,b = generate_confidence(c,80,mw1,mw2)
plot!(0.2:0.01:0.5,b)


0.2:0.01:0.5

plot(0.2:0.01:0.5,b)

df = DataFrame(Stim = Float64[], batch = Int[], mean_confidence = Float64[], percentage_rejet=Float64[])


for i=1:100
	a,b = generate_confidence(c,i,mw1,mw2)
for j=1:length(a)
push!(df,[0.2+0.01*(j-1),i,b[j],a[j]])
end

end

CSV.write(datadir(fname,"11-optimized-confidence-batch.csv"),df)


using Plots
pyplot()






#<<<<<<<<<<<<<<



#%% Save scripts as with end of learning at 10k, but at 2k for specifc values of parameters

#%% Juste les perfs
args_temp = Dict(
    "N" => [10,20,30], # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "box_conv" => [0.1,1.0],
	"α"=>[0.25],
	"batch"=>[20])
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt =Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => dict["type_TC"],
		"α" => dict["α"],#αlist, # largeur des catégories gaussiennes
		"box_conv" => 0.1
	)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
	dummy = data[:accuracy]
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"]])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"]])

	end
end

CSV.write(datadir("notebooks","11-Batch20-Perf-20Hz.csv"),df)
#<<<<<<<

#%% COmparaison between coding

# Opti - Uniform
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

αlist = [i for i=10:40]
αlist = αlist./100
thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
percentagelist20 = threshold_to_percentage(thlist,20.0)

args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => [0.1],
	"α"=>[0.2,0.3,0.4],
	"batch"=>20)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = dict["batch"]
	dictt = Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => "Optimized",
		"box_conv" => dict["box_conv"],
		"α"=>dict["α"])
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
		"N" => dict["N"], # nombre de neurones dans le codage
		"T_simu" => 10000,
		"ΔTsave" => 100, # pas de temps de sauvegarde
		"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
		"seuil_decision" => 20.0, # seuil de la décision
		"type_TC" => "Uniform",
		"box_conv" => dict["box_conv"],
		"α"=>dict["α"])
		simulations = create_simu_from_dict(dict_temp_Uniform )

		fnametemp = foldername_create_accuracy(simulations)
		sname = savename(simulations,"bson")
		fname = "sims\\06-script\\from_script03\\Batch=$indice\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy] .- dataUniform[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i]])
push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i]])
	end
end
CSV.write(datadir("notebooks","11-Batch20-Delta-Optimized-Uniform-20Hz.csv"),df)

#<<<<<<<<

#<<<<<<<<<<<
