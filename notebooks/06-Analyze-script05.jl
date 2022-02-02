## In this file I analyze the performance obtained after learning of the task
## This file will save the results in order to be able to observe the fig. with R
#%% Premabule
using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using Distributed
using BSON
using DataFrames
  include(srcdir("training_network.jl"))

#<<<<<<<


#%% Default values of savenames
#########

DrWatson.default_prefix(e::ListParameters_θconfidenceSimulation) = ""
# dossier de type
# sims/scriptname/results/Modele/TypeCat/Seuil
DrWatson.allaccess(::ListParameters_θconfidenceSimulation) = (:N,
:T_simu,
:ΔTsave,
:Seuil_confidence,
:seuil_decision,
:type_TC,
:α,
:box_conv)

  function foldername_create(dummy::ListParameters_θconfidenceSimulation,box::Float64)
    seuil = dummy.seuil_decision
    Modele = "Decision_Seuil=$seuil"
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = dummy.α#tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    fname = "$Modele\\$tcatname\\"
    return fname
end

  function foldername_create_accuracy(dummy::ListParameters_θconfidenceSimulation)
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
  function create_simu_from_dict(diclist)
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

  function mean_learning_Gaussian_Cat(diclist,box=1.0,rep=10)
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


  function save_meanlearning_Gauss_script03(diclist)
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


  function save_meanlearning_Gauss_script04(diclist)
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

  function create_simu_Weibull_from_dict(diclist)
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


  function mean_learning_Weibull_Cat(diclist,box=1.0,rep=10)
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


  function save_meanlearning_Weib(diclist)

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

  function generate_mean_weights_script04(simulation)
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
			mw1[j,k] = norm_w*mw1[j,k]/temp_w1
			mw2[j,k] = norm_w*mw2[j,k]/temp_w2

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

  function generate_mean_weights_script04_Weibull(simulation)
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
			mw1[j,k] = norm_w*mw1[j,k]/temp_w1
			mw2[j,k] = norm_w*mw2[j,k]/temp_w2

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

  function generate_mean_weights_script03(simulation)
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


  function true_performance(simulationstemp,
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
  function save_accuracy_from_script03(dict,
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

  function save_accuracy_from_script04(dict,
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
  function save_accuracy_from_script04_Weibull(dict,
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


#%% Analyze scirpt 04 performances
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
for i=1:11
print(percentagelist20[i])
print("\n")
end
percentagelist20 = threshold_to_percentage(thlist,20.0)

general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Hebbian","Uniform","Optimized"],
    "box_conv" => [1.0,0.1],
	"α"=>0.257)

dicts = dict_list(general_args)

### On veut
# performances pour chaque type de TC
# Différence de performances également

df_Uniform = DataFrame(N=Int[],Seuil_confidence=Float64[],percentageConfidence=Float64[],
box_conv=Float64[],Stim=Float64[],Perf=Float64[])


dict_Uniform = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Uniform",
    "box_conv" => [1.0,0.1],
	"α"=>0.257)

	list_Uniform = dict_list(dict_Uniform)

for dict in list_Uniform
	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script04\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df_Uniform,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],data[:stimlist][i],
						data[:accuracy][i]])
    push!(df_Uniform,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
                    						dataDict["box_conv"],1.0-data[:stimlist][i],
                    						data[:accuracy][i]])
	end
end

df_Uniform

using CSV
CSV.write(datadir("notebooks","06-Script04-TC-Uniform.csv"),df_Uniform)

## Optimized TC

df_Opti = DataFrame(N=Int[],Seuil_confidence=Float64[],percentageConfidence=Float64[],
box_conv=Float64[],Stim=Float64[],Perf=Float64[])
dict_Opti = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => [1.0,0.1],
	"α"=>0.257)

	list_Opti = dict_list(dict_Opti)
for dict in list_Opti
	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script04\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df_Opti,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],data[:stimlist][i],
						data[:accuracy][i]])
	push!(df_Opti,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-data[:stimlist][i],
											data[:accuracy][i]])
	end
end
CSV.write(datadir("notebooks","06-Script04-TC-Optimized.csv"),df_Opti)

## Hebbian TC

df_Opti = DataFrame(N=Int[],Seuil_confidence=Float64[],percentageConfidence=Float64[],
box_conv=Float64[],Stim=Float64[],Perf=Float64[])
dict_Opti = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Hebbian",
    "box_conv" => [1.0,0.1],
	"α"=>0.257)

	list_Opti = dict_list(dict_Opti)
for dict in list_Opti
	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script04\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df_Opti,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],data[:stimlist][i],
						data[:accuracy][i]])
		push!(df_Opti,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-data[:stimlist][i],
											data[:accuracy][i]])
	end
end
CSV.write(datadir("notebooks","06-Script04-TC-Hebbian.csv"),df_Opti)


df_Opti .== df_Uniform

## Delta between uniform and Optimized
df_Delta_OU = DataFrame(N=Int[],Seuil_confidence=Float64[],percentageConfidence=Float64[],
box_conv=Float64[],Stim=Float64[],Perf=Float64[])
dict_Delta_OU = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => [1.0,0.1],
	"α"=>0.257)

	list_Delta_OU = dict_list(dict_Delta_OU)
for dict in list_Delta_OU
	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script04\\$fnametemp"
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
		"α"=>0.257)
		simulations = create_simu_from_dict(dict_temp_Uniform )

		fnametemp = foldername_create_accuracy(simulations)
		sname = savename(simulations,"bson")
		fname = "sims\\05-script\\from_script04\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy] .- dataUniform[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df_Delta_OU,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],data[:stimlist][i],
						dummy[i]])
	push!(df_Delta_OU,[dataDict["N"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-data[:stimlist][i],
											dummy[i]])
	end
end
CSV.write(datadir("notebooks","06-Script04-TC-Delta-OU.csv"),df_Delta_OU)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<

#%% Functions pour les performances après optimisation non linéaire

## function that load the weights
## function that compute the performance (approximation)
function load_optimality(dict,a=true)
	# right now we only ocmpare to the best possible network
	N = Int(dict["N"])
	alpha = dict["α"]
	if a
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
#%% Analyse des données du script 03
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
    "type_TC" => "Optimized",
    "box_conv" => [0.1,1.0],
	"α"=>αlist)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts

	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script03\\$fnametemp"
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
						dataDict["box_conv"],dummy[i],data[:stimlist][i]])
push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i]])
	end
end
CSV.write(datadir("notebooks","06-Script03-Delta-Optimized-Uniform-20Hz.csv"),df)




# Best - Opti
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => 0.1,
	"α"=>αlist)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts

	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script03\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	temp = load_optimality(dict)
	dummy = temp[!,:Column1]./(sqrt(sum(temp[!,:Column1].*temp[!,:Column1]))).*0.25
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)
	codinglayer.weightsto1 = dummy
dummy_opti = zeros(101)
j=0
	for x in [i*0.005 for i=0:100]
		j=j+1
		dummy_opti[j] = performance_approx(x,codinglayer)
	end
		dummy =  dummy_opti .- data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

		for i=1:length(data[:stimlist])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
							dataDict["box_conv"],dummy[i],data[:stimlist][i]])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
												dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i]])
		end
end
CSV.write(datadir("notebooks","06-Script03-Delta-Best-Optimized-20Hz.csv"),df)

# best - Uni
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
				box_conv = Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts

	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script03\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	temp = load_optimality(dict)
	dummy = temp[!,:Column1]./(sqrt(sum(temp[!,:Column1].*temp[!,:Column1]))).*0.25

	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)
	codinglayer.weightsto1 = dummy
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
							dataDict["box_conv"],dummy[i],data[:stimlist][i]])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
												dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i]])

		end
end
CSV.write(datadir("notebooks","06-Script03-Delta-Best-Uniform-20Hz.csv"),df)



## Juste les perfs
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "box_conv" => [0.1,1.0],
	"α"=>αlist)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[])
for dict in dicts

	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script03\\$fnametemp"
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
		dummy = data[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"]])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"]])

	end
end
CSV.write(datadir("notebooks","06-Script03-Perf-20Hz.csv"),df)

# Simple perf with deduction from no confidence
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[],Stim = Float64[],TC=String[])
for dict in dicts

	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script03\\$fnametemp"
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
		fname = "sims\\05-script\\from_script03\\$fnametemp"
		dataUniform = load(datadir(fname,sname))
		dataDictUniform = dataUniform[:dict]
		dummy = data[:accuracy] - dataUniform[:accuracy]
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
						dataDict["box_conv"],dummy[i],data[:stimlist][i],dataDict["type_TC"]])
	push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],1.0-dummy[i],1.0-data[:stimlist][i],dataDict["type_TC"]])

	end
end
CSV.write(datadir("notebooks","06-Script03-Delta-Perf-20Hz.csv"),df)


#<<<<<<<<<<<<<
#######################################################################


#%% Analysis of integrale of probability of performances

# Best - Opti
args_temp = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Optimized",
    "box_conv" => 0.1,
	"α"=>αlist)
dicts = dict_list(args_temp)
df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
				box_conv = Float64[],Perf=Float64[])
				for dict in dicts

					box = dict["box_conv"]
					simulations = create_simu_from_dict(dict)
					sname = savename(simulations,"bson")
					fnametemp = foldername_create_accuracy(simulations)
					fname = "sims\\05-script\\from_script03\\$fnametemp"
					data = load(datadir(fname,sname))
					dataDict = data[:dict]

					temp = load_optimality(dict)
					dummy = temp[!,:Column1]./(sqrt(sum(temp[!,:Column1].*temp[!,:Column1]))).*0.25
					cats = MyCatGauss(width=[dict["α"],dict["α"]])
					codinglayer = init_layer(Int(dict["N"]),cats)
					codinglayer.weightsto1 = dummy
				dummy_opti = zeros(101)
				j=0
				dummya = 0.0
					for x in [i*0.005 for i=0:100]
						j=j+1
if j>40
						dummy_opti[j] = performance_approx(x,codinglayer)
						dummya = dummya + 0.005*(dummy_opti[j] - data[:accuracy][j])
end
					end
					dummy2 = 0.0
					temp = load_optimality(dict,false)
					dummy = temp[!,:Column1]./(sqrt(sum(temp[!,:Column1].*temp[!,:Column1]))).*0.25

					cats = MyCatGauss(width=[dict["α"],dict["α"]])
					codinglayer = init_layer(Int(dict["N"]))
					codinglayer.weightsto1 = dummy
					j=0
						for x in [i*0.005 for i=0:100]
							j=j+1
							if j>40
							dummy_opti[j] = performance_approx(x,codinglayer)
							dummy2 = dummy2 + 0.005*(dummy_opti[j] - data[:accuracy][j])
end
						end
						dummyb = dummya
						percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

						push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
											dataDict["box_conv"],dummyb])

				end
CSV.write(datadir("notebooks","06-Script03-Intergale-Delta-Best-Optimized-20Hz.csv"),df)



# Best - Uni
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
				box_conv = Float64[],Perf=Float64[])
for dict in dicts

	box = dict["box_conv"]
	simulations = create_simu_from_dict(dict)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\05-script\\from_script03\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	temp = load_optimality(dict)
	dummy = temp[!,:Column1]./(sqrt(sum(temp[!,:Column1].*temp[!,:Column1]))).*0.25
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)
	codinglayer.weightsto1 = dummy
dummy_opti = zeros(101)
j=0
dummya = 0.0
	for x in [i*0.005 for i=0:100]
		j=j+1
if j>40
		dummy_opti[j] = performance_approx(x,codinglayer)
		dummya = dummya + 0.005*(dummy_opti[j] - data[:accuracy][j])
end
	end
	dummy2 = 0.0
	temp = load_optimality(dict,false)
	dummy = temp[!,:Column1]./(sqrt(sum(temp[!,:Column1].*temp[!,:Column1]))).*0.25

	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]))
	codinglayer.weightsto1 = dummy
	j=0
		for x in [i*0.005 for i=0:100]
			j=j+1
if j>40
			dummy_opti[j] = performance_approx(x,codinglayer)
			dummy2 = dummy2 + 0.005*(dummy_opti[j] - data[:accuracy][j])
end
		end
		dummyb = max(dummya,dummy2)
		percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
							dataDict["box_conv"],dummyb])

end
CSV.write(datadir("notebooks","06-Script03-Intergale-Delta-Best-Uniform-20Hz.csv"),df)


#<<<<<<<<<<


#%%% Test de quelques figures
df = CSV.read(datadir("notebooks","06-Script03-Delta-Optimized-Uniform-20Hz.csv"))


df

using Gadfly
using ColorSchemes

font_panel = Theme(
	major_label_font="Arial",
	minor_label_font="Arial",
	major_label_font_size=16pt,
	minor_label_font_size=14pt
)

#%% alpha 0.2 avec conf et sans conf
t = df[df[:alpha] .== 0.2,:]
dfplot = t[t[:Seuil_confidence] .== 15.0,:]
dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
M = maximum(dfplot[:Perf])
m = minimum(dfplot[:Perf])
mtot = max(M,-m)

p=plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5))
img = SVG(plotsdir("05-Batch100","Delta-Best-Uni-Conf15-alpha02.svg"), 6inch, 4inch)
draw(img, p)

t = df[df[:alpha] .== 0.2,:]
dfplot = t[t[:Seuil_confidence] .== 20.0,:]
dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
M = maximum(dfplot[:Perf])
m = minimum(dfplot[:Perf])
mtot = max(M,-m)

p=plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5))
img = SVG(plotsdir("05-Batch100","Delta-Best-Uni-Conf20-alpha02.svg"), 6inch, 4inch)
draw(img, p)

#<<<<<<<<<

#%% alpha 0.3 avec conf et sans conf
t = df[df[:alpha] .== 0.3,:]
dfplot = t[t[:Seuil_confidence] .== 15.0,:]
dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
M = maximum(dfplot[:Perf])
m = minimum(dfplot[:Perf])
mtot = max(M,-m)

p=plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5))
img = SVG(plotsdir("05-Batch100","Delta-Best-Uni-Conf15-alpha03.svg"), 6inch, 4inch)
draw(img, p)

t = df[df[:alpha] .== 0.3,:]
dfplot = t[t[:Seuil_confidence] .== 20.0,:]
dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
M = maximum(dfplot[:Perf])
m = minimum(dfplot[:Perf])
mtot = max(M,-m)

p=plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5))
img = SVG(plotsdir("05-Batch100","Delta-Best-Uni-Conf20-alpha03.svg"), 6inch, 4inch)
draw(img, p)

#<<<<<<<<<


#%% alpha 0.4 avec conf et sans conf
t = df[df[:alpha] .== 0.4,:]
dfplot = t[t[:Seuil_confidence] .== 15.0,:]
dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
M = maximum(dfplot[:Perf])
m = minimum(dfplot[:Perf])
mtot = max(M,-m)

p=plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5))
img = SVG(plotsdir("05-Batch100","Delta-Best-Uni-Conf15-alpha04.svg"), 6inch, 4inch)
draw(img, p)

t = df[df[:alpha] .== 0.4,:]
dfplot = t[t[:Seuil_confidence] .== 20.0,:]
dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
M = maximum(dfplot[:Perf])
m = minimum(dfplot[:Perf])
mtot = max(M,-m)

p=plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5))
img = SVG(plotsdir("05-Batch100","Delta-Best-Uni-Conf20-alpha04.svg"), 6inch, 4inch)
draw(img, p)

#<<<<<<<<<


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
