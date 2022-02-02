#%% Preambule
using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using Distributed
using BSON
using DataFrames
using CSV
include(srcdir("training_network.jl"))

using Plots
pyplot()
#<<<<<<<<<<<<<

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

function foldername_create(dummy::ListParameters_θconfidenceSimulation)
    seuil = dummy.seuil_decision
    Modele = "Decision_Seuil=$seuil"
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = dummy.α # this is because alpha is truncated when saved
    tcatname = "$tcatnametemp-Alpha=$alpha"


    f = "$Modele\\$tcatname\\"
    return f
end
#<<<<<<<<<<<<<<<<<<<<<<<<<<<


#%% Functions used for prouce or load the data
function create_simu_from_dict(diclist)
	@unpack N,T_simu,ΔTsave,Seuil_confidence,seuil_decision,type_TC,box_conv,α = diclist

    mycatGauss = MyCatGauss(width=[0.257,0.257])
    simulations = ListParameters_θconfidenceSimulation(N=N,
    T_simu=T_simu,
    ΔTsave=ΔTsave,
    Seuil_confidence=Seuil_confidence,
    seuil_decision=seuil_decision,
    type_TC=type_TC,
    α = 0.257,
    cats=mycatGauss,
	box_conv = box_conv)

    return simulations
end

function create_simu_Weibull_from_dict(diclist)
	@unpack N,T_simu,ΔTsave,Seuil_confidence,seuil_decision,type_TC,box_conv,α = diclist


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
    α = 0.257,#alphag, Approximation
    cats=mycatWeib,
	box_conv=box_conv)

    return simulations
end
function mean_learning_Gaussian_Cat(diclist,box=1.0,rep=10)
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

function save_meanlearning_Gauss(diclist)
    box = diclist["box_conv"]

    dummy = mean_learning_Gaussian_Cat(diclist,box)
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


#<<<<<<<<<<<<<


#%%% Functions to obtain the mean_weights

function generate_mean_weights(simulation)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL)
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

#<<<<<<<<<


#%% Compute Performances network
function true_performance(simulationstemp,stimlist,tempnumber = 1000)
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
		mw1,mw2,d1,d2 = generate_mean_weights(simulationstemp)
		codinglayer.weightsto1 = copy(d1)
		codinglayer.weightsto2 = copy(d2)

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

#<<<<<<<<<<<

#%%%

Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end


thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]

general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Hebbian","Uniform","Optimized"],
    "box_conv" => [1.0,0.1],
	"α"=>0.257)


    general_args = Dict(
        "N" => 20, # nombre de neurones dans le codage
        "T_simu" => 10000,
        "ΔTsave" => 100, # pas de temps de sauvegarde
        "Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
        "seuil_decision" => 20.0, # seuil de la décision
        "type_TC" => "Uniform",
        "box_conv" => 0.1,
        "α"=>0.257)
        simulations = create_simu_from_dict(general_args)

        mw1,mw2,a,b = generate_mean_weights(simulations)

        using VegaLite
        
        df_temp = DataFrame(Poids=Float64[], N=Int64[], epoch = Int64[],Trials=Int64[])
        
        for i = 1:length(mw1[:,1])
        for j=1:length(mw1[end,:])
        
            push!(df_temp,(mw1[i,j],j,i,i*100))
        
        end
        end
        df_temp
        
        df_temp |> @vlplot( width=800,
        height=500,
        :line,
        x={:N, axis= {tickCount = 10}},
        y ={:Poids, axis= {tickCount = 10}}, 
        color={:Trials,scale={scheme="plasma",reverse=true}})
        
        
dicts = dict_list(general_args)

dicts[1]["N"]

length(dicts)
#<<<<<<<<<


#%% Test analyse performances
function save_data(dicts, stimlist)
	df = DataFrame(stim=Float64[],Perf=Float64[],N = Int[], Seuil_confidence = Float64[], type_TC = String[], box_conv = Float64[] )

	for d in dicts
	a=true_performance(d,stimlist,2000)
	@simd for i=1:length(stimlist)
		push!(df,[stimlist[i],a[i],d["N"],d["Seuil_confidence"],d["type_TC"],d["box_conv"]])
	end
	end

	return df
end

stimlist = [0.3 + 0.005*i for i=0:40]

UniformTC=[dicts[i]["type_TC"] .== "Uniform" for i=1:length(dicts)]
dicts_Uniform = dicts[UniformTC]
@time df = save_data(dicts_Uniform,stimlist)
CSV.write(datadir("notebooks","04-dataUniform-20Hz.csv"),df)

OptimizedTC=[dicts[i]["type_TC"] .== "Optimized" for i=1:length(dicts)]
dicts_Optimized = dicts[OptimizedTC]
@time df = save_data(dicts_Optimized,stimlist)
CSV.write(datadir("notebooks","04-dataOptimized-20Hz.csv"),df)


HebbianTC=[dicts[i]["type_TC"] .== "Hebbian" for i=1:length(dicts)]
dicts_Hebbian = dicts[HebbianTC]
@time df = save_data(dicts_Hebbian,stimlist)
CSV.write(datadir("notebooks","04-dataHebbian-20Hz.csv"),df)
#<<<<<<<
