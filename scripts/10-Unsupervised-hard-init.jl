# Here I run the simulations for unsupervised learning
# there will be 2 unsupervised (random or with structure)
# influence of coding layer
# firs try will be a structure (symmetry) but nothing else



using DrWatson
quickactivate("C:\\Users\\kevin\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
#ubuntu load
#quickactivate("/home/kevin/Documents/Boulot/1 - PhD Projects/4-project_coding_categories/","4-project_coding_categories")

# now can load Packages
using Distributed
using BSON
using DataFrames

#addprocs(1)
@everywhere using DrWatson
@everywhere quickactivate("C:\\Users\\kevin\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
#ubuntu load
#@everywhere quickactivate("/home/kevin/Documents/Boulot/1 - PhD Projects/4-project_coding_categories/","4-project_coding_categories")

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
:box_conv,
:learning_rate)

@everywhere function foldername_create(dummy::ListParameters_θconfidenceSimulation,box::Float64)
    seuil = dummy.seuil_decision
    Modele = "Decision_Seuil=$seuil"
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    f = "$Modele/"
    return f
end

@everywhere function foldername_create_accuracy(dummy::ListParameters_θconfidenceSimulation)
    seuil = dummy.seuil_decision
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = dummy.α#tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    fname = "$tcatname/"
    return fname
end
#<<<<<<<<<<<<<<<<<<<<<<<<<<<

######
# Function to run the simulations
############
@everywhere function create_simu_from_dict(diclist)
	@unpack N,T_simu,ΔTsave,Seuil_confidence,seuil_decision,type_TC,α,box_conv,learning_rate = diclist

    mycatGauss = MyCatGauss(width=[α,α])
    simulations = ListParameters_θconfidenceSimulation(N=N,
    T_simu=T_simu,
    ΔTsave=ΔTsave,
    Seuil_confidence=Seuil_confidence,
    seuil_decision=seuil_decision,
    type_TC=type_TC,
    α = α,
    cats=mycatGauss,
	box_conv = box_conv,learning_rate=learning_rate)

    return simulations
end

@everywhere function mean_learning_Gaussian_Cat(diclist,box=1.0,rep=10)
    list_weights_to_unit1 = [] # chaque elt correspond à une simu
    list_weights_to_unit2 = [] # chaque elt correspond à une simu
	box = diclist["box_conv"]
    simulations = create_simu_from_dict(diclist)

    for i=1:rep
        Weightsto1,Weightsto2 = learning_poisson_input_unsupervised_hard_init(simulations,box)
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
    fname = "sims\\10-script\\results\\$fnametemp"
    fnameparam = "sims\\10-script\\param\\"

    mkpath(datadir(fname))
    mkpath(datadir(fnameparam))

    save(datadir(fname,sname),dummy)
    save(datadir(fnameparam,sname),diclist)
    return true
end

@everywhere function generate_mean_weights_script10(simulation,indice)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL,box)
	fname = "sims\\10-script\\results\\$fnametemp"
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


@everywhere function performance_approx(simulationstemp,
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
				x = stimlist[k]
		dummy = 0
		for j=1:length(codinglayer.centers)
			dummy = dummy .+codinglayer.weightsto1[j].*fi(x,codinglayer.centers[j],codinglayer.widths[j])
		end

		beta = 1.35#1.25
		alpha_rescaled = 0.049#0.072
		dummy = dummy[1]
			if sign(dummy)<0
				perf_list[k] =0.5*exp(-exp(beta.*log.(-dummy/alpha_rescaled)))
			else
				perf_list[k] =(1.0 - 0.5*exp(-exp(beta.*log.(dummy/alpha_rescaled))))
			end

end

			return perf_list


end



@everywhere function save_accuracy_from_script09_batch(dicttemp,
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
				    "box_conv" => 0.1,
					"learning_rate"=>dicttemp["learning_rate"]
				)
		w1,w2,mw1,mw2=	generate_mean_weights_script09(dict,indice)
		#accuracy=true_performance(dict,stimlist,mw1,mw2,number)
		accuracy=performance_approx(dict,stimlist,mw1,mw2,number)

	    box = dict["box_conv"]
		simulations = create_simu_from_dict(dict)
 	    sname = savename(simulations,"bson")
	    fnametemp = foldername_create_accuracy(simulations)
	    fname = "sims\\10-script\\P_App\\Batch=$indice\\$fnametemp" # P_App means performance_approx
	    #fnameparam = "sims\\05-script\\param\\"

	    mkpath(datadir(fname))
	    #mkpath(datadir(fnameparam))

		data = @dict stimlist accuracy dict number
	    save(datadir(fname,sname),data)
	    #save(datadir(fnameparam,sname),diclist)
end


## We have a defualt learning rate of 0.005
# a default normlisation cst of 0.5

# initialize lsit of number of neurones
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

Nlist = zeros(7)

for j=1:7
   
    Nlist[j] = 9+(j-1)*5+1
end




αlist = [i for i=10:40]
αlist = αlist./100


thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
thlist=[10.0,12.0,14.0,16.0,18.0,20.0]

general_args = Dict(
    "N" => [25,30,35,40], # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" =>thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.20,0.25,0.3,0.4], # largeur des catégories gaussiennes
    "box_conv" => 0.1,
	"learning_rate" =>[0.05]
)

dicts = dict_list(general_args)
 @time save_meanlearning_Gauss(dicts[1])
@time pmap(save_meanlearning_Gauss,dicts)

general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 1000,
    "ΔTsave" => 5, # pas de temps de sauvegarde
    "Seuil_confidence" =>thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.20,0.25,0.3,0.4], # largeur des catégories gaussiennes
    "box_conv" => 0.1,
	"learning_rate" =>[0.05]
)

dicts = dict_list(general_args)
@time save_meanlearning_Gauss(dicts[1])
@time pmap(save_meanlearning_Gauss,dicts)


general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.20,0.25,0.3,0.4], # largeur des catégories gaussiennes
    "box_conv" => 0.1,
	"learning_rate" =>[0.05],#[0.05,0.5,0.1]#[1.0,0.1]
	"batch" => [10,20,100]
)

dicts = dict_list(general_args)


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
