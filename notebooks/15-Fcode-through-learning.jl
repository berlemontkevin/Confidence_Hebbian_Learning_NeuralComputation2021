#' '# On fait les fig en pyplot vite fait et on refait dans script en plotlyjs
#'  en enregistrant dans _reseacrh
using DrWatson
quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
using BSON
using DataFrames
include(srcdir("load_data.jl"))
include(srcdir("training_network.jl"))
include(srcdir("information-theory.jl"))

#' #TODO: src files with load functions for data files
path = "figures\\notebook15\\"
using Plots
using ColorSchemes
using Colors
plotlyjs()
# This notebook will study the effect of learning on the Fisher inforation

#%% Functions for Gaussian Categories
function create_simu_from_dict(diclist)
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
function generate_mean_weights_script07(simulation,indice,meanindice)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL,box)
	fname = "sims\\07-script\\results\\$fnametemp"
	file = produce_or_load(
	datadir("$fname"), # path
	simulation, # container
	mean_learning_Gaussian_Cat, # function
	prefix = "", # prefix for savename
	suffix = "bson")


	weights1_temp = file[1][:list_weights_to_unit1]
	weights2_temp = file[1][:list_weights_to_unit2]

	mw1 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)
	mw2 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)

	for i=1:10

		for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)-(Int(meanindice-1))
			for k=1:simulationL.N
				for q=0:Int(meanindice-1)
			#	mw1[j,k] += max(0.0,weights1_temp[i][j+q,k]) .- min(0.0,weights1_temp[i][j+q,end+1-k])
			#	mw2[j,k] += min(0.0,weights2_temp[i][j+q,k]) .- max(0.0,weights2_temp[i][j+q,end+1-k])
				mw1[j,k] +=weights1_temp[i][j+q,k] .- weights1_temp[i][j+q,end+1-k]
				mw2[j,k] += weights2_temp[i][j+q,k] .- weights2_temp[i][j+q,end+1-k]

			end
		end
		end

	end
	mw1 = mw1 .- mw2
	mw2 = -mw1
	mw1 = mw1./(40*meanindice)
	mw2 = mw2./(40*meanindice)

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

function mean_learning_Gaussian_Cat_Engel(diclist,box=1.0,rep=10)
    list_weights_to_unit1 = [] # chaque elt correspond à une simu
    list_weights_to_unit2 = [] # chaque elt correspond à une simu
	box = diclist["box_conv"]
    simulations = create_simu_from_dict(diclist)

    for i=1:rep
        Weightsto1,Weightsto2 = learning_poisson_input_reward_modulated(simulations,box)
        push!(list_weights_to_unit1,Weightsto1)
        push!(list_weights_to_unit2,Weightsto2)
        Weightsto1 = nothing
        Weightsto2 = nothing
    end



    data = @dict list_weights_to_unit1 list_weights_to_unit2 simulations

    #save(datadir(fname,sname),data)
    return data
end

function generate_mean_weights_script_layer(simulation,indice,meanindice)
	# generate the mean wieghts after learning
	#
	simulationL = create_simu_from_dict(simulation)
	box = simulation["box_conv"]
	norm_w = simulationL.cst_normalisation

	fnametemp = foldername_create(simulationL,box)
	fname = "sims\\08-script\\results\\$fnametemp"


	file = produce_or_load(
		datadir("$fname"), # path
		simulation, # container
		mean_learning_Gaussian_Cat_Engel, # function
		prefix = "", # prefix for savename
		suffix="bson")


	weights1_temp = file[1][:list_weights_to_unit1]
	weights2_temp = file[1][:list_weights_to_unit2]

	mw1 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)
	mw2 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)

	for i=1:10

		for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)-(Int(meanindice-1))
			for k=1:simulationL.N
				for q=0:Int(meanindice-1)
			#	mw1[j,k] += max(0.0,weights1_temp[i][j+q,k]) .- min(0.0,weights1_temp[i][j+q,end+1-k])
			#	mw2[j,k] += min(0.0,weights2_temp[i][j+q,k]) .- max(0.0,weights2_temp[i][j+q,end+1-k])
				mw1[j,k] +=weights1_temp[i][j+q,k] .- weights1_temp[i][j+q,end+1-k]
				mw2[j,k] += weights2_temp[i][j+q,k] .- weights2_temp[i][j+q,end+1-k]

			end
		end
		end

	end
	mw1 = mw1 .- mw2
	mw2 = -mw1
	mw1 = mw1./(40*meanindice)
	mw2 = mw2./(40*meanindice)

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

function generate_Fisher(N,α,stimlist,weights,type_coding)
    #this function generate the Fcode of the coding layer with a specific gain

    mycatGauss = MyCatGauss(width=[α,α])

    if type_coding == "Uniform"
        codinglayer = init_layer(N)

    else
        codinglayer = init_layer(N,mycatGauss)

    end
    FcodeL = zeros(length(stimlist))

    for j=1:length(stimlist)
        x = stimlist[j]
        FcodeL[j] =Fcode(x,codinglayer.centers,codinglayer.widths,abs.(weights))
        end


    return FcodeL
end


function Fisher_from_dict(dict,indice,stimlist)
    w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,10)

	FcodeL = generate_Fisher(dict["N"],dict["α"],stimlist,mw1,dict["type_TC"])

    return FcodeL
end

function Fisher_from_dict_layer(dict,indice,stimlist)
    w1,w2,mw1,mw2=	generate_mean_weights_script_layer(dict,indice,10)

	FcodeL = generate_Fisher(dict["N"],dict["α"],stimlist,mw1,dict["type_TC"])

    return FcodeL
end


function plot_Fisher(dict,listbatch,stimlist)
	# this function plot Fisher information for both coding layer
	fig = plot(linewidth=2)
	dict["type_TC"] = "Optimized"
	dict["Seuil_confidence"] = 15



	FcodeO = zeros(length(stimlist),length(listbatch))
	for i=1:length(listbatch)
		indice = listbatch[i]
	FcodeL= Fisher_from_dict(dict,indice,stimlist)
	FcodeO[:,i] .= FcodeL

	end

	plot!(fig,stimlist,FcodeO,line_z=(0:6)',color=:blues,colorbar=false,linewidth=2,
	label=string.(listbatch.*dict["ΔTsave"]))

	FcodeL_U = generate_Fisher_Uniform(dict["N"],dict["α"],stimlist,dict["type_TC"])
 #	plot!(fig,stimlist,FcodeL_U,color=:black,linewidth=2,
 #	label="")


	dict["type_TC"] = "Optimized"
	dict["Seuil_confidence"] = 20

	FcodeO = zeros(length(stimlist),length(listbatch))
	for i=1:length(listbatch)
		indice = listbatch[i]
	FcodeL= Fisher_from_dict(dict,indice,stimlist)
	FcodeO[:,i] .= FcodeL
	#push!(labellist," Batch $indice")

	end
	plot!(fig,stimlist,FcodeO,line_z=(0:6)',color=:blues,colorbar=false,linewidth=2,linestyle=:dash,label="")




	dict["type_TC"] = "Optimized"
	dict["Seuil_confidence"] = 20

	FcodeO = zeros(length(stimlist),length(listbatch))
	for i=1:length(listbatch)
		indice = listbatch[i]
	FcodeL= Fisher_from_dict_layer(dict,indice,stimlist)
	FcodeO[:,i] .= FcodeL
	#push!(labellist," Batch $indice")

	end

	plot!(fig,stimlist,FcodeO,line_z=(0:6)',color=:greens,colorbar=false,linewidth=2,label="")


#	plot!(fig,stimlist,FcodeO,palette=[ColorGradient(:blues)[z] for z=LinRange(0,1,length(listbatch))],linewidth=2,
#	linestyle=:dash,label="")



	labellist = String[]
	dict["type_TC"] = "Uniform"
	dict["Seuil_confidence"] = 15

	FcodeO = zeros(length(stimlist),length(listbatch))
	for i=1:length(listbatch)
		indice = listbatch[i]
	FcodeL= Fisher_from_dict(dict,indice,stimlist)
	FcodeO[:,i] .= FcodeL

	end
	plot!(fig,stimlist,FcodeO,line_z=(0:6)',color=:reds,colorbar=false,linewidth=2,
	label="")

	labellist = String[]
	dict["type_TC"] = "Uniform"
	dict["Seuil_confidence"] = 20
	FcodeO = zeros(length(stimlist),length(listbatch))
	for i=1:length(listbatch)
		indice = listbatch[i]
	FcodeL= Fisher_from_dict(dict,indice,stimlist)
	FcodeO[:,i] .= FcodeL
	push!(labellist," Batch $indice")
	end
	plot!(fig,stimlist,FcodeO,line_z=(0:6)',color=:reds,colorbar=false,linewidth=2,
	linestyle=:dash,label="")

	FcodeL_U = generate_Fisher_Uniform(dict["N"],dict["α"],stimlist,dict["type_TC"])
	plot!(fig,stimlist,FcodeL_U,color=:black,linewidth=2,
	label="")

	return fig
end


function generate_Fisher_Uniform(N,α,stimlist,type_coding)
    #this function generate the Fcode of the coding layer with a specific gain

    mycatGauss = MyCatGauss(width=[α,α])

    if type_coding == "Uniform"
        codinglayer = init_layer(N)

    else
        codinglayer = init_layer(N,mycatGauss)

    end
    FcodeL = zeros(length(stimlist))

    for j=1:length(stimlist)
        x = stimlist[j]
        FcodeL[j] =Fcode(x,codinglayer.centers,codinglayer.widths,0.25./sqrt(N).*ones(N))
        end


    return FcodeL
end

#%% Test functions
dict = Dict(
"N" => 20, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.2,
"learning_rate"=>0.05)
listbatch=[20 30 40 50 100 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)

fig

#test of uniform TC
a,b,mw1,d = generate_mean_weights_script07(dict,20,5)
sum(abs.(mw1).*abs.(mw1))
sqrt(15)*0.25
sum(abs.(mw1))

plot(mw1)

plot!(0.25./sqrt(35).*ones(35))

#<<<<<<<<<


#%% Some values of N

#%% N35
dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
listbatch=[20 30 40 50 100 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=35 alpha=0.25")
fig

savefig(fig,plotsdir("notebooks","notebook15","RM-First-steps-N35-Alpha25.svg"))



dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
listbatch=[10 20 30 40 50 90]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=35 alpha=0.25")
fig

savefig(fig,plotsdir("notebooks","notebook15","RM-Long-steps-N35-Alpha25.svg"))

#<<<<<<<<<

#%% N15
dict = Dict(
"N" => 15, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
listbatch=[20 30 40 50 100 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=15 alpha=0.25")
fig

savefig(fig,plotsdir("notebooks","notebook15","RM-First-steps-N15-Alpha25.svg"))



dict = Dict(
"N" => 15, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
listbatch=[10 20 30 40 50 90]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=15 alpha=0.25")

fig
savefig(fig,plotsdir("notebooks","notebook15","RM-Long-steps-N15-Alpha25.svg"))


#<<<<<<<<<


#%% N15
dict = Dict(
"N" => 15, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.4,
"learning_rate"=>0.05)
listbatch=[20 30 40 50 100 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=15 alpha=0.4")

fig
savefig(fig,plotsdir("notebooks","notebook15","RM-First-steps-N15-Alpha4"))




dict = Dict(
"N" => 15, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.4,
"learning_rate"=>0.05)
listbatch=[10 20 30 40 50 90]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=15 alpha=0.4")

fig
savefig(fig,plotsdir("notebooks","notebook15","RM-Long-steps-N15-Alpha4"))



#<<<<<<<<<


#%% N35
dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.4,
"learning_rate"=>0.05)
listbatch=[20 30 40 50 100 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=35 alpha=0.4")

fig
savefig(fig,plotsdir("notebooks","notebook15","RM-First-steps-N35-Alpha4"))



dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.4,
"learning_rate"=>0.05)
listbatch=[10 20 30 40 50 90]
stimlist=[i*0.005 for i=0:100]
fig = plot_Fisher(dict,listbatch,stimlist)
title!("Fcode N=35 alpha=0.4")

fig
savefig(fig,plotsdir("notebooks","notebook15","RM-Long-steps-N35-Alpha4"))



#<<<<<<<<<


#<<<<<<<<<<


#%% Performances totales du réseau au ocurs de l'apprentissage


function performance_approx(simulationstemp,
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


pyplot()

dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
stimlist=[i*0.01 for i=40:49]

listperf = Float64[]

for indice=1:160
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
test = performance_approx(dict,stimlist,mw1,mw2)
push!(listperf,mean(test))
end


plot(listperf,color=:blue,linewidth=2)

dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 20.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
stimlist=[i*0.01 for i=40:49]

listperf = Float64[]

for indice=1:160
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
test = performance_approx(dict,stimlist,mw1,mw2)
push!(listperf,mean(test))
end


plot!(listperf,linewidth=2,color=:red)
plot!(legend=false)
yaxis!((0.5,1.0))


dict = Dict(
"N" => 25, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.2,
"learning_rate"=>0.05)
listbatch=LinRange(0,1000,200)
stimlist=[i*0.01 for i=40:41]
indice = listbatch[1]

listperf = Float64[]

for indice=1:80
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,10)
test = performance_approx(dict,stimlist,mw1,mw2)
push!(listperf,mean(test))
end

plot(listperf)


#<<<<<<<<<


#%% Figure évolution perf au cours du temps

dict = Dict(
"N" => 25, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 20, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
#stimlist=[i*0.01 for i=40:49]
listtotal=[]
listperf = Float64[]
for indice=1:5:20
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
listperf=Float64[]
for i=10:2:90
	stimlist=[i*0.01]
test = performance_approx(dict,stimlist,mw1,mw2)
push!(listperf,mean(test))
end
push!(listtotal,listperf)
end

listbatch = 1:5:20.0
liststim = 10:2:90
plot(liststim.*0.01,listtotal,palette=[ColorGradient(:blues)[z] for z=LinRange(0,1,4)],linewidth=2,
label=string.(Int.(listbatch'.*dict["ΔTsave"])))
savefig(plotsdir("notebooks","notebook15","Performance-N35-alpha025.svg"))




dict = Dict(
"N" => 20, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Uniform",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
stimlist=[i*0.01 for i=40:49]
listtotal=[]
listperf = Float64[]
for indice=1:5:100
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
listperf=Float64[]
for i=20:2:80
	stimlist=i*0.01
codinglayer = init_layer(dict["N"])
codinglayer.weightsto1 = mw1
codinglayer.weightsto2 = mw2

test = performance_approx(stimlist,codinglayer)
push!(listperf,mean(test))
end
push!(listtotal,listperf)
end

plotlyjs()
listbatch = 1:5:20.0
liststim = 20:2:80
fig = plot()
plot!(fig,liststim.*0.01,listtotal,linewidth=2,line_z=(100:1000:10000)',color=:blues,lw=3,linealpha=0.7)
plot!(legend=false)
fig
savefig(plotsdir("notebooks","notebook15","Performance-Uniform-N20-alpha025.svg"))



dict = Dict(
"N" => 20, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
stimlist=[i*0.01 for i=40:49]
listtotal=[]
listperf = Float64[]
for indice=1:5:100
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
listperf=Float64[]
for i=20:1:80
	stimlist=i*0.01
codinglayer = init_layer(dict["N"], MyCatGauss(width=[dict["α"],dict["α"]]))
codinglayer.weightsto1 = mw1
codinglayer.weightsto2 = mw2

test = performance_approx(stimlist,codinglayer)
push!(listperf,mean(test))
end
push!(listtotal,listperf)
end

plotlyjs()
pyplot()
listbatch = 1:5:20.0
liststim = 20:1:80
fig = plot()
plot!(fig,liststim.*0.01,listtotal,linewidth=2,line_z=(100:1000:10000)',color=:blues,lw=3,linealpha=0.7)
plot!(legend=false)
fig
savefig(plotsdir("notebooks","notebook15","Performance-Opti-N20-alpha025.svg"))


#<<<<<<<<<


#%%%%%%%%%% weights through batch


dict = Dict(
"N" => 20, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
stimlist=[i*0.01 for i=40:49]
listtotal=[]
listperf = Float64[]
for indice=1:5:100
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
listperf=Float64[]
#for i=20:1:80
#stimlist=i*0.01
codinglayer = init_layer(dict["N"], MyCatGauss(width=[dict["α"],dict["α"]]))
codinglayer.weightsto1 = mw1
codinglayer.weightsto2 = mw2

#test = performance_approx(stimlist,codinglayer)
#push!(listperf,mean(test))
#end
push!(listtotal,mw1)
end

plotlyjs()
pyplot()
listbatch = 1:5:20.0
liststim = 20:1:80
fig = plot()
codinglayer = init_layer(dict["N"], MyCatGauss(width=[dict["α"],dict["α"]]))

plot!(fig,codinglayer.centers,listtotal,linewidth=2,line_z=(100:1000:10000)',color=:blues,lw=3,linealpha=0.7)
plot!(legend=false)
fig
savefig(plotsdir("notebooks","notebook15","Weights-Opti-N20-alpha025.svg"))




dict = Dict(
"N" => 20, # nombre de neurones dans le codage
"T_simu" => 10000,
"ΔTsave" => 100, # pas de temps de sauvegarde
"Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Uniform",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
stimlist=[i*0.01 for i=40:49]
listtotal=[]
listperf = Float64[]
for indice=1:5:100
w1,w2,mw1,mw2=	generate_mean_weights_script07(dict,indice,1)
listperf=Float64[]
#for i=20:1:80
#stimlist=i*0.01
codinglayer = init_layer(dict["N"])
codinglayer.weightsto1 = mw1
codinglayer.weightsto2 = mw2

#test = performance_approx(stimlist,codinglayer)
#push!(listperf,mean(test))
#end
push!(listtotal,mw1)
end

plotlyjs()
pyplot()
listbatch = 1:5:20.0
liststim = 20:1:80
fig = plot()
codinglayer = init_layer(dict["N"], MyCatGauss(width=[dict["α"],dict["α"]]))

plot!(fig,codinglayer.centers,listtotal,linewidth=2,line_z=(100:1000:10000)',color=:blues,lw=3,linealpha=0.7)
plot!(legend=false)
fig
savefig(plotsdir("notebooks","notebook15","Weights-Uniform-N20-alpha025.svg"))







#%%% early steps for fisher information

#%% N35
dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
listbatch=[2 5 10 15 20]# 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
#title!("Fcode N=35 alpha=0.25")
savefig(plotsdir("notebooks","notebook15","Fisher-early-N35-alpha025.svg"))

fig


dict = Dict(
"N" => 15, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)
listbatch=[2 5 10 15 20]# 150 180]
stimlist=[i*0.005 for i=0:100]
pyplot()
fig = plot_Fisher(dict,listbatch,stimlist)
#title!("Fcode N=35 alpha=0.25")
savefig(plotsdir("notebooks","notebook15","Fisher-early-N15-alpha025.svg"))

fig



##%%% generate TUning curves
N=20
codinglayer = init_layer(N)

dict = Dict(
"N" => 35, # nombre de neurones dans le codage
"T_simu" => 1000,
"ΔTsave" => 5, # pas de temps de sauvegarde
"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
"seuil_decision" => 20.0, # seuil de la décision
"type_TC" => "Optimized",
"box_conv" => 0.1,
"α"=>0.25,
"learning_rate"=>0.05)


codinglayer = init_layer(N,mycatGauss)

#<<<<<<<<<<<<<<<<<
