#' '# On fait les fig en pyplot vite fait et on refait dans script en plotlyjs
#'  en enregistrant dans _reseacrh
using DrWatson
quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
using BSON
using DataFrames
include(srcdir("load_data.jl"))
#' #TODO: src files with load functions for data files
path = "figures\\notebook14\\"
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
:learning_rate)#,
#:corrupted)

@everywhere function foldername_create(dummy::ListParameters_θconfidenceSimulation,box::Float64)
    seuil = dummy.seuil_decision
    Modele = "Decision_Seuil=$seuil"
    tcat = dummy.cats
    tcatnametemp = tcat.name
    alpha = tcat.width[1]
    tcatname = "$tcatnametemp-Alpha=$alpha"


    f = "$Modele\\"
    return f
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


#%% Parameters of the coding layer
Nlist = zeros(31)
for i=1:31
    Nlist[i] = 9+i
end

αlist = [i for i=10:40]
αlist = αlist./100


thlist=[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0]
percentagelist20 = threshold_to_percentage(thlist,20.0)

#<<<<<<<<<<<


#%% Parameters of simulations - Comparison between Conf modualted and rewardm odualted
general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.20,0.3,0.4], # largeur des catégories gaussiennes
    "box_conv" => 0.1,
	"learning_rate" =>[0.05],#[0.05,0.5,0.1]#[1.0,0.1]
	"batch" => [20,100]
)

dicts = dict_list(general_args)


df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
box_conv = Float64[],Perf=Float64[],Stim = Float64[],type_TC = String[])
for dict in dicts

	box = dict["box_conv"]
	indice = 20#dict["batch"]
	dictt = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => dict["type_TC"],
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\07-script\\performances\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => 20.0,#dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => dict["type_TC"],
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dict_temp_Uniform )
	fnametemp = foldername_create_accuracy(simulations)
	sname = savename(simulations,"bson")
	fname = "sims\\08-script\\P_App\\Batch=$indice\\$fnametemp"
	dataUniform = load(datadir(fname,sname))
	dataDictUniform = dataUniform[:dict]
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)

	dummy = -(dataUniform[:accuracy] .- data[:accuracy])
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		dataDict["box_conv"],dummy[i],data[:stimlist][i],dict["type_TC"]])
		push!(df,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i],dict["type_TC"]])
	end
end



df_layer = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
box_conv = Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = 100#dict["batch"]
	dictt = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Optimized",
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\08-script\\P_App\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => 20.0,#dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Uniform",
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dict_temp_Uniform )
	fnametemp = foldername_create_accuracy(simulations)
	sname = savename(simulations,"bson")
	fname = "sims\\08-script\\P_App\\Batch=$indice\\$fnametemp"
	dataUniform = load(datadir(fname,sname))
	dataDictUniform = dataUniform[:dict]
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)

	dummy = -(dataUniform[:accuracy] .- data[:accuracy])
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
		push!(df_layer,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		dataDict["box_conv"],dummy[i],data[:stimlist][i]])
		push!(df_layer,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i]])
	end
end

#<<<<<<<<<<<<


#%% Figure for Uniform tuning curves
using Gadfly
using ColorSchemes

font_panel = Theme(
	major_label_font="Arial",
	minor_label_font="Arial",
	major_label_font_size=16pt,
	minor_label_font_size=14pt
)
function plot_diff_learning(df,α::Float64,type_TC::String,seuil::Float64)
	t = df[df[:type_TC] .== type_TC,:]

	t_alpha = t[t[:alpha] .== α,:]
	dfplot = t_alpha[t_alpha[:Seuil_confidence] .== seuil,:]
	dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
	dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
	M = maximum(dfplot[:Perf])
	m = minimum(dfplot[:Perf])
	mtot = max(M,-m)
	p=Gadfly.plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
	Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
	Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5),Guide.title("Reward - Confidence"))
	return p
end

function plot_diff_learning_layer(t,α::Float64,type_TC::String,seuil::Float64)

	t_alpha = t[t[:alpha] .== α,:]
	dfplot = t_alpha[t_alpha[:Seuil_confidence] .== seuil,:]
	dfplot = dfplot[dfplot[:Stim] .< 0.5002,:]
	dfplot = dfplot[dfplot[:Stim] .> 0.1902,:]
	M = maximum(dfplot[:Perf])
	m = minimum(dfplot[:Perf])
	mtot = max(M,-m)
	p=Gadfly.plot(dfplot, x=:Stim , y=:N, color=:Perf, Geom.rectbin, font_panel,
	Scale.ContinuousColorScale(minvalue=-mtot,maxvalue=mtot,p -> get(ColorSchemes.balance, p)),
	Coord.Cartesian(ymin=10,ymax=40,xmin=0.2,xmax=0.5),Guide.title("Uniform - Opti"))
	return p
end

mkdir(plotsdir("notebooks\\notebook14"))

p = plot_diff_learning(df,0.2,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","fig1.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.2,"Optimized",20.0)
img = SVG(plotsdir("notebooks","notebook14","fig2.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.3,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","fig3.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.3,"Optimized",20.0)
img = SVG(plotsdir("notebooks","notebook14","fig4.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.4,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","fig5.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.4,"Optimized",20.0)
img = SVG(plotsdir("notebooks","notebook14","fig6.svg"), 6inch, 4inch)
draw(img, p)


#<<<<<<


#%% Confidence modulation
p = plot_diff_learning(df,0.2,"Uniform",15.0)
img = SVG(plotsdir("notebooks","notebook14","fig7.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.2,"Optimized",15.0)
img = SVG(plotsdir("notebooks","notebook14","fig8.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.3,"Uniform",15.0)
img = SVG(plotsdir("notebooks","notebook14","fig9.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.3,"Optimized",15.0)
img = SVG(plotsdir("notebooks","notebook14","fig10.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.4,"Uniform",15.0)
img = SVG(plotsdir("notebooks","notebook14","fig11.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.4,"Optimized",15.0)
img = SVG(plotsdir("notebooks","notebook14","fig12.svg"), 6inch, 4inch)
draw(img, p)

#<<<<

#%% Confidence modulation and not for RM
p = plot_diff_learning(df,0.2,"Optimized",20.0)
img = SVG(plotsdir("notebooks","notebook14","Half-conf-fig13.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.2,"Optimized",15.0)
img = SVG(plotsdir("notebooks","notebook14","Half-conf-fig14.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.3,"Optimized",20.0)
img = SVG(plotsdir("notebooks","notebook14","Half-conf-fig15.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.3,"Optimized",15.0)
img = SVG(plotsdir("notebooks","notebook14","Half-conf-fig16.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.4,"Optimized",20.0)
img = SVG(plotsdir("notebooks","notebook14","Half-conf-fig17.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning(df,0.4,"Optimized",15.0)
img = SVG(plotsdir("notebooks","notebook14","Half-conf-fig18.svg"), 6inch, 4inch)
draw(img, p)

#<<<<

#%% Diff Layer for remward modulated
p = plot_diff_learning_layer(df_layer,0.2,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Layer-fig1.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning_layer(df_layer,0.3,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Layer-fig2.svg"), 6inch, 4inch)
draw(img, p)


p = plot_diff_learning_layer(df_layer,0.4,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Layer-fig3.svg"), 6inch, 4inch)
draw(img, p)


p = plot_diff_learning_layer(df_layer,0.2,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Long-Layer-fig1.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning_layer(df_layer,0.3,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Long-Layer-fig2.svg"), 6inch, 4inch)
draw(img, p)


p = plot_diff_learning_layer(df_layer,0.4,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Long-Layer-fig3.svg"), 6inch, 4inch)
draw(img, p)

#<<<<<<<<<<<<











#%%%  Long stage learning for PHL and CML




df_PHL = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
box_conv = Float64[],Perf=Float64[],Stim = Float64[])
for dict in dicts

	box = dict["box_conv"]
	indice = 100#dict["batch"]
	dictt = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Optimized",
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\07-script\\performances\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]

	dict_temp_Uniform = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => 20.0,#dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Uniform",
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dict_temp_Uniform )
	fnametemp = foldername_create_accuracy(simulations)
	sname = savename(simulations,"bson")
	fname = "sims\\07-script\\performances\\Batch=$indice\\$fnametemp"
	dataUniform = load(datadir(fname,sname))
	dataDictUniform = dataUniform[:dict]
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)

	dummy = -(dataUniform[:accuracy] .- data[:accuracy])
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=1:length(data[:stimlist])
		push!(df_PHL,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		dataDict["box_conv"],dummy[i],data[:stimlist][i]])
		push!(df_PHL,[dataDict["N"],dataDict["α"],dataDict["Seuil_confidence"],percentagetemp,
		dataDict["box_conv"],-dummy[i],1.0-data[:stimlist][i]])
	end
end

p = plot_diff_learning_layer(df_PHL,0.2,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Long-PHL-fig1.svg"), 6inch, 4inch)
draw(img, p)

p = plot_diff_learning_layer(df_PHL,0.3,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Long-PHL-fig2.svg"), 6inch, 4inch)
draw(img, p)


p = plot_diff_learning_layer(df_PHL,0.4,"Uniform",20.0)
img = SVG(plotsdir("notebooks","notebook14","Long-PHL-fig3.svg"), 6inch, 4inch)
draw(img, p)

#<<<<<<<<<<<<



#%% Learning speed: confidence and reward modulated (probabilistic approximation)

###First let's start with confidence modulated
using Plots
plotlyjs()
general_args = Dict(
    "N" => Nlist, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => thlist, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.20,0.3,0.4], # largeur des catégories gaussiennes
    "box_conv" => 0.1,
	"learning_rate" =>[0.05],#[0.05,0.5,0.1]#[1.0,0.1]
	"batch" => [1,2,3,4,5,6,7,8,9,10,20,50,100,200]
)

dicts = dict_list(general_args)
listbatch = [1,2,3,4,5,6,7,8,9,10,20,50,100,200]




listperf = zeros(length(listbatch))
for i=1:length(listbatch)
	box = 0.1
	indice = listbatch[i]
	dictt = Dict(
	"N" => 15, # nombre de neurones dans le codage
	"T_simu" => 1000,
	"ΔTsave" => 5, # pas de temps de sauvegarde
	"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Uniform",
	"box_conv" => 0.1,
	"α"=>0.2,
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\07-script\\P_App\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
listperf[i] = data[:accuracy][95]
#print(data[:stimlist][95])
end
listperf20 = zeros(length(listbatch))
for i=1:length(listbatch)
	box = 0.1
	indice = listbatch[i]
	dictt = Dict(
	"N" => 15, # nombre de neurones dans le codage
	"T_simu" => 1000,
	"ΔTsave" => 5, # pas de temps de sauvegarde
	"Seuil_confidence" => 20.0, # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Uniform",
	"box_conv" => 0.1,
	"α"=>0.2,
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\07-script\\P_App\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
listperf20[i] = data[:accuracy][95]
#print(data[:stimlist][95])
end



plot(listbatch,listperf)
plot!(listbatch,listperf20)




listperf = zeros(length(listbatch))
for i=1:length(listbatch)
	box = 0.1
	indice = listbatch[i]
	dictt = Dict(
	"N" =>35, # nombre de neurones dans le codage
	"T_simu" => 1000,
	"ΔTsave" => 5, # pas de temps de sauvegarde
	"Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Optimized",
	"box_conv" => 0.1,
	"α"=>0.2,
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\07-script\\P_App\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
listperf[i] = data[:accuracy][75]
#print(data[:stimlist][95])
end
listperf20 = zeros(length(listbatch))
for i=1:length(listbatch)
	box = 0.1
	indice = listbatch[i]
	dictt = Dict(
	"N" => 35, # nombre de neurones dans le codage
	"T_simu" => 1000,
	"ΔTsave" => 5, # pas de temps de sauvegarde
	"Seuil_confidence" => 20.0, # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Optimized",
	"box_conv" => 0.1,
	"α"=>0.2,
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	sname = savename(simulations,"bson")
	fnametemp = foldername_create_accuracy(simulations)
	fname = "sims\\07-script\\P_App\\Batch=$indice\\$fnametemp"
	data = load(datadir(fname,sname))
	dataDict = data[:dict]
listperf20[i] = data[:accuracy][75]
#print(data[:stimlist][95])
end



plot(listbatch,listperf)
plot!(listbatch,listperf20)

# faire code 100 premiers puis par pas de 100


function generate_Fisher_learning(N,α,stimlist,weights)
    mycatGauss = MyCatGauss(width=[α,α])

    codinglayer = init_layer(N,mycatGauss)
    df = DataFrame(stim=Float64[],Delta_FCode=Float64[])

    FcodeO = zeros(length(stimlist))

    for i=1:length(stimlist)
        x = stimlist[i]
    FcodeO[i] =Fcode(x,codinglayer.centers,codinglayer.widths,weights)


end

    return FcodeO
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
	prefix = "" # prefix for savename
	)


	weights1_temp = file[1][:list_weights_to_unit1]
	weights2_temp = file[1][:list_weights_to_unit2]

	mw1 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)
	mw2 = zeros(Int(simulationL.T_simu/simulationL.ΔTsave),simulationL.N)

	for i=1:10

		for j=1:Int(simulationL.T_simu/simulationL.ΔTsave)-(Int(meanindice-1))
			for k=1:simulationL.N
				for q=0:Int(meanindice-1)
				mw1[j,k] += weights1_temp[i][j+q,k] .- weights1_temp[i][j+q,end+1-k]
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


stimlist=[i*0.005 for i=0:100]

test = generate_Fisher_learning(15,0.2,stimlist,ones(15))

plot(stimlist,test)


using PlotThemes
theme(:default)

fig = plot()
fig2 = plot()
listbatch=[2,5,10,15,18]
FcodeO = zeros(length(stimlist),length(listbatch))
N=15
McodeO = zeros(N,length(listbatch))
for i=1:length(listbatch)
	box = 0.1
	indice = listbatch[i]
	dictt = Dict(
	"N" => N, # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => 20.0, # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20.0, # seuil de la décision
	"type_TC" => "Optimized",
	"box_conv" => 0.1,
	"α"=>0.25,
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dictt)
	w1,w2,mw1,mw2=	generate_mean_weights_script07(dictt,indice,10)
	mycatGauss = MyCatGauss(width=[dictt["α"],dictt["α"]])
	#codinglayer = init_layer(dictt["N"])
	codinglayer = init_layer(dictt["N"],mycatGauss)

	for j=1:length(stimlist)
        x = stimlist[j]
    FcodeO[j,i] =Fcode(x,codinglayer.centers,codinglayer.widths,abs.(mw1))
end
McodeO[:,i] = mw1

end

plot!(fig,stimlist,FcodeO,palette=[ColorGradient(:blues)[1.0-z] for z=LinRange(0,1,length(listbatch))])
fig

plot!(fig2,McodeO,palette=[ColorGradient(:inferno)[1.0-z] for z=LinRange(0,1,length(listbatch))])
fig2

plot!(fig2,mw1)

fig

fig2

using PerceptualColourMaps
using ColorSchemes
clibraries()
ColorSchemes.ice.colors
#<<<<<<<<<




########
# test fig 2 pour RMHL

#%% Parameters of simulations - Comparison between Conf modualted and rewardm odualted
general_args = Dict(
    "N" => [10,20,30], # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => 20, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20, # seuil de la décision
    "type_TC" => ["Uniform"],
    "α" => [0.25], # largeur des catégories gaussiennes
    "box_conv" => 0.1,
	"learning_rate" =>[0.05],#[0.05,0.5,0.1]#[1.0,0.1]
	"batch" => [20]
)

dicts = dict_list(general_args)


df = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
box_conv = Float64[],Perf=Float64[],Stim = Float64[],type_TC = String[])
for dict in dicts

	box = dict["box_conv"]
	indice = 20#dict["batch"]


	dict_temp_Uniform = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => 20,#dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20, # seuil de la décision
	"type_TC" => dict["type_TC"],
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dict_temp_Uniform )
	fnametemp = foldername_create_accuracy(simulations)
	sname = savename(simulations,"bson")
	s2 = string("N=",string(dict["N"]),"_Seuil_confidence=20_T_simu=10000_box_conv=0.1_learning_rate=0.05_seuil_decision=20_type_TC=",dict["type_TC"],"_ΔTsave=100_α=0.25.bson")
	fname = "sims\\08-script\\P_App\\Batch=$indice\\$fnametemp"
	dataUniform = load(datadir(fname,s2))
	dataDictUniform = dataUniform[:dict]
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)

	dummy = dataUniform[:accuracy] 
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=30:length(dataUniform[:stimlist])
		push!(df,[dataDictUniform["N"],dataDictUniform["α"],dataDictUniform["Seuil_confidence"],percentagetemp,
		dataDictUniform["box_conv"],dummy[i],dataUniform[:stimlist][i],dict["type_TC"]])
		push!(df,[dataDictUniform["N"],dataDictUniform["α"],dataDictUniform["Seuil_confidence"],percentagetemp,
		dataDictUniform["box_conv"],1.0-dummy[i],1.0-dataUniform[:stimlist][i],dict["type_TC"]])
	end
end



df


scatter(df[!,:Stim][df[!,:N] .== 30 ],df[!,:Perf][df[!,:N] .== 30 ], color=:green)
scatter!(df[!,:Stim][df[!,:N] .== 20 ],df[!,:Perf][df[!,:N] .== 20 ], color=:green,alpha=0.5)
scatter!(df[!,:Stim][df[!,:N] .== 10 ],df[!,:Perf][df[!,:N] .== 10 ], color=:green,alpha=0.3)



df2 = DataFrame(N=Int[],alpha=Float64[],Seuil_confidence = Float64[],percentage_confidence = Float64[],
box_conv = Float64[],Perf=Float64[],Stim = Float64[],type_TC = String[])
for dict in dicts

	box = dict["box_conv"]
	indice = 20#dict["batch"]


	dict_temp_Uniform = Dict(
	"N" => dict["N"], # nombre de neurones dans le codage
	"T_simu" => 10000,
	"ΔTsave" => 100, # pas de temps de sauvegarde
	"Seuil_confidence" => 20,#dict["Seuil_confidence"], # seuil sur la confiancep our l'appr.
	"seuil_decision" => 20, # seuil de la décision
	"type_TC" => "Optimized",
	"box_conv" => dict["box_conv"],
	"α"=>dict["α"],
	"learning_rate"=>0.05)
	simulations = create_simu_from_dict(dict_temp_Uniform )
	fnametemp = foldername_create_accuracy(simulations)
	sname = savename(simulations,"bson")
	s2 = string("N=",string(dict["N"]),"_Seuil_confidence=20_T_simu=10000_box_conv=0.1_learning_rate=0.05_seuil_decision=20_type_TC=",dict_temp_Uniform["type_TC"],"_ΔTsave=100_α=0.25.bson")
	fname = "sims\\08-script\\P_App\\Batch=$indice\\$fnametemp"
	dataUniform = load(datadir(fname,s2))
	dataDictUniform = dataUniform[:dict]
	cats = MyCatGauss(width=[dict["α"],dict["α"]])
	codinglayer = init_layer(Int(dict["N"]),cats)

	dummy = dataUniform[:accuracy] 
	percentagetemp = percentagelist20[findfirst(x -> x==dict["Seuil_confidence"],thlist)]

	for i=30:length(dataUniform[:stimlist])
		push!(df2,[dataDictUniform["N"],dataDictUniform["α"],dataDictUniform["Seuil_confidence"],percentagetemp,
		dataDictUniform["box_conv"],dummy[i],dataUniform[:stimlist][i],dict["type_TC"]])
		push!(df2,[dataDictUniform["N"],dataDictUniform["α"],dataDictUniform["Seuil_confidence"],percentagetemp,
		dataDictUniform["box_conv"],1.0-dummy[i],1.0-dataUniform[:stimlist][i],dict["type_TC"]])
	end
end

sort!(df,[:Stim])
fig=plot()
plot!(fig,df[!,:Stim][df[!,:N] .== 30 ],df[!,:Perf][df[!,:N] .== 30 ], color=:green, w=4, label = "Uniform coding 30")
plot!(fig,df[!,:Stim][df[!,:N] .== 20 ],df[!,:Perf][df[!,:N] .== 20 ], color=:green,alpha=0.5, w=4, label = "Uniform coding 20")
plot!(fig,df[!,:Stim][df[!,:N] .== 10 ],df[!,:Perf][df[!,:N] .== 10 ], color=:green,alpha=0.3, w=4, label = "Uniform coding 10")

sort!(df2,[:Stim])
plot!(fig,df2[!,:Stim][df2[!,:N] .== 30 ],df2[!,:Perf][df2[!,:N] .== 30 ], color=:gray, w=4, label = "Optimized coding 30")
plot!(fig,df2[!,:Stim][df2[!,:N] .== 20 ],df2[!,:Perf][df2[!,:N] .== 20 ], color=:gray,alpha=0.5, w=4, label = "Optimized coding 20")
plot!(fig,df2[!,:Stim][df2[!,:N] .== 10 ],df2[!,:Perf][df2[!,:N] .== 10 ], color=:gray,alpha=0.3, w=4, label = "Optimized coding 10")


savefig(plotsdir("notebooks","notebook14","RMHL-performances.svg"))


