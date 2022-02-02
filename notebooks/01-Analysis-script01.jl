using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using BSON
using DataFrames
using CSV
include(srcdir("training_network.jl"))
#using Makie
#using AbstractPlotting

# In this file: various parameters are tested
# We save various data TODO

### Important
# We do not change the value of α for the Gaussian ctaegorie
# It stays constant through the learning
# We save the data for Uniform, Hebbian and Optimized tuning curves

#########



#%% Plots for temp fig
using Plots
pyplot()
#<<<<<<<<<<<

############
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
:α)

function foldername_create(dummy::ListParameters_θconfidenceSimulation)
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

#%% Function to save data (useful forp roduce or load)

######
# Function to run the simulations
############
 function create_simu_from_dict(diclist)
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

function mean_learning_Gaussian_Cat(diclist,rep=10)
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


 function save_meanlearning_Gauss(diclist)

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


#<<<<<<<<<<<
function generate_mean_weights(simulation)
    # generate the mean wieghts after learning
    #
    fnametemp = foldername_create(simulation)
    fname = "sims\\01-script\\results\\$fnametemp"
    file = produce_or_load(
         datadir("$fname"), # path
         simulation, # container
         mean_learning_Gaussian_Cat, # function
         prefix = "" # prefix for savename
      )


    weights1_temp = file[1][:list_weights_to_unit1]
    weights2_temp = file[1][:list_weights_to_unit2]

    mw1 = zeros(100,simulation.N)
    mw2 = zeros(100,simulation.N)

    for i=1:10
        mw1 = weights1_temp[i] .- weights1_temp[i][:,end:-1:1]
        mw2 = weights2_temp[i] .- weights2_temp[i][:,end:-1:1]
        mw1 = mw1 .- mw2
        mw2 = -mw1

    end
    mw1 = mw1./40
    mw2 = mw2./40

    return mw1,mw2
end



simulationstemp = Dict(
    "N" => 20, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => 15.0, # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Uniform",
    "α" => 0.25 # largeur des catégories gaussiennes
)
simulations = create_simu_from_dict(simulationstemp)


mw1,mw2 = generate_mean_weights(simulations)

using VegaLite

df_temp = DataFrame(Poids=Float64[], N=Int64[], epoch = Int64[],Trials=Int64[])

for i = 10:length(mw1[:,1])
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






plot(mw1[end,:])


#%% key functions to analyze the data


function performance(domaine,centers,widths,gain)
	coherencelist = zeros(length(domaine))
	perflist = zeros(length(coherencelist))
	gain2 = copy(gain)
	for i=1:length(domaine)
		dummy = 0
		for j=1:length(centers)
			dummy = dummy .+gain2[j].*fi(domaine[i],centers[j],widths[j])
		end
		coherencelist[i] = dummy[1]

	end
	beta = 1.29 # this value is for ??
	alpha_rescaled = 0.05 ##
	for i = 1:length(coherencelist)
		if sign(coherencelist[i])<0
			perflist[i] =0.5*exp(-exp(beta.*log.(-coherencelist[i]/alpha_rescaled)))
		else
			perflist[i] =(1.0 - 0.5*exp(-exp(beta.*log.(coherencelist[i]/alpha_rescaled))))
		end

	end
	return perflist,coherencelist
end

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


percentagelist20 = threshold_to_percentage([10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0],20.0)



#<<<<<<<

#%%% Obtain performances at fixed stim

function fixed_stim(stimlist,diclist)
    # diclist is the list of dictionnary we are interested in
    # we save everything inside a dataFrame

    ## the dataframe has a field perf, theta,N, type TCs
	percentagelist20 = threshold_to_percentage([10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0],20.0)
	percentagelist15 = threshold_to_percentage([5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0],15.0)

	df = DataFrame([Float64,Float64,Int,String,Float64,Float64,Float64],[:Delta_Perf,:Theta,:N,:Type_TC,:Alpha,:Percentage,:Stim])
	for stim in stimlist
	for d in diclist
	N= d["N"]
	simulations = create_simu_from_dict(d)
	mw1,mw2 = generate_mean_weights(simulations)
	dtemp = copy(d)
	dtemp["Seuil_confidence"]=dtemp["seuil_decision"]
	simulations = create_simu_from_dict(dtemp)
	mw1SC,mw2SC = generate_mean_weights(simulations)

	if d["type_TC"]=="Uniform"
		layer_G = init_layer(N)
		layer_U = init_layer(N)
	else
		layer_G = init_layer(N,simulations.cats)
		layer_U = init_layer(N)
	end


Δperf =  performance(stim,layer_G.centers,layer_G.widths,(mw1[end,:]))[1] - performance(stim,layer_G.centers,layer_G.widths,(mw1SC[end,:]))[1]
percentagetemp = percentagelist20[findfirst(x -> x==d["Seuil_confidence"],[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0])]
#percentagetemp = percentagelist15[findfirst(x -> x==d["Seuil_confidence"],[5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0])]
push!(df,[Δperf,d["Seuil_confidence"],d["N"],d["type_TC"],d["α"],percentagetemp,stim])

end
end

	return df
end


function full_data(diclist)
	 # return a dataframe
	 # performance, typeTC, stim,N,seuil,percentage,alpha
	 percentagelist20 = threshold_to_percentage([10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0],20.0)
	 percentagelist15 = threshold_to_percentage([5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0],15.0)

	 df = DataFrame([Float64,Float64,Float64,Int,String,Float64,Float64,Float64],[:Delta_Perf,:Coherence,:Theta,:N,:Type_TC,:Alpha,:Percentage,:Stim])
	 for d in diclist

	 N= d["N"]
	 simulations = create_simu_from_dict(d)
	 mw1,mw2 = generate_mean_weights(simulations)
	 mw1 = mw1./0.5.*2.0

	 if d["type_TC"]=="Uniform"
	 	layer_G = init_layer(N)
	 	layer_U = init_layer(N)
	 else
	 	layer_G = init_layer(N,simulations.cats)
	 	layer_U = init_layer(N)
	 end
	# percentagetemp = percentagelist20[findfirst(x -> x==d["Seuil_confidence"],[10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0])]
	 percentagetemp = percentagelist15[findfirst(x -> x==d["Seuil_confidence"],[5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0])]

	 for stim in 0.0:0.01:0.5

	 Δp,coh =  performance(stim,layer_G.centers,layer_G.widths,(mw1[end,:]))
	 perf = Δp[1]
	 cohf = coh[1]
	 push!(df,[perf,cohf,d["Seuil_confidence"],d["N"],d["type_TC"],d["α"],percentagetemp,stim])
	 push!(df,[perf,cohf,d["Seuil_confidence"],d["N"],d["type_TC"],d["α"],percentagetemp,1.0-stim])

	 end
	 end

	return df
end

#<<<<<<<<<<<<




simulationstemp = Dict(
    "N" => [i for i=10:40], # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => [10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0], # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
    "α" => [0.10,0.15,0.20,0.22,0.25,0.27,0.30,0.34,0.40] # largeur des catégories gaussiennes
)


diclist = dict_list(simulationstemp)


dftest=fixed_stim([0.3,0.35,0.4,0.45,0.48,0.5],diclist)


df_fulldata = full_data(diclist)

CSV.write(datadir("notebooks","01-temptest.csv"),dftest)

CSV.write(datadir("notebooks","01-fulldata-20Hz.csv"),df_fulldata)


simulationstemp = Dict(
	"N" => [i for i=10:40], # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => [5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0], # seuil sur la confiancep our l'appr.
    "seuil_decision" => 15.0, # seuil de la décision
    "type_TC" => ["Uniform","Optimized"],
	"α" => [0.10,0.15,0.20,0.22,0.25,0.27,0.30,0.34,0.40] # largeur des catégories gaussiennes
)


diclist = dict_list(simulationstemp)
df_fulldata = full_data(diclist)
CSV.write(datadir("notebooks","01-fulldata-15Hz.csv"),df_fulldata)

dftest=fixed_stim([0.3,0.35,0.4,0.45,0.48,0.5],diclist)


CSV.write(datadir("notebooks","01-temptest2.csv"),dftest)


#%%% Optimal weights

#<<<<<<<<

#%%% true performance (pmap version)


#<<<<<

using Colors

diclist = dict_list(simulationstemp)
print("test")
templist = []
pyplot()
fig = plot()
clist = colormap("Blues",20)
for (i,d) in enumerate(diclist)
    simulations = create_simu_from_dict(d)
    mw1,mw2 = generate_mean_weights(simulations)
    seuil = d["Seuil_confidence"]
    #temp =lines!(1:1:simulations.N,mw1[end,:],color = colorlist[i],linewidth=5)[end]
    #push!(templist,temp)
    plot!(fig,1:1:simulations.N,mw1[end,:],color=clist[5+i],label="$i",linewidth=4*upscale)
end
fig
a = plotsdir("01-Analysis-script01")
savefig(fig,"$a\\test3.svg")

scene
lg = AbstractPlotting.legend(templist,["Seuil $i" for i in 1:10])
st = Stepper(vbox(scene,lg),"$a\\test.svg")
l = lg[end]
l[:textsize]= 12
step!(st)
st
test=vbox(scene,lg)
center!(test)
using CairoMakie, FileIO
AbstractPlotting.current_backend[] = CairoMakie.CairoBackend("path/to/svg")

open("$a\\test2.svg", "w") do io
     show(io, MIME"image/svg+xml"(), test)
end


Makie.save("$a\\test2.png",test)
a = plotsdir("01-Analysis-script01")




function performance(domaine,centers,widths,gain)
coherencelist = zeros(length(domaine))
perflist = zeros(length(coherencelist))
gain2 = copy(gain)
#gain2[centers.>0.5] = -gain2[centers.>0.5]
for i=1:length(domaine)
	dummy = 0
	for j=1:length(centers)
		dummy = dummy .+gain2[j].*fi(domaine[i],centers[j],widths[j])
end
coherencelist[i] = dummy[1]

end
    beta = 1.29
    alpha_rescaled = 0.05
    for i = 1:length(coherencelist)
        if sign(coherencelist[i])<0
            perflist[i] =0.5*exp(-exp(beta.*log.(-coherencelist[i]/alpha_rescaled)))
        else
            perflist[i] =(1.0 - 0.5*exp(-exp(beta.*log.(coherencelist[i]/alpha_rescaled))))
        end

    end
    return perflist
end

function threshold_to_percentage(θlist)

	network = MyAttractorNetwork(0.02,0.3255,0.2609,0.0497,0.0497,0.2609,100.0,0.000641,2.0,20.0)

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


percentagelist = threshold_to_percentage([10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0])


Nlist = 11:1:40
domaine = 0:0.01:1.0
perf_Gc_Uc =  zeros(length(Nlist),length(domaine))
perf_G_U =  zeros(length(percentagelist))
m=0.2
typeO = "Optimized"
typeU = "Uniform"
typeH="Hebbian"
alpha_list=[0.15,0.18,0.25,0.30]#0.35]#,0.30]#,0.20]#,0.25,0.3,0.35]



simulationstemp = Dict(
    "N" => 30, # nombre de neurones dans le codage
    "T_simu" => 10000,
    "ΔTsave" => 100, # pas de temps de sauvegarde
    "Seuil_confidence" => [10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0], # seuil sur la confiancep our l'appr.
    "seuil_decision" => 20.0, # seuil de la décision
    "type_TC" => "Uniform",
    "α" => 0.35 # largeur des catégories gaussiennes
)
diclist=dict_list(simulationstemp)
fig=plot()
for (y,N) in enumerate([15,20,25,35])
	for (j,d) in enumerate(diclist)

        d["N"]=N
        simulations = create_simu_from_dict(d)
        mw1,mw2 = generate_mean_weights(simulations)
        dtemp = copy(d)
        dtemp["Seuil_confidence"]=20.0
        simulations = create_simu_from_dict(dtemp)
        mw1SC,mw2SC = generate_mean_weights(simulations)


		layer_G = init_layer(N,simulations.cats)
		layer_U = init_layer(N)




	FcodelistGC = Float64[]
	FcodelistG = Float64[]
	FcodelistUC = Float64[]
	FcodelistU = Float64[]
	x=0.48
				push!(FcodelistGC,performance(x,layer_U.centers,layer_U.widths,(mw1[end,:]))[1])
				push!(FcodelistG,performance(x,layer_U.centers,layer_U.widths,(mw1SC[end,:]))[1])

	perf_G_U[j] = FcodelistGC[1] .- FcodelistG[1]
	end
	scatter!(fig,percentagelist,perf_G_U,label="N=$N",marker=:dot,markersize=4*upscale)
	print("\n")
end
fig
