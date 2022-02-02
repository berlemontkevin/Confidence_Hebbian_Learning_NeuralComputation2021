using DrWatson
#quickactivate("@__DIR__","4-project_coding_categories")
# now can load Packages

# In this file we can find the necessary function to compute some information theory related quantities
# as well as the optimization of tuning curves following Bonasse Gahot and Nadal
# there is noise in the prob distribution
# I will use automatic differentiation
include(srcdir("structures.jl"))
include(srcdir("dynamical_system.jl"))
using ForwardDiff


function P_mu_x_corrupted_noise(x::Float64,cats::MyCatGauss)
    N = length(cats.centers)

    return 1.0/N
end


function pmux_gaussian(stim::Float64,cats::MyCatGauss,indice::Int)
    N = length(cats.centers)
     pxmutotal = zeros(N);
     Px = 0.0;
    pxmutotal = (1.0./(cats.width.*sqrt(2.0*pi))).*exp.(-((stim.-cats.centers).^2)./(2.0.*(cats.width.^2))) # (nbc,N) vector of P( x|mu) for all	categories
   
      Px = sum(cats.proba.*pxmutotal[:])#  probability of x, whatever the categorie
     z = pxmutotal[indice].*cats.proba[indice]
    return z
end


function pmux(stim::Float64,cats::MyCatGauss,indice::Int,p_corrupted::Float64)
    p_c = P_mu_x_corrupted_noise(stim,cats)

    N = length(cats.centers)
    pxmutotal = zeros(N);
    Px = 0.0;
   pxmutotal = (1.0./(cats.width.*sqrt(2.0*pi))).*exp.(-((stim.-cats.centers).^2)./(2.0.*(cats.width.^2))) # (nbc,N) vector of P( x|mu) for all	categories
  
     Px = sum(cats.proba.*pxmutotal[:])#  probability of x, whatever the categorie
    z = pxmutotal[indice].*cats.proba[indice]
    y = p_corrupted*p_c + (1.0-p_corrupted).*z
    ysum = sum(p_corrupted*p_c .+ (1.0-p_corrupted).*pxmutotal[:].*cats.proba[:])
        return y./ysum
end

function cumulative_integral(l::Array{Float64})
    # compute the cumulative integral
    sl = zeros(length(l))
    for i=1:length(l)
        sl[i] = sum(l[1:i])
    end
    return sl
end


function Fcat(x::Float64,cats::MyCatGauss,p_corrupted::Float64)
    """# Compute the Fcat of our gaussian categories
    # We compute Fcat using P'^2/P
    # x : input where Fcat is computing
    # xmean : mean position of the several gaussian categories
    # width : list of the width of the several categories
    # Pmu : different categories probability
    """
    nbC = length(cats.centers) # number of categories
    fcate = 0 # initialization
    th = 1e-6  #threshold to know when stop the integration

    for k=1:nbC # sum on all the categories
        pmukx = pmux(x,cats,k,p_corrupted)

        if (pmukx >th)
            a = gradient(x -> pmux(x, cats,k,p_corrupted), x)
            difpmukx = a[1]
            fcate = fcate .+ 1.0./pmukx.*difpmukx.*difpmukx  # Direct calcul of the 2-dif logP

        end
    end

    return fcate[1]
end



function constant_gain(N::Int,cat::MyCatGauss,p_corrupted::Float64)
	# return a coding layer structure for a Weibull distribution of cat
	x = range(0.0,length=400,1.0)
	domaine = Float64[]
	for a in x
	push!(domaine,a)
	end
	Fcatlist = Float64[]
	for xi in x
				push!(Fcatlist,Fcat(xi,cat,p_corrupted))
			end

	p = 1.0./400.0.*ones(400) # uniform prob. of stim

	dstemp = (p.*Fcatlist).^(1/3) #
	ds = dstemp./(sum(dstemp.*1.0./400)).*(N-1)
	D = cumulative_integral(ds./400)
	width = 1.0./ds

	centers = zeros(N)
	wf = zeros(N)
	for i=1:N-1
		ind = minimum(findall(x->x>=i-0.005,D))
	centers[i+1],wf[i+1] = domaine[ind],width[ind]
	end
	wf[1]=wf[N] # adding symmetry
	centers[1] = 1.0-centers[N]

	mylayer = MyCodingLayer(centers,wf,zeros(N),zeros(N),1.0.*ones(N))


	return mylayer,Fcatlist
end


using Zygote
N = 20
α=0.25
mycatGauss = MyCatGauss(width=[α,α])

codinglayer = init_layer(N,mycatGauss)

function TC_value(stim,center,width)
    TC = zeros(length(stim))
    for i=1:length(stim)
        TC[i] = exp(-(stim[i]-center)*(stim[i]-center)/(2*width*width))
    end

    return TC
end
p_corrupted = 0.0
function generate_TC(N,α,p_corrupted)
    mycatGauss = MyCatGauss(width=[α,α])

    codinglayer,t = constant_gain(N,mycatGauss,p_corrupted)
    df = DataFrame(stim=Float64[],N=Int[],TC=Float64[],type=String[])
    stimlist = [0.001*i for i=0:1000]

    for j=1:length(codinglayer.centers)
             TC = TC_value(stimlist,codinglayer.centers[j],codinglayer.widths[j])
             for i=1:length(stimlist)
                 push!(df,[stimlist[i],j,TC[i],"Optimized"])
             end
    end
 
    return df
end

using StatsPlots


N = 20
α=0.2
p_corrupted = 0.2

df = generate_TC(N,α,p_corrupted)

plot()
@df df scatter(:stim,:TC,group=:type,linewidth=2)

function fi(x,position,width)

    # Compute the response curve of the neuron (tuning curve)
    # In this case, gaussian one
    # x: iput on the neuron
    # position : mean of the tuning curve of the neuron (gaussian)
    # width : width of the neurons tuning curve
   
    f = zeros(length(x))
    for k=1:length(x)
          f[k] = exp.(-((x[k].-position).^2)./(2.0.*width.^2))
    end
   
    return f
end
   
function fisher_dwidth(x,pos,width,indice)
    # Compute the gradient relativ to width
    # indice : numéro du layer
    # x : input stim
    # pos : mean of the gaussians tuning curve
    # width : width of the gaussian tuning curve
   
    F = Fcode(x,pos,width)
    fxi = fi(x,pos[indice],width[indice])  # compute the response of the ieme layer
    dfix = -(x-pos[indice])/(width[indice]^2) # Differentiation of the formula respect to pos
    dwix = (x - pos[indice])^2/(width[indice]^3) # differention term for w
    dfxidwxi = 2.0*(x-pos[indice])/(width[indice]^3)# second order differentiation
    DFxdwi  = dfix*fxi*(2.0*dfxidwxi+dfix*dwix) # compute the analytical derivative of F respect to wi
    ft = -DFxdwi/(F.^2)
   
     return ft
end
   
   
function fisher_dpos(x,pos,width,indice)
    # Compute the gradient relativ to position of the fisher information
    # indice : numéro du layer
    # x : input stim
    # pos : mean of the gaussians tuning curve
    # width : width of the gaussian tuning curve
   
    F = Fcode(x,pos,width)
   
    fxi = fi(x,pos[indice],width[indice]) # compute the response of the ieme layer
   
    dfix = -(x-pos[indice])/(width[indice]^2) # partial Differentiation of the formula
   
    DFxdxi  = dfix*fxi*(2.0/(width[indice]^2)-dfix^2) # compute the analytical derivative of F respect to xi ( in ordert o do gradient descent)
   
    ft = -DFxdxi/(F.^2)
   
   
    return ft
end
   
   
function Fcode(x,position,width,gain=ones(length(position)))
    # Compute the Fcode of the neural code
    # We use the formula fi'^2/fi
    # x : input on the neuron
    # position : mean of the gaussian tuning curve
    # width : width of the gaussian tuning curve
   
    fcod = 0 ; #initialization
    N = length(position);  # number of codig layer
   
    for j=1:N # sum on all the neurons
   
        fjx = gain[j].*fi(x,position[j],width[j]);
        dfjx = .-(x.-position[j])./(width[j]^2); # Differentiation of the formula
        fcod = fcod .+ (dfjx.^2).*fjx ; # we ahve performed an integration over r which are poissonian
   
    end
   
    return fcod[1]
end



function plot_Fcat(α::Float64)
    mycatGauss = MyCatGauss(width=[α,α])
    N=20
    fig =plot(color_palette=:viridis)
    for p_corrupted in [0.0, 0.1,0.2,0.5,0.7,0.9]
        
        codinglayer,t = constant_gain(N,mycatGauss,p_corrupted)
        #plot!(fig,t./maximum(t),line_z=p_corrupted,color=:viridis,
        #linewidth=2, label="p_c = $p_corrupted")
        plot!(fig,t,line_z=p_corrupted,color=:viridis,
        linewidth=2, label="p_c = $p_corrupted")
    end
    title!("Normalized Fcat alpha=$α")
    return fig
end

