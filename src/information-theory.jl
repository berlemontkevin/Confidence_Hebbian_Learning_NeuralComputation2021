using DrWatson
#quickactivate("@__DIR__","4-project_coding_categories")
# now can load Packages

# In this file we can find the necessary function to compute some information theory related quantities
# as well as the optimization of tuning curves following Bonasse Gahot and Nadal
include(srcdir("structures.jl"))
include(srcdir("dynamical_system.jl"))


using Random
Random.seed!(123)


#%%All functions for Gaussian categories
function pmux(stim::Float64,cats::MyCatGauss,indice::Int)
    N = length(cats.centers)
     pxmutotal = zeros(N);
     Px = 0.0;
     for k=1:N
      pxmutotal[k] = (1.0./(cats.width[k].*sqrt(2.0*pi))).*exp(-((stim-cats.centers[k]).^2)./(2.0.*(cats.width[k]^2))) # (nbc,N) vector of P( x|mu) for all	categories
     end
      Px = sum(cats.proba.*pxmutotal[:])#  probability of x, whatever the categorie
     z = pxmutotal[indice].*cats.proba[indice]./Px
    return z
end

function pxmu(x::Float64,c::Float64,w::Float64)
    """# we suppose the categories are gaussian"""
    dummy = (1.0./(w.*sqrt(2.0*pi))).*exp(-((x-c).^2)./(2.0.*(w^2)))
    return dummy
end

function px(x::Float64,cats::MyCatGauss)
    dummy = 0.0
    for i=1:length(cats.centers)
        dummy += cats.proba*pxmu(x,cats.centers[i],cats.width[i])
    end
    return dummy
end
function cumulative_integral(l::Array{Float64})
    # compute the cumulative integral
    sl = zeros(length(l))
    for i=1:length(l)
        sl[i] = sum(l[1:i])
    end
    return sl
end

function Fcat(x::Float64,cats::MyCatGauss)
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
        pmukx = pmux(x,cats,k)

        if (pmukx >th)
            difpmukx = difpmux(x,k,cats)
            fcate = fcate .+ 1.0./pmukx.*difpmukx.*difpmukx  # Direct calcul of the 2-dif logP

        end
    end

    return fcate[1]
end
function difpmux(x::Float64,n::Int64,cats::MyCatGauss)
    """ Compute the differential of Pmux : use for Fcat
    # We compute first Pxmu (easier) and then we use Bayes rule to finish the calcul
    # n : categorie of interest
    # xmean : mean position of the several gaussian categories
    # width : list of the width of the several categories
    # Pmu : different categories probability
    """

    nbC = length(cats.centers) # number of categories
    pxmuk=zeros(nbC)
    difpxmuk=zeros(nbC)
    for k=1:nbC
        pxmuk[k] = pxmu(x,cats.centers[k],cats.width[k])
        difpxmuk[k] = - (x - cats.centers[k])./(cats.width[k]^2).*pxmuk[k] # Direct computation of the analytic expression for the differential

    end

    difpx = 0.0 ; px = 0.0;
    for k=1:nbC
        difpx = difpx .+ cats.proba[k].*difpxmuk[k]  # differential of the probability of x
        px = px .+ cats.proba[k].*pxmuk[k]
    end

    p = cats.proba[n].*difpxmuk[n]./px - cats.proba[n].*pxmuk[n]./(px.^2).*difpx # direct computation of the quantity
    return p
end
#<<<<<<<<<<<<<<<<<


#%%All functions for Weibull categories
function pmux(x::Float64,cat::MyCatWeibull,n::Int64)
    """# Compute P(mu |x) for the Weibull distribution
    # n : indice of the categorie we want to comput
    # x : 1,N vectors of the points
    # Pmu : proba of each categorie
    # xmean : list of the centers of the gaussian distribution function
    # width : list of the width of all gaussian distribution function
    """
    if x < cat.centers[n] && n == 1
        Px = 1.0 - 0.5*exp(-(abs((x - cat.centers[n]))/cat.alpha)^(cat.beta))
        return Px

    elseif x > cat.centers[n] && n == 2
        Px = 1.0 - 0.5*exp(-(abs((x - cat.centers[n]))/cat.alpha)^(cat.beta))
        return Px


    elseif x > cat.centers[n] && n == 1
        Px =  0.5*exp(-(abs((x - cat.centers[n]))/cat.alpha)^(cat.beta))
        return Px


    else x < cat.centers[n] && n == 2
        Px =  0.5*exp(-(abs((x - cat.centers[n]))/cat.alpha)^(cat.beta))
        return Px

    end

end

function Fcat(x::Float64,cats::MyCatWeibull)
    """# Compute the Fcat of our Weibull categories
    # We compute Fcat using P'^2/P
    # x : input where Fcat is computing
    """
    nbC = length(cats.centers) # number of categories
    fcate = 0 # initialization
    th = 1e-6  #threshold to know when stop the integration

    for k=1:nbC # sum on all the categories
        pmukx = pmux(x,cats,k)
        if (pmukx >th)
            difpmukx = difpmux(x,k,cats)
            fcate = fcate + 1.0/pmukx*difpmukx*difpmukx  # Direct calcul of the 2-dif logP
        end
    end
    return fcate
end

function difpmux(x::Float64,n::Int64,cats::MyCatWeibull)
    """ Compute the differential of Pmux : use for Fcat for Weibull
    # We compute first Pxmu (easier) and then we use Bayes rule to finish the calcul
    # n : categorie of interest
    # We compute in absolute value the derivative
    """
        Px = 0.5*exp(-(abs((x - cats.centers[n]))/cats.alpha)^(cats.beta))*(abs((x - cats.centers[n]))/cats.alpha)^(cats.beta-1)
        return Px
end


#<<<<<<<<<<<<<<<<<


#%% Compute Fcode for gaussian tuning curves
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


#<<<<<<<<<<
