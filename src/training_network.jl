
include(srcdir("structures.jl"))
include(srcdir("dynamical_system.jl"))
include(srcdir("information-theory.jl"))
include(srcdir("optimal-tuning-curves.jl"))


function init_layer(N::Int)
    centers = LinRange(0.0,1.0,N)
    widths = 1.0./N.*ones(N)
    wt1 = 1.0.-2.0.*rand(N)
    wt2 = 1.0.-2.0.*rand(N)
    codinglayer = MyCodingLayer(centers,widths,wt1,wt2,ones(N))
    return codinglayer
end

function init_layer(N::Int,cats::MyCatGauss)
    codinglayer,c=constant_gain(N,cats)
    codinglayer.weightsto1 = 1.0.-2.0.*rand(N)
    codinglayer.weightsto2 = 1.0.-2.0.*rand(N)
    return codinglayer
end

function init_layer(N::Int,cats::MyCatWeibull)
    codinglayer,c=constant_gain(N,cats)
    codinglayer.weightsto1 = 1.0.-2.0.*rand(N)
    codinglayer.weightsto2 = 1.0.-2.0.*rand(N)
    return codinglayer
end

function init_layer_Hebbian(N::Int,cats::MyCatWeibull)
    codinglayer,c=gain_Hebbian(N,cats)
    codinglayer.weightsto1 = 1.0.-2.0.*rand(N)
    codinglayer.weightsto2 = 1.0.-2.0.*rand(N)
    return codinglayer
end

function init_layer_Hebbian(N::Int,cats::MyCatGauss)
    codinglayer,c=gain_Hebbian(N,cats)
    codinglayer.weightsto1 = 1.0.-2.0.*rand(N)
    codinglayer.weightsto2 = 1.0.-2.0.*rand(N)
    return codinglayer
end


function responsetostim(stim::Float64,codinglayer::MyCodingLayer)
    response = zeros(length(codinglayer.centers))
    for i=1:length(codinglayer.centers)
        response[i] = exp(-(stim-codinglayer.centers[i])*(stim-codinglayer.centers[i])/
        (2*codinglayer.widths[i]*codinglayer.widths[i]))
    end
        return response
end

function inputCodingtoNetwork(codinglayer::MyCodingLayer,response::Array{Float64})
    return sum(codinglayer.weightsto1.*response),sum(codinglayer.weightsto2.*response)
end


function iscategorie(stim::Float64,cats::MyCatGauss,indice::Int)
    # indice: numéro de la catégorie dont on souhaite la proba
    N = length(cats.centers)
     pxmutotal = zeros(N);
     Px = 0.0;
     for k=1:N
      pxmutotal[k] = (1.0./(cats.width[k].*sqrt(2.0*pi))).*exp(-((stim-cats.centers[k]).^2)./(2.0.*(cats.width[k]^2))) # (nbc,N) vector of P( x|mu) for all	categories
     end
      Px = sum(cats.proba.*pxmutotal[:])#  probability of x, whatever the categorie
     z = pxmutotal[indice].*cats.proba[indice]./Px

     bool = 0
     if rand()<z
         bool = 1
     else
         bool = 2
     end
     return bool
end


function iscategorie(stim::Float64,cats::MyCatWeibull,indice::Int)
    # indice: numéro de la catégorie dont on souhaite la proba

     z = pmux(stim,cats,indice)

     bool = 0
     if rand()<z
         bool = 1
     else
         bool = 2
     end
     return bool
end

function normalization!(codinglayer::MyCodingLayer,rate::Float64)

    # produit sclaaire
    sum_w1 = BLAS.dot(length(codinglayer.weightsto1),codinglayer.weightsto1,1,codinglayer.weightsto1,1)
    sum_w2 = BLAS.dot(length(codinglayer.weightsto2),codinglayer.weightsto2,1,codinglayer.weightsto2,1)

    sum_w1 = sqrt(sum_w1)
    sum_w2 = sqrt(sum_w2)

    #normalization
    BLAS.axpby!(0.0,codinglayer.weightsto1,rate/sum_w1,codinglayer.weightsto1)
    BLAS.axpby!(0.0,codinglayer.weightsto2,rate/sum_w2,codinglayer.weightsto2)

end




function learning(simulations::ListParameters_θconfidenceSimulation)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,70)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

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
    for i=1:simulations.T_simu

        indice = rand(1:1:70)
        stim = stimrange[indice]
        response = responsetostim(stim,codinglayer)
        i1,i2=inputCodingtoNetwork(codinglayer,response)
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(i1,i2,30.0,0.00052)
        r1,r2 = simu(unit1,unit2,inputToUnits,network,timeparameters,tf)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
         @views Δr = mean( r1 .- r2)
        aΔr = abs(Δr)
            if aΔr>simulations.Seuil_confidence
            else
                if Δr>0.0 && c==1
                    reward = 1.0
                  codinglayer.weightsto1 +=simulations.learning_rate.*(reward).*response
                  codinglayer.weightsto2 -=simulations.learning_rate.*(reward).*response

                elseif Δr < 0.0 && c==2
                    reward = 1.0
                  codinglayer.weightsto1 -=simulations.learning_rate.*(reward).*response
                  codinglayer.weightsto2 +=simulations.learning_rate.*(reward).*response
                elseif Δr < 0.0
                    reward = -1.0
                  codinglayer.weightsto1 -=simulations.learning_rate.*(reward).*response
                  codinglayer.weightsto2 +=simulations.learning_rate.*(reward).*response
                else
                    reward = -1.0
                  codinglayer.weightsto1 +=simulations.learning_rate.*(reward).*response
                  codinglayer.weightsto2 -=simulations.learning_rate.*(reward).*response
                end
                normalization!(codinglayer,simulations.cst_normalisation)
            end

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end


function learning_poisson_input(simulations::ListParameters_θconfidenceSimulation,box=1.0)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,100)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

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
    for i=1:simulations.T_simu

        indice = rand(1:1:100)
        stim = stimrange[indice]
        #response = responsetostim(stim,codinglayer)
        #i1,i2=inputCodingtoNetwork(codinglayer,response)
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(0.0,0.0,30.0,0.00052) 
        r1,r2,response = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
         Δr = mean( r1 .- r2)
        aΔr = abs(Δr)
            if aΔr>simulations.Seuil_confidence
            else
                if Δr>0.0 && c==1
                    reward = 1.0
                    BLAS.axpy!(simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)


                elseif Δr < 0.0 && c==2
                    reward = 1.0
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

                elseif Δr < 0.0
                    reward = -1.0
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

                else
                    reward = -1.0
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

                end
                normalization!(codinglayer,simulations.cst_normalisation)
            end

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            normalization!(codinglayer,simulations.cst_normalisation)

            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end


function learning(time::Int,confidenceth::Float64,λ::Float64,timesave::Int,N::Int,cats::MyCatWeibull,rate::Float64,type::String)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,70)
    timelistto1 = zeros(Int(time/timesave),N)
    timelistto2 = zeros(Int(time/timesave),N)

    network = MyAttractorNetwork(0.02,0.3255,0.2609,0.0497,0.0497,0.2609,100.0,0.000641,2.0,15.0)
    timeparameters = MyEulerParameters(0.5,4.0,4.0)

    tf = MyTransferFunction(270.0,108.0,0.1540)

if type == "Uniform"
    codinglayer =init_layer(N)
else
    codinglayer =init_layer(N,cats)

end
    for i=1:time

        indice = rand(1:1:70)
        stim = stimrange[indice]
        response = responsetostim(stim,codinglayer)
        i1,i2=inputCodingtoNetwork(codinglayer,response)
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(i1,i2,30.0,0.00052)
        r1,r2 = simu(unit1,unit2,inputToUnits,network,timeparameters,tf)
        reward = 0.0
        c = iscategorie(stim,cats,1)
            if abs.(r1[end] - r2[end])>confidenceth
            else
                if r1[end]>r2[end] && c==1
                    reward = 1.0
                    codinglayer.weightsto1 +=λ.*(reward).*response
                    codinglayer.weightsto2 -=λ.*(reward).*response

                elseif r1[end]<r2[end] && c==2
                    reward = 1.0
                    codinglayer.weightsto1 -=λ.*(reward).*response
                    codinglayer.weightsto2 +=λ.*(reward).*response
                elseif r1[end]<r2[end]
                    reward = -1.0
                    codinglayer.weightsto1 -=λ.*(reward).*response
                    codinglayer.weightsto2 +=λ.*(reward).*response
                else
                    reward = -1.0
                    codinglayer.weightsto1 +=λ.*(reward).*response
                    codinglayer.weightsto2 -=λ.*(reward).*response
                end
                normalization!(codinglayer,rate)
            end

        r1 = nothing
        r2 = nothing
        if i%timesave ==0
            timelistto1[Int(i/timesave),:] = codinglayer.weightsto1
            timelistto2[Int(i/timesave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end
### Input stochastique (Poisson train)
using Distributions

function gen_input!(response::Array{Float64},codinglayer::MyCodingLayer,input::MyInput)
    # la fonction retourne un vecteur correspondant à un nombre de spike
    #l'unité est le ms
    # response : mean firing rate of the ocding layer
        N=length(codinglayer.centers)
    input.i1 = sum(response.*codinglayer.weightsto1)
    input.i2 = sum(response.*codinglayer.weightsto2)

end

function rate_feedback(unit1::MyUnit,unit2::MyUnit,stim::Float64,codinglayer::MyCodingLayer)
    # la fonction retourne un vecteur correspondant à un nombre de spike
    #l'unité est le ms
    # response : mean firing rate of the ocding layer
        N=length(codinglayer.centers)
        res = responsetostim(stim,codinglayer)
        meanrate = zeros(N)
        for i=1:N
            meanrate[i] = max(0.0,res[i] + (codinglayer.weightsto1[i]*unit1.rate + codinglayer.weightsto2[i]*unit2.rate)/10)
        end
        return meanrate
end


function poisson_train(meanrate::Array{Float64},fmin=0.001,fmax=1.0)
r = zeros(length(meanrate))
    for i=1:length(meanrate)
        r[i] = rand(Poisson(fmin.+(fmax-fmin).*meanrate[i]),1)[1]
    end
    return r
end

function stochastic_attractor_network(stim::Float64,unit1::MyUnit,unit2::MyUnit,input::MyInput,network::MyAttractorNetwork,timeparameters::MyEulerParameters,transferparameters::MyTransferFunction,codinglayer::MyCodingLayer,feedback=false)
      decision = false
      rate1 = Float64[]
      rate2 = Float64[]
      savingrate = []
      t = 0
      if feedback
          while decision == false && t < 4000
            res = rate_feedback(unit1,unit2,stim,codinglayer)
            response = poisson_train(res)
            gen_input!(response,codinglayer,input)
            update_isyn!(unit1,unit2,network,input)
            update_noise!(unit1,timeparameters.dt,network.tampa,network.i0,network.noiseamp)
            update_noise!(unit2,timeparameters.dt,network.tampa,network.i0,network.noiseamp)
            update_s!(unit1,unit2,network,timeparameters,transferparameters)
            decision = check_threshold(unit1,unit2,network)
            push!(rate1,unit1.rate)
            push!(rate2,unit2.rate)
            t +=timeparameters.dt
            push!(savingrate,res)
          end
          return rate1,rate2,savingrate

      else
      while decision == false && t < 4000
        res = responsetostim(stim,codinglayer)
        response = poisson_train(res)
        gen_input!(response,codinglayer,input)
        update_isyn!(unit1,unit2,network,input)
        update_noise!(unit1,timeparameters.dt,network.tampa,network.i0,network.noiseamp)
        update_noise!(unit2,timeparameters.dt,network.tampa,network.i0,network.noiseamp)
        update_s!(unit1,unit2,network,timeparameters,transferparameters)
        decision = check_threshold(unit1,unit2,network)
        push!(rate1,unit1.rate)
        push!(rate2,unit2.rate)
        t +=timeparameters.dt
      end
        return rate1,rate2
  end


end




function learning_poisson_input_reward_modulated(simulations::ListParameters_θconfidenceSimulation,box=1.0)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation
    # Modulation of reward following Engel et al 2015
    # still consider Normalization to compare with our model

    stimrange = LinRange(0.0,1.0,100)
    reward_list = zeros(100)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

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
    for i=1:simulations.T_simu

        indice = rand(1:1:100)
        stim = stimrange[indice]
     
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(0.0,0.0,30.0,0.00052) 
        r1,r2,response = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
         Δr = mean( r1 .- r2)
        aΔr = abs(Δr)

                if Δr>0.0 && c==1
                    reward = 1.0 - reward_list[indice]
                    BLAS.axpy!(simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)
                    reward_list[indice] = reward_list[indice] + (reward - reward_list[indice] )/5.0


                elseif Δr < 0.0 && c==2
                    reward = 1.0 - reward_list[indice]
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)
                    reward_list[indice] = reward_list[indice] + (reward - reward_list[indice] )/5.0

          
                elseif Δr < 0.0
                    reward = -1.0 - reward_list[indice]
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)
                    reward_list[indice] = reward_list[indice] + (reward - reward_list[indice] )/5.0

             
                else
                    reward = -1.0 - reward_list[indice]
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)
                    reward_list[indice] = reward_list[indice] + (reward - reward_list[indice] )/5.0

                 
                end
                normalization!(codinglayer,simulations.cst_normalisation)

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            normalization!(codinglayer,simulations.cst_normalisation)

            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end




function learning_poisson_input_unsupervised(simulations::ListParameters_θconfidenceSimulation,box=1.0)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,100)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

    network = MyAttractorNetwork(threshold=simulations.seuil_decision)
    timeparameters = MyEulerParameters(0.5,4.0,4.0)
    tf = MyTransferFunction(270.0,108.0,0.1540)

# initialisation a changer

    if simulations.type_TC == "Uniform"
        codinglayer =init_layer(simulations.N)
    elseif simulations.type_TC == "Hebbian"
        codinglayer =init_layer_Hebbian(simulations.N,simulations.cats)

    else
        codinglayer =init_layer(simulations.N,simulations.cats)

    end

    # symmetrie des poids
    codinglayer.weightsto1 = (codinglayer.weightsto1 .+ reverse(codinglayer.weightsto1))./2.0
    codinglayer.weightsto2 = -1.0.*codinglayer.weightsto1


    for i=1:simulations.T_simu

        indice = rand(1:1:100)
        stim = stimrange[indice]
        #response = responsetostim(stim,codinglayer)
        #i1,i2=inputCodingtoNetwork(codinglayer,response)
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(0.0,0.0,30.0,0.00052) 
        r1,r2,response = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
         Δr = mean( r1 .- r2)
        aΔr = abs(Δr)
 

            if aΔr>simulations.Seuil_confidence
            else
               
                    reward = 1.0
                    BLAS.axpy!(simulations.learning_rate.*mean(r1)./20.0,response,codinglayer.weightsto1)
                    BLAS.axpy!(simulations.learning_rate.*mean(r2)./20.0,response,codinglayer.weightsto2)
                normalization!(codinglayer,simulations.cst_normalisation)
            end

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            normalization!(codinglayer,simulations.cst_normalisation)

            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end




function learning_poisson_input_unsupervised_hard_init(simulations::ListParameters_θconfidenceSimulation,box=1.0)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,100)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

    network = MyAttractorNetwork(threshold=simulations.seuil_decision)
    timeparameters = MyEulerParameters(0.5,4.0,4.0)
    tf = MyTransferFunction(270.0,108.0,0.1540)

# initialisation a changer

    if simulations.type_TC == "Uniform"
        codinglayer =init_layer(simulations.N)
    elseif simulations.type_TC == "Hebbian"
        codinglayer =init_layer_Hebbian(simulations.N,simulations.cats)

    else
        codinglayer =init_layer(simulations.N,simulations.cats)

    end

    # symmetrie des poids
    codinglayer.weightsto1 = (abs.(codinglayer.weightsto1) .+ abs.(reverse(codinglayer.weightsto1)))./2.0
  
    codinglayer.weightsto1[Int(round((simulations.N)./2.0+0.45)):end] = -1.0.*codinglayer.weightsto1[1:round(Int, (simulations.N)./2.0)]
    codinglayer.weightsto2 = -1.0.*codinglayer.weightsto1


    for i=1:simulations.T_simu

        indice = rand(1:1:100)
        stim = stimrange[indice]
        #response = responsetostim(stim,codinglayer)
        #i1,i2=inputCodingtoNetwork(codinglayer,response)
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(0.0,0.0,30.0,0.00052) 
        r1,r2,response = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
         Δr = mean( r1 .- r2)
        aΔr = abs(Δr)
   

            if aΔr>simulations.Seuil_confidence
            else
               
                    reward = 1.0
                    BLAS.axpy!(simulations.learning_rate.*mean(r1)./20.0,response,codinglayer.weightsto1)
                    BLAS.axpy!(simulations.learning_rate.*mean(r2)./20.0,response,codinglayer.weightsto2)
                normalization!(codinglayer,simulations.cst_normalisation)
            end

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            normalization!(codinglayer,simulations.cst_normalisation)

            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end



function learning_poisson_input_corrupted(simulations::ListParameters_θconfidenceSimulation,box=1.0)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,100)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

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
    for i=1:simulations.T_simu

        indice = rand(1:1:100)
        stim = stimrange[indice]
     
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(0.0,0.0,30.0,0.00052) 
        r1,r2,response = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
        test = rand()
        if  test < simulations.corrupted
            c = rand(1:2)
        end
         Δr = mean( r1 .- r2)
        aΔr = abs(Δr)
            if aΔr>simulations.Seuil_confidence
            else
                if Δr>0.0 && c==1
                    reward = 1.0
                    BLAS.axpy!(simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)



                elseif Δr < 0.0 && c==2
                    reward = 1.0
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)


                elseif Δr < 0.0
                    reward = -1.0
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

              
                else
                    reward = -1.0
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

                
                end
                normalization!(codinglayer,simulations.cst_normalisation)
            end

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            normalization!(codinglayer,simulations.cst_normalisation)

            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end



function learning_poisson_input_corrupted_constant(simulations::ListParameters_θconfidenceSimulation,box=1.0)
    # confidence denotes the threshold of confidence level in order
    # to still perform the learning step
    # for now only a learning with normalisation

    stimrange = LinRange(0.0,1.0,100)
    timelistto1 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)
    timelistto2 = zeros(Int(simulations.T_simu/simulations.ΔTsave),simulations.N)

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
    for i=1:simulations.T_simu

        indice = rand(1:1:100)
        stim = stimrange[indice]
        #response = responsetostim(stim,codinglayer)
        #i1,i2=inputCodingtoNetwork(codinglayer,response)
        unit1 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        unit2 = MyUnit(2.0,0.1,0.0,network.noiseamp*randn(), [0.0,0.0,0.0,0.0])
        inputToUnits = MyInput(0.0,0.0,30.0,0.00052)
        r1,r2,response = simu_stochastic(unit1,unit2,inputToUnits,network,timeparameters,tf,codinglayer,stim,box)
        reward = 0.0
        c = iscategorie(stim,simulations.cats,1)
        test = rand()
        if stim <0.5
        z = pmux(stim,simulations.cats,2)
        else
            z = pmux(stim,simulations.cats,1)

        end
        if  test < 2.0*( simulations.corrupted - z)
            c = rand(1:2)
        end
         Δr = mean( r1 .- r2)
        aΔr = abs(Δr)
            if aΔr>simulations.Seuil_confidence
            else
                if Δr>0.0 && c==1
                    reward = 1.0
                    BLAS.axpy!(simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)


                elseif Δr < 0.0 && c==2
                    reward = 1.0
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

                elseif Δr < 0.0
                    reward = -1.0
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

              
                else
                    reward = -1.0
                    BLAS.axpy!(1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto1)
                    BLAS.axpy!(-1.0*simulations.learning_rate.*reward,response,codinglayer.weightsto2)

                  
                end
                normalization!(codinglayer,simulations.cst_normalisation)
            end

        r1 = nothing
        r2 = nothing
        if i%simulations.ΔTsave ==0
            normalization!(codinglayer,simulations.cst_normalisation)

            timelistto1[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto1
            timelistto2[Int(i/simulations.ΔTsave),:] = codinglayer.weightsto2

        end


    end
    return timelistto1,timelistto2
end
