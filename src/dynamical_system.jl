using Distributions,DataFrames
using LinearAlgebra
using LinearAlgebra.BLAS
include(srcdir("structures.jl"))



function  f(isyn::Float64,parameters::MyTransferFunction)
  dummy= (parameters.a*isyn-parameters.b)/(1.0-exp(-parameters.d*(parameters.a*isyn-parameters.b)))
  if dummy <0 # firing rates can't be negativ
    dummy = 0.0
  end
  return dummy
end


function update_noise!(unit::MyUnit,dt::Float64,Tampa::Float64,iresting::Float64,noiseamp::Float64)
   unit.inoise += (dt/Tampa)*(iresting-unit.inoise) + sqrt(dt/Tampa)*noiseamp*randn()
  return
end


function update_isyn!(unit1::MyUnit,unit2::MyUnit,network::MyAttractorNetwork,input::MyInput)

  tounit1 = input.jaext*input.μ0*(1.0+max(input.i1,-1.0))
  tounit2 = input.jaext*input.μ0*(1.0+max(input.i2,-1.0))
  unit1.isyn_coming = network.jn11*unit1.s - network.jn21*unit2.s + tounit1 + unit1.inoise
  unit2.isyn_coming = network.jn22*unit2.s - network.jn12*unit1.s + tounit2 + unit2.inoise
end

function update_s!(unit1::MyUnit,unit2::MyUnit,network::MyAttractorNetwork,timeparameters::MyEulerParameters,transferparameters::MyTransferFunction,input::MyInput)

  # update syn current
  dummy_curr = input.jaext*input.μ0
  tounit1 = dummy_curr + dummy_curr*max(input.i1,-1.0)
  tounit2 = dummy_curr + dummy_curr*max(input.i2,-1.0)
  unit1.isyn_coming = network.jn11*unit1.s - network.jn21*unit2.s + tounit1 + unit1.inoise
  unit2.isyn_coming = network.jn22*unit2.s - network.jn12*unit1.s + tounit2 + unit2.inoise

  # update noise
  dummy_time = timeparameters.dt/network.tampa
  unit1.inoise += (dummy_time)*(network.i0-unit1.inoise) + sqrt(dummy_time)*network.noiseamp*randn()
  unit2.inoise += (dummy_time)*(network.i0-unit2.inoise) + sqrt(dummy_time)*network.noiseamp*randn()

  # update S variable and rate
  unit1.rate = f(unit1.isyn_coming,transferparameters)
  unit2.rate = f(unit2.isyn_coming,transferparameters)

  unit1.s += timeparameters.dt*(-unit1.s/network.tnmda + (1.0-unit1.s)*network.γ*unit1.rate)
  unit2.s += timeparameters.dt*(-unit2.s/network.tnmda + (1.0-unit2.s)*network.γ*unit2.rate)


@inbounds @simd  for i=1:3
  unit1.meanrate[i]=unit1.meanrate[i+1]
  unit2.meanrate[i]=unit2.meanrate[i+1]
  end
  unit1.meanrate[4] = unit1.rate
  unit2.meanrate[4] = unit2.rate
end


function check_threshold(unit1::MyUnit,unit2::MyUnit,network::MyAttractorNetwork)
  if mean(unit1.meanrate) > network.threshold
    bool = true
  elseif mean(unit2.meanrate) > network.threshold
    bool = true
  else
    bool = false
  end
  return bool
end

function simu(unit1::MyUnit,unit2::MyUnit,input::MyInput,network::MyAttractorNetwork,timeparameters::MyEulerParameters,transferparameters::MyTransferFunction)
  decision = false
  rate1 = Float64[]
  rate2 = Float64[]
  t = 0
  while decision == false && t < 4000
 
    update_s!(unit1,unit2,network,timeparameters,transferparameters,input)
    decision = check_threshold(unit1,unit2,network)
    push!(rate1,unit1.rate)
    push!(rate2,unit2.rate)
    t +=timeparameters.dt
  end
  return rate1,rate2
end

function update_noise_variance!(unit::MyUnit,dt::Float64,Tampa::Float64,iresting::Float64,noiseamp::Float64,input::MyInput,c::Float64)
  unit.inoise += (dt/Tampa)*(iresting-unit.inoise) + sqrt(dt/Tampa)*sqrt(noiseamp*noiseamp)*randn()
end

function update_isyn2!(unit1::MyUnit,unit2::MyUnit,network::MyAttractorNetwork,input::MyInput)

  tounit1 = input.jaext*input.μ0*(max(input.i1,-0.0))
  tounit2 = input.jaext*input.μ0*(max(input.i2,-0.0))
  unit1.isyn_coming = network.jn11*unit1.s - network.jn21*unit2.s + tounit1 + unit1.inoise
  unit2.isyn_coming = network.jn22*unit2.s - network.jn12*unit1.s + tounit2 + unit2.inoise
end

function update_spikes!(rates::Array{Float64},codinglayer::MyCodingLayer,input::MyInput,poissons_train::Array{Float64,1},box=1.0)
  # poissons_train : represent the trian of poisson spikes for each neuron of coding layer
  # it is updated at each time step ()
  # on divise par μ0 par soucis de normalisation des paramètres biologiques
  # defualt is 1s for box size of spike train
@inbounds @simd  for i=1:length(codinglayer.centers)
    poissons_train[i] = 1.0/(input.μ0*box)*rand(Poisson(rates[i]*box))
  end
end

function spikes_to_input!(input::MyInput,codinglayer::MyCodingLayer,response::Array{Float64})
    input.i1 = BLAS.dot(length(codinglayer.weightsto1),codinglayer.weightsto1,1,response,1)
    input.i2 = BLAS.dot(length(codinglayer.weightsto1),codinglayer.weightsto2,1,response,1)
end


function simu_stochastic(unit1::MyUnit,unit2::MyUnit,input::MyInput,network::MyAttractorNetwork,timeparameters::MyEulerParameters,transferparameters::MyTransferFunction,codinglayer::MyCodingLayer,stim::Float64,box=1.0)
  # we have change of variance depending on stimulus
  decision = false
  rate1 = Float64[]
  rate2 = Float64[]
  t = 0
  response = zeros(length(codinglayer.centers))
  poisson_train = zeros(length(codinglayer.centers))

  @inbounds @simd  for i=1:length(codinglayer.centers)
      response[i] = input.μ0*exp(-(stim-codinglayer.centers[i])*(stim-codinglayer.centers[i])/
      (2*codinglayer.widths[i]*codinglayer.widths[i]))
    end



  while decision == false && t < 4000
    if t%2 == 0.0
      update_spikes!(response,codinglayer,input,poisson_train)
      spikes_to_input!(input,codinglayer,poisson_train)
    end

    update_s!(unit1,unit2,network,timeparameters,transferparameters,input)
    decision = check_threshold(unit1,unit2,network)


    t +=timeparameters.dt
  end
  push!(rate1,unit1.rate)
  push!(rate2,unit2.rate)
  return unit1.meanrate,unit2.meanrate,response./input.μ0
end


end
