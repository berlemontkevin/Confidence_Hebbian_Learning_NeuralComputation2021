using DrWatson



include(srcdir("structures.jl"))
include(srcdir("dynamical_system.jl"))
include(srcdir("information-theory.jl"))

### In this file we compute the optimal curves for the tuning cruves
## We add several possibility of normalization


function constant_gain(N::Int,cat::MyCatWeibull)
	# return a coding layer structure for a Weibull distribution of cat
	x = range(0.0,length=400,1.0)
	domaine = Float64[]
	for a in x
	push!(domaine,a)
	end
	Fcatlist = Float64[]
	for xi in x
				push!(Fcatlist,Fcat(xi,cat))
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


function constant_gain(N::Int,cat::MyCatGauss)
	# return a coding layer structure for a Weibull distribution of cat
	x = range(0.0,length=400,1.0)
	domaine = Float64[]
	for a in x
	push!(domaine,a)
	end
	Fcatlist = Float64[]
	for xi in x
				push!(Fcatlist,Fcat(xi,cat))
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

function gain_Hebbian(N::Int,cat::MyCatGauss)
	## optimal TCs when the gain is fixed via Hebbian learning
	x = range(0.0,length=400,1.0)
	domaine = Float64[]
	for a in x
		push!(domaine,a)
	end
	Fcatlist = Float64[]
	for xi in x
		push!(Fcatlist,Fcat(xi,cat))
	end

	g = zeros(length(domaine))
	for i=1:length(domaine)

		g[i] = abs(pmux(domaine[i],cat,1) - pmux(domaine[i],cat,2))
	end
	p = 1.0./400.0.*ones(400) # uniform prob. of stim

	g = g./sum(g.*p)


	dstemp = (p.*Fcatlist./g).^(1/3) #
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

	return mylayer,g
end

function gain_Hebbian(N::Int,cat::MyCatWeibull)
	## optimal TCs when the gain is fixed via Hebbian learning
	x = range(0.0,length=400,1.0)
	domaine = Float64[]
	for a in x
	push!(domaine,a)
	end
	Fcatlist = Float64[]
	for xi in x
				push!(Fcatlist,Fcat(xi,cat))
			end

	p = 1.0./400.0.*ones(400) # uniform prob. of stim

		g = zeros(length(domaine))
		for i=1:length(domaine)

			g[i] = abs(pmux(domaine[i],cat,1) - pmux(domaine[i],cat,2))
		end
		p = 1.0./400.0.*ones(400) # uniform prob. of stim

		g = g./sum(g.*p)
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

	return mylayer,g
end

function with_gain(N::Int,cats::MyCatGauss,normGain::Float64,IndiceNorme::Int)
  # compute the optimal tuning curves with constant gain
  x = range(0.0,length=400,1.0)
  domaine = Float64[]
  for a in x
    push!(domaine,a)
  end
  Fcatlist = Float64[]

  for xi in x
  push!(Fcatlist,Fcat(xi,cats)
  )
  end

  p = 1.0./400.0.*ones(400)
  dstemp = (Fcatlist).^(IndiceNorme/(3*IndiceNorme+1)).*(p.^((IndiceNorme+1)/(3*IndiceNorme+1)))
  ds = dstemp./(sum(dstemp.*1.0./400)).*(N-1)
  D = cumulative_integral(ds./400)
  width = 1.0./ds

  gstemp = (Fcatlist).^(IndiceNorme/(3*IndiceNorme+1)).*(p.^(-2/(3*IndiceNorme+1)))
  gs = gstemp./(sum(p.*(gstemp.^IndiceNorme).*1.0./400)).*(normGain^(1/IndiceNorme))
  G = gs

  centers = zeros(N)
  wf = zeros(N)
  gain = zeros(N)
  for i=1:N-1
      ind = minimum(findall(x->x>=i-0.005,D))
  centers[i+1],wf[i+1],gain[i+1] = domaine[ind],width[ind],G[ind]
  end
  wf[1]=wf[N]
  centers[1] = 1.0-centers[N]
  gain[1] = gain[N]

  mylayer = MyCodingLayer(centers,wf,zeros(N),zeros(N),gain)

  	return mylayer
end


function with_gain(N::Int,cats::MyCatWeibull,normGain::Float64,IndiceNorme::Int)
  # compute the optimal tuning curves with constant gain
  x = range(0.0,length=400,1.0)
  domaine = Float64[]
  for a in x
    push!(domaine,a)
  end
  Fcatlist = Float64[]

  for xi in x
  push!(Fcatlist,Fcat(xi,cats)
  )
  end

  p = 1.0./400.0.*ones(400)
  dstemp = (Fcatlist).^(IndiceNorme/(3*IndiceNorme+1)).*(p.^((IndiceNorme+1)/(3*IndiceNorme+1)))
  ds = dstemp./(sum(dstemp.*1.0./400)).*(N-1)
  D = cumulative_integral(ds./400)
  width = 1.0./ds

  gstemp = (Fcatlist).^(IndiceNorme/(3*IndiceNorme+1)).*(p.^(-2/(3*IndiceNorme+1)))
	gs = gstemp./(sum(p.*(gstemp.^IndiceNorme).*1.0./400)).*(normGain^(1/IndiceNorme))
  G = gs

  centers = zeros(N)
  wf = zeros(N)
  gain = zeros(N)
  for i=1:N-1
      ind = minimum(findall(x->x>=i-0.005,D))
  centers[i+1],wf[i+1],gain[i+1] = domaine[ind],width[ind],G[ind]
  end
  wf[1]=wf[N]
  centers[1] = 1.0-centers[N]
  gain[1] = gain[N]

  mylayer = MyCodingLayer(centers,wf,zeros(N),zeros(N),gain)

  	return mylayer
end
