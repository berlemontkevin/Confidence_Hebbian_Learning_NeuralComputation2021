
#%% Premabule
using DrWatson
quickactivate("C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\","4-project_coding_categories")
# now can load Packages
using Distributed
using BSON
using DataFrames
include(srcdir("training_network.jl"))

#<<<<<<<


#%% Generate a set of old and new tuning curves at fixed parameters

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

function generate_TC(N,α)
    mycatGauss = MyCatGauss(width=[α,α])

    codinglayer = init_layer(N,mycatGauss)
    df = DataFrame(stim=Float64[],N=Int[],TC=Float64[],type=String[])
    stimlist = [0.001*i for i=0:1000]

    for j=1:length(codinglayer.centers)
             TC = TC_value(stimlist,codinglayer.centers[j],codinglayer.widths[j])
             for i=1:length(stimlist)
                 push!(df,[stimlist[i],j,TC[i],"Optimized"])
             end
    end
    codinglayer = init_layer(N)
    stimlist = [0.001*i for i=0:1000]

    for j=1:length(codinglayer.centers)
             TC = TC_value(stimlist,codinglayer.centers[j],codinglayer.widths[j])
             for i=1:length(stimlist)
                 push!(df,[stimlist[i],j,TC[i],"Uniform"])
             end
    end
    return df
end

#<<<<<<<<<<<<
N=20
df=generate_TC(N,α)
using CSV
CSV.write(datadir("notebooks","05-TuningCurves.csv"),df)


#%% Fisher Information

function generate_Fisher(N,α,stimlist)
    mycatGauss = MyCatGauss(width=[α,α])

    codinglayer = init_layer(N,mycatGauss)
    df = DataFrame(stim=Float64[],Delta_FCode=Float64[])
    codinglayerU = init_layer(N)

    FcodeO = zeros(length(stimlist))
    FcodeU = zeros(length(stimlist))

    for i=1:length(stimlist)
        x = stimlist[i]
    dummy =Fcode(x,codinglayer.centers,codinglayer.widths)
    dummy2 =Fcode(x,codinglayerU.centers,codinglayerU.widths)
    push!(df,[x,dummy-dummy2])

end

    return df
end

df = generate_Fisher(N,α,[0.001*i for i=0:1000])
CSV.write(datadir("notebooks","05-FisherInformation.csv"),df)



function generate_Cat(N,α,stimlist)
    mycatGauss = MyCatGauss(width=[α,α])
    df = DataFrame(stim=Float64[],prob_Cat=Float64[],Categorie=Int[])

    for x in stimlist
        cat1 = pxmu(x,mycatGauss.centers[1],mycatGauss.width[1])
        cat2 = pxmu(x,mycatGauss.centers[2],mycatGauss.width[2])
        push!(df,[x,cat1, 1])

        push!(df,[x,cat2,2])

    end


    return df
end

df = generate_Cat(N,α,[0.001*i for i=0:1000])
CSV.write(datadir("notebooks","05-Categories.csv"),df)

#<<<<<
