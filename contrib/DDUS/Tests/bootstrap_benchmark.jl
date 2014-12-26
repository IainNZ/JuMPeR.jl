#Does a rough timing test of the bootstrapping calculation for the LCX Set

include("../helpers.jl")

function f1(boot_sample::Matrix, data::Matrix, numSamples::Int, a::Vector)
    Gamma::Float64 = 0.0
    for i = 1:numSamples
        randn!(a)
        a /= norm(a)
        const as = data * a
        const as_sample = boot_sample * a
        for b in as
            const lhs::Float64 = mean(max(as-b, 0)) - mean(max(as_sample-b, 0))   #bottleneck line
            if Gamma < lhs
                Gamma = lhs
            end
        end
    end
    Gamma
end


#Returns indx of last instance of val, assumes sorted
#0 if not a member
#Must have vector[start] >= val
function findlast_sort(val::Float64, vector::Vector; TOL::Float64=1e-10, start::Int64=1)
    i::Int64 = 1
    for i = int(start):length(vector)
        if abs(val - vector[i]) > TOL
            break
        end
    end
    i-1
end

function singlepass!(zetas::Vector, zetahats::Vector)
    sort!(zetas)  
    sort!(zetahats)

    vstar::Float64 = mean(zetas) - zetas[1] 
    vb::Float64    = mean(zetahats) - zetas[1]

    Gamma::Float64 = vstar - vb
    const N::Int64 = length(zetas)
    pbar::Float64    = 1.0 
    hat_indx::Int64  = 1
    hat_indx_::Int64 = 0
    
    for k = 2:length(zetas)
        vstar += (zetas[k-1] - zetas[k]) * (N-k+1)/N
        hat_indx = findlast_sort(zetas[k-1], zetahats, start=hat_indx_ + 1)
        pbar -=  (hat_indx - hat_indx_)/N
        hat_indx_ = hat_indx
        vb  += (zetas[k-1] - zetas[k]) * pbar
        Gamma = max(Gamma, vstar - vb)
    end
    Gamma
end

#a::Vector is used as storage between bootstraps
function f2(boot_sample::Matrix, data::Matrix, numSamples::Int, a::Vector)
    Gamma::Float64 = 0.
    for i = 1:numSamples
        randn!(a)
        a /= norm(a)
        Gamma = max(Gamma, singlepass!(data*a, boot_sample*a))
    end
    Gamma
end

#Approximates the threshold by sampling a bunch of abs for each bootstrap rep
#VG Could be optimized
function calc_ab_thresh(data::Matrix, delta::Float64, numBoots::Int, numSamples::Int)
    const N = size(data, 1)
    const d = size(data, 2)
    a::Vector{Float64} = zeros(Float64, d)
    boot(data, f2, 1-delta, numBoots, data, numSamples, a)
end

##################

data = randn(int(ARGS[1]), 2)
calc_ab_thresh(data, .1, 10, 10)  # one for the money!

Profile.init(10^6, .01)
tic()
@profile calc_ab_thresh(data, .1, 10, 1000)  #two for the show
println("Elapsed Time:\t", toc())
Profile.print()


#Baseline 47 s
