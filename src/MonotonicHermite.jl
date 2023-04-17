module MonotonicHermite

using SearchSortedNearest

export MonotonicCubicHermite

struct MonotonicCubicHermite{TF<:Real, TT<:NTuple{4,TF}, TX<:AbstractVector{TF}}
    polys::Vector{TT}
    xs::TX
end

# implement Dougherty et al 1989
function MonotonicCubicHermite(
    xs::AbstractVector{TF}, # must be sorted
    ys::AbstractVector{TF},
    zs::AbstractVector{TF}
    ) where {TF <: Real}
    nknots = length(xs)
    @assert nknots >= 2
    @assert nknots == length(ys) == length(zs)
    N = nknots - 1 # number of intervals
    C = Vector{NTuple{4,TF}}(undef, N)
    
    # init loop
    xₗ = xs[1]
    yₗ = ys[1]
    zₗ = zs[1]
    for i in eachindex(C) # iterate intervals by index of leftmost knot
        xᵣ = xs[i+1]
        yᵣ = ys[i+1]
        zᵣ = zs[i+1]
        dx = xᵣ - xₗ      # length of current interval
        dx > zero(TF) || error("grid is not sorted: xs[$i] ≤ xs[$(i-1)]")
        dy = yᵣ - yₗ
        S  = dy/dx        # slope in the current interval

        # # interior derivative correction formula (eq 4.3)
        # if i<n
        #     σ  = sign(znew)
        #     sm = min(abs(s₋), abs(s₊))
        #     y' = σ>0 ? min(max(0,z),  3sm) : max(min(0,z), -3sm)
        # end

        # compute coefficients (eq 2.1)
        c₀   = yₗ
        c₁   = zₗ
        c₂   = (3S - zᵣ - 2zₗ) / dx
        c₃   = -(2S - zᵣ - zₗ) / (dx*dx)
        C[i] = (c₀, c₁, c₂, c₃)
        
        # update caches
        xₗ = xᵣ
        yₗ = yᵣ
        zₗ = zᵣ
    end
    MonotonicCubicHermite(C,xs)
end

function (mch::MonotonicCubicHermite)(x)
    xs = mch.xs
    Base.isbetween(first(xs), x, last(xs)) || throw(ArgumentError("x is outside of the range."))
    i = searchsortedprevious(xs, x)
    evalpoly(x-xs[i], mch.polys[i])
end

end # module MonotonicHermite
