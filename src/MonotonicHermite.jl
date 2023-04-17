module MonotonicHermite

using SearchSortedNearest: searchsortedprevious

export MonotonicCubicHermite

struct MonotonicCubicHermite{TF<:Real, TT<:NTuple{4,TF}, TX<:AbstractVector{TF}}
    polys::Vector{TT}
    xs::TX
end

"""
    MonotonicCubicHermite(xs, ys, zs)

Compute a monotonic cubic spline Hermite interpolation from a sorted grid `xs`,
with corresponding function values `ys` and derivatives `zs`.

# Examples
```julia-repl
julia> xs = 10. .^ range(-2,0,6);
julia> ys = log.(xs);
julia> zs = inv.(xs);
julia> mch = MonotonicCubicHermite(xs,ys,zs);
```
"""
function MonotonicCubicHermite(
    xs::AbstractVector{TF},                                # must be sorted
    ys::AbstractVector{TF},
    zs::AbstractVector{TF}
    ) where {TF <: Real}
    nknots = length(xs)
    @assert nknots >= 2
    @assert nknots == length(ys) == length(zs)
    N = nknots - 1                                         # number of intervals
    C = Vector{NTuple{4,TF}}(undef, N)
    
    # init loop
    xₗ = first(xs)
    yₗ = first(ys)
    zₗ = first(zs)
    @inbounds for i in eachindex(C)                        # iterate intervals by index of leftmost knot
        xᵣ = xs[i+1]
        yᵣ = ys[i+1]
        zᵣ = zs[i+1]
        dx = xᵣ - xₗ                                       # length of current interval
        dx > zero(TF) || throw(ArgumentError("grid is not sorted: xs[$(i+1)] ≤ xs[$i]"))
        dy = yᵣ - yₗ
        S  = dy/dx                                         # slope in the current interval
        
        # adjust derivatives
        i == 1 && (zₗ = derivative_correction(zₗ, S, S))   # zₗ is already fixed (==zᵣ of previous step) except for the 1st iteration. uses the same slope for both (see first paragraph after eq 4.2)
        S₊ = i == N ? S : (ys[i+2] - yᵣ) / (xs[i+2] - xᵣ)  # slope of the next interval, default to current if at extreme (see first paragraph after eq 4.2)
        zᵣ = derivative_correction(zᵣ, S, S₊)

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

# derivative correction formula (eqs 4.2-4.3)
# note: z is derivative at knot whose interval to the left(right) has slope Sₗ(Sᵣ)
function derivative_correction(z::TF, Sₗ::TF, Sᵣ::TF) where {TF<:Real}
    σ  = sign(Sₗ) == sign(Sᵣ) ? sign(Sₗ) : sign(z)
    Sₘ = min(abs(Sₗ), abs(Sᵣ))
    return σ>0 ? min(max(zero(TF),z), 3Sₘ) : max(min(zero(TF),z), -3Sₘ)
end

"""
    (mch::MonotonicCubicHermite)(x)

Evaluate the `MonotonicCubicHermite` interpolator at a point `x`.

# Examples
```julia-repl
julia> xs = 10. .^ range(-2,0,6);
julia> ys = log.(xs);
julia> zs = inv.(xs);
julia> mch = MonotonicCubicHermite(xs,ys,zs);
julia> mch(0.5)
-0.6884921421564123
```
"""
function (mch::MonotonicCubicHermite)(x::Real)
    xs = mch.xs
    Base.isbetween(first(xs), x, last(xs)) || throw(ArgumentError("x is outside of the range."))
    i = searchsortedprevious(xs, x)
    evalpoly(x-xs[i], mch.polys[i])
end

end # module MonotonicHermite
