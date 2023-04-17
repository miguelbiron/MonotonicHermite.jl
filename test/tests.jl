using MonotonicHermite, Plots

xs=cumsum([1e-2;rand(5)])
ys=log.(xs)
zs=inv.(xs)

mch = MonotonicCubicHermite(xs,ys,zs)

scatter(xs,ys)
plot!(x->mch(x))

using Distributions

d  = Beta(1/2,1/2)
xs = 10. .^ range(-3,-1)
xs = [0.;xs;reverse(1 .- xs);1.]
ys = cdf(d, xs)
zs = map(x->pdf(d,x), xs)
mch = MonotonicCubicHermite(xs,ys,zs)
scatter(xs,ys)
plot!(x->mch(x))
