using MonotonicHermite, Plots

xs=cumsum([1e-2;rand(5)])
ys=log.(xs)
zs=inv.(xs)

mch = MonotonicCubicHermite(xs,ys,zs)

scatter(xs,ys)
plot!(x->mch(x))
