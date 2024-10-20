using DiscreteVariableRepresentation, Printf

ng = 200
x0 = 0.0
l = 14.0
dx = l / ng
xg = x0 - l / 2 + dx / 2 .+ dx * [i for i in 0:ng-1]

dvr_obj = DVR(xg, v=Morse(d=3.0, a=0.5, x0=0.0), mu=1.0)
e, u = solve(dvr_obj, 5)

e_sol = e
e_ref = [-2.41887756, -1.44413068, -0.71844164, -0.19817524, 0.3261174]
err = abs.(e_sol - e_ref)
@assert all(err .< 1e-6)

import DiscreteVariableRepresentation: plot
import Plots: savefig
plot(dvr_obj, e, u, xmin=-2.0, xmax=6.0, ymin=-3.6, ymax=0.5)
savefig("../dvr-morse.png")
