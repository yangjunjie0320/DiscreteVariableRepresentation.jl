using DiscreteVariableRepresentation, Printf

ng = 200
x0 = 0.0
l = 14.0
dx = l / ng
xg = x0 - l / 2 + dx / 2 .+ dx * [i for i in 0:ng-1]

dvr_obj = DVR(xg, v=SquareWell(x1=-5.0, x2=5.0, depth=4.5))
e, u = solve(dvr_obj, 5)

e_sol = e
e_ref = [0.04383565, 0.17523424, 0.39385678, 0.69908726, 1.08993682]
err = abs.(e_sol - e_ref)
@assert all(err .< 1e-6)

import DiscreteVariableRepresentation: plot
import Plots: savefig
plot(dvr_obj, e, u, xmin=-10.0, xmax=10.0, ymin=-0.25, ymax=2.0)
savefig("../dvr-sw.png")
