using DiscreteVariableRepresentation, Printf

ng = 200
x0 = 0.0
l = 14.0
dx = l / ng
xg = x0 - l / 2 + dx / 2 .+ dx * [i for i in 0:ng-1]

dvr_obj = DVR(xg, v=SimpleHarmonicOscillator(k=1.0, x0=0.0), mu=1.0)
e, u = solve(dvr_obj, 5)

e_sol = e
e_ref = [0.5, 1.5, 2.5, 3.5, 4.5]
err = abs.(e_sol - e_ref)
@assert all(err .< 1e-6)

import DiscreteVariableRepresentation: plot
import Plots: savefig
plot(dvr_obj, e, u, xmin=-3.5, xmax=3.5, ymin=-0.25, ymax=6.0)
savefig("../dvr-sho.png")
