using DiscreteVariableRepresentation, Printf

ng = 200
x0 = 0.0
l = 14.0
dx = l / ng
xg = x0 - l / 2 + dx / 2 .+ dx * [i for i in 0:ng-1]

dvr_obj = DVR(xg, v=DoubleWell(), mu=1.0)
e, u = solve(dvr_obj, 5)

e_sol = e
e_ref = [0.61742138, 0.75148846, 1.18262479, 1.2231278, 1.58168559]
err = abs.(e_sol - e_ref)
@assert all(err .< 1e-6)

import DiscreteVariableRepresentation: plot
import Plots: savefig
plot(dvr_obj, e, u, xmin=-3.5, xmax=3.5, ymin=-0.1, ymax=2.0)
savefig("../dvr-dw.png")
