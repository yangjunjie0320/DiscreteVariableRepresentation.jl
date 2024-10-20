using DiscreteVariableRepresentation, Printf

ng = 200
x0 = 0.0
l = 14.0
dx = l / ng
xg = x0 - l / 2 + dx / 2 .+ dx * [i for i in 0:ng-1]
v = SimpleHarmonicOscillator(1.0, x0)
b = SincBasis()

dvr_obj = DVR(xg, 1.0, b, v)
e, u = solve(dvr_obj, 5)

@printf("e = [")
for i in 1:length(e)
    @printf("% 10.4f ", e[i])
end
@printf("]\n")

import DiscreteVariableRepresentation: plot
import Plots: savefig
plot(dvr_obj, e, u; xmin=-5.0, xmax=5.0, ymin=-0.1, ymax=10.0, nroots=20)
savefig("dvr.png")
