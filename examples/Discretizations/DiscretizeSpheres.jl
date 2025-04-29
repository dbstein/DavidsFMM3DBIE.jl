
using Revise
using DavidsFMM3DBIE
using Makie
using GLMakie

c1 = (0.0, 0.0, 0.0)
c2 = (1.0, 0.0, 0.0)
c3 = (0.0, 3.0, 0.0)

sphere1 = DiscretizeSphere(1.0, c1, 0.2)
sphere2 = DiscretizeSphere(0.5, c2, 0.1; order=4, panel_type=:QuadrangularGL)
sphere3 = DiscretizeSphere(2.0, c3, 0.125; order=8, panel_type=:QuadrangularChebyshev)

scatter(getCoordinates(sphere1), markersize=2, color=:blue)
scatter!(getCoordinates(sphere2), markersize=2, color=:red)
scatter!(getCoordinates(sphere3), markersize=2, color=:green)

