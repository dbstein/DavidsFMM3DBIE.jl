
using Revise
using DavidsFMM3DBIE
using Makie
using GLMakie

abc1 = (2.0, 0.5, 0.5)
abc2 = (0.5, 1.0, 2.0)
c1 = (0.0, 0.0, 0.0)
c2 = (1.0, 0.0, 0.0)

ellipsoid1 = DiscretizeEllipsoid(
                abc1,
                c1,
                0.1
            )
ellipsoid2 = DiscretizeEllipsoid(
                abc2,
                c2,
                0.0125;
                order=12,
                panel_type=:QuadrangularChebyshev
            )

scatter(getCoordinates(ellipsoid1), markersize=2, color=:blue)
scatter!(getCoordinates(ellipsoid2), markersize=2, color=:red)
