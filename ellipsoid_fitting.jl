using DataFrames
using CSV: read
#using Convex: Semidefinite, tr, minimize, solve
using Convex
using SCS
using DelimitedFiles
using Rotations
using LinearAlgebra
using StaticArrays

"""
    findMinEllipsoid(points)

    Find the minimum volume ellipsoid containing a set of points
    return the min ellipsoid represented by
    (x-x0)^t P (x-x0) = 1
"""
function findMinEllipsoid(points)
    # solve norm(P*x - x0, 2) <= 1 for all points p
    # then, the volume of the ellipsoid is 4/3 pi / det(X)
    # minimizing the volume is equivalent to minimizing -logdet(P) which is convex
    P = Semidefinite(3)
    x0 = Variable(3)
    problem = minimize(-logdet(P))
    problem.constraints += [P >= 0]

    for row in eachrow(points)
        x = Vector(row)
        problem.constraints += [norm(P*x-x0,2) <= 1]
    end
    solver = () -> SCS.Optimizer(verbose=true)
    solve!(problem, solver)

    # return x0 and P in the correct representation
    return reshape(inv(P.value)*x0.value,3), P.value^2
end

include("qhull.jl")
include("ellipsoid.jl")

ellipsoidVolume(P::AbstractMatrix) = 4/3 * pi / det(P)

"""
    findMinEllipsoidF(file; display_convex_hull=false, display_points=false)

    Find the minimum volume ellipsoid containing a set of points in a text file.
    Display the ellipsoid, and if queried, the set of points or the convex hull.

    Return the scene
"""
function findMinEllipsoidF(file; display_convex_hull=true, display_points=true)
    df = readdlm(file, ' ', Float64)
    x0, P = findMinEllipsoid(df)
    println("Min volume: ", ellipsoidVolume(P), "m^3")
    scene = displayEllipsoid(x0, P; alpha=0.1, transparency=true)

    if display_convex_hull
        ch = chull(df)
        cpoints, csimplices = getMesh(ch)
        mesh!(scene, cpoints, csimplices)
    end
    if display_points
        ox, oy, oz = df[:,1], df[:,2], df[:,3]
        scatter!(scene, ox, oy, oz, markersize=0.005)
    end
    return scene
end
