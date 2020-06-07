# Helper functions to display super ellipsoids/ellipsoids

include("geometry.jl")

using LinearAlgebra: norm

"""
    getSuperEllipsoidMesh(x0::AbstractVector, E::AbstractMatrix, p::Real)

    Return an array of points representing the super ellipsoid defined by norm(E*(x-x0),p) = 1
    E represents both the semi axes and the orientation of the super ellipsoid
"""
function getSuperEllipsoidMesh(x0::AbstractVector, E::AbstractMatrix, p::Real)
    # get a circle mesh
    θ = range(0,pi,length=50)
    ϕ = permutedims(range(0,2pi,length=100))

    v = spherical_to_cartesian.(1.0, θ, ϕ)

    # compute the value of the ellipsoid at the point and return the points on the surface
    r = norm.(Ref(E) .* v, p)
    ps = v ./ r

    # Add the offset
    ps = Ref(reshape(x0,3)) .+ ps

    return ps
end

using Rotations

"""
    getSuperEllipsoidMesh(x0::AbstractVector, E::AbstractMatrix, p::Real)

    Return an array of points representing the super ellipsoid defined by norm(D*R*(x-x0),p) = 1
    D represents the semi axes and R the orientation of the super ellipsoid
"""
function getSuperEllipsoidMesh(x0::AbstractVector, D::AbstractMatrix, R::RotMatrix, p::Real)
    return getSuperEllipsoidMesh(x0, D*R, p)
end

using Makie

"""
    displaySuperEllipsoid(x0::AbstractVector, D::AbstractMatrix, R::RotMatrix, p::Real=2; options...)

    Display a super ellipsoid represented by the equation norm(D*R*(x-x0),p) = 1, options are passed to the display
"""
function displaySuperEllipsoid(x0::AbstractVector, D::AbstractMatrix, R::RotMatrix, p::Real=2; options...)
    ps = getSuperEllipsoidMesh(x0, D, R, p)

    x, y, z = split_coords(ps)
    scene = surface(x, y, z; options...)
    return scene
end

"""
    displayEllipsoid(x0::AbstractVector, P::AbstractMatrix; options...)

    Display an ellipsoid represented by the equation (x-x0)^t * P * (x-x0) = 1, options are passed to the display
"""
function displayEllipsoid(x0::AbstractVector, P::AbstractMatrix; options...)
    D = Diagonal(sqrt.(eigvals(P)))
    R = RotMatrix(transpose(SMatrix{3,3}((eigvecs(P)))))
    return displaySuperEllipsoid(x0, D, R, 2; options...)
end
