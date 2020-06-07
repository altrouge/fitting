#
# Convert spherical coordinates to Cartesian coordinates.
# Note that θ,ϕ are named according to the "physicist convention"
# Note that using Vec3f0 here is efficient, but you could just use a
# normal array and everything would still work.
spherical_to_cartesian(r,θ,ϕ) = Vec3f0(r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ))

# Utility to split an array-of-length-3-arrays into x,y,z components
split_coords(ps) = getindex.(ps,1), getindex.(ps,2), getindex.(ps,3)
