using DelimitedFiles
using StaticArrays

using QHull

"""
    getMesh(convex_hull::Chull)
    
    Small function return the points and the hull simplices in the correct display format.
    Then a direct call to mesh(points simplices) display the hull
"""
function getMesh(convex_hull::Chull)
    simplices = vcat(convex_hull.simplices...)
    return convex_hull.points, simplices
end
