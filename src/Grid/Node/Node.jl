Base.@kwdef mutable struct Node
    P = nothing
    N = nothing
    E = nothing
    F = nothing
end

function Node(Point, Pos)
  Node(; N = Pos, P = Point')
end
