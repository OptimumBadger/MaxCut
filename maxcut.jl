using JuMP, Gurobi, SCS, LinearAlgebra

# =============================================================================
# Graph Reading
# =============================================================================

"""
    read_graph(filename) -> (n, edges)

Reads a Max-Cut instance file. Supports both .rud and .mc formats.
Format: first line is "n_vertices n_edges", then each subsequent line is "i j w".
Returns number of vertices n and a vector of (i, j, w) edge tuples.
"""
function read_graph(filename::String)
    open(filename) do f
        parts = split(strip(readline(f)))
        n = parse(Int, parts[1])
        m = parse(Int, parts[2])
        edges = Tuple{Int,Int,Float64}[]
        for _ in 1:m
            row = split(strip(readline(f)))
            i = parse(Int, row[1])
            j = parse(Int, row[2])
            w = parse(Float64, row[3])
            push!(edges, (i, j, w))
        end
        return n, edges
    end
end

# =============================================================================
# SDP Relaxation  (C^SDP)
# =============================================================================

"""
    solve_sdp(n, edges) -> objective_value

Solves the SDP relaxation of Max-Cut:

    max  Σ_{ij ∈ E} w_ij (x_i + x_j - 2 X_ij)
    s.t. X ⪰ 0                  (n×n PSD matrix)
         X_ii = x_i  ∀ i        (diagonal linking constraint)
         x ∈ [0,1]^n

Uses SCS solver (supports semidefinite constraints).
"""
function solve_sdp(n::Int, edges::Vector{Tuple{Int,Int,Float64}})
    model = Model(SCS.Optimizer)
    #set_silent(model)

    @variable(model, 0 <= x[1:n] <= 1)

    # Lifted (n+1)×(n+1) PSD matrix M = [1, x'; x, X]
    # Schur complement forces X - xx' ⪰ 0, preventing negative X_ij
    @variable(model, M[1:n+1, 1:n+1], Symmetric)
    @constraint(model, M in PSDCone())

    # Top-left entry = 1
    @constraint(model, M[1, 1] == 1)

    # First row/col = x
    for i in 1:n
        @constraint(model, M[1, i+1] == x[i])
    end

    # Diagonal of X block = x_i
    for i in 1:n
        @constraint(model, M[i+1, i+1] == x[i])
    end

    # M[i+1, j+1] plays the role of X_ij
    @objective(model, Max,
        sum(w * (x[i] + x[j] - 2 * M[i+1, j+1]) for (i, j, w) in edges))

    optimize!(model)
    return objective_value(model)
end

# =============================================================================
# McCormick + Triangle Relaxation  (C^{MC+Tri})
# =============================================================================

"""
    solve_mc_tri(n, edges) -> objective_value

Solves the LP relaxation of Max-Cut with McCormick and triangle inequalities:

    max  Σ_{ij ∈ E} w_ij (x_i + x_j - 2 X_ij)
    s.t. McCormick inequalities for all pairs (i,j):
             X_ij ≥ 0
             X_ij ≥ x_i + x_j - 1
             X_ij ≤ x_i
             X_ij ≤ x_j
         Triangle inequalities for all triples (i,j,k):
             X_ij + X_ik ≤ x_i + X_jk
             X_ij + X_jk ≤ x_j + X_ik
             X_ik + X_jk ≤ x_k + X_ij
             x_i + x_j + x_k - X_ij - X_ik - X_jk ≤ 1
         x ∈ [0,1]^n,  X_ij ∈ [0,1]

NOTE: X variables are needed for ALL pairs (not just edges) due to triangle
inequalities. For large n this becomes O(n^3) constraints — use only on
small-to-medium instances (n ≲ 300).

Uses Gurobi solver (pure LP).
"""
function solve_mc_tri(n::Int, edges::Vector{Tuple{Int,Int,Float64}})
    model = Model(Gurobi.Optimizer)
    #set_silent(model)

    # Build edge set and adjacency list for triangle detection
    edge_set = Set((min(i,j), max(i,j)) for (i,j,_) in edges)
    adj = [Set{Int}() for _ in 1:n]
    for (i, j, _) in edges
        push!(adj[i], j)
        push!(adj[j], i)
    end

    @variable(model, 0 <= x[1:n] <= 1)
    # X variables only for edges (i,j) ∈ E
    @variable(model, 0 <= X[e in edge_set] <= 1)

    # Helper: retrieve X variable for edge (i,j)
    Xvar(i, j) = X[(min(i,j), max(i,j))]

    # McCormick inequalities only for edges
    for (i, j, _) in edges
        @constraint(model, Xvar(i,j) >= x[i] + x[j] - 1)
        @constraint(model, Xvar(i,j) <= x[i])
        @constraint(model, Xvar(i,j) <= x[j])
    end

    # Triangle inequalities only for triangles (i,j,k) where all 3 pairs are edges
    triangles_seen = Set{Tuple{Int,Int,Int}}()
    for (i, j, _) in edges
        for k in intersect(adj[i], adj[j])
            t = Tuple(sort([i, j, k]))
            t in triangles_seen && continue
            push!(triangles_seen, t)
            a, b, c = t
            @constraint(model, Xvar(a,b) + Xvar(a,c) <= x[a] + Xvar(b,c))
            @constraint(model, Xvar(a,b) + Xvar(b,c) <= x[b] + Xvar(a,c))
            @constraint(model, Xvar(a,c) + Xvar(b,c) <= x[c] + Xvar(a,b))
            @constraint(model, x[a] + x[b] + x[c] - Xvar(a,b) - Xvar(a,c) - Xvar(b,c) <= 1)
        end
    end

    @objective(model, Max,
        sum(w * (x[i] + x[j] - 2 * Xvar(i, j)) for (i, j, w) in edges))

    optimize!(model)
    return objective_value(model)
end

# =============================================================================
# McCormick + Triangle + SDP Relaxation  (C^{MC+Tri+SDP})
# =============================================================================

"""
    solve_mc_tri_sdp(n, edges) -> objective_value

Solves the combined SDP + McCormick + Triangle relaxation of Max-Cut:

    max  Σ_{ij ∈ E} w_ij (x_i + x_j - 2 X_ij)
    s.t. X ⪰ 0
         X_ii = x_i  ∀ i
         McCormick inequalities for all pairs (i,j)
         Triangle inequalities for all triples (i,j,k)
         x ∈ [0,1]^n

NOTE: Same O(n^3) constraint scaling as solve_mc_tri. Practical for n ≲ 200.

Uses SCS solver.
"""
function solve_mc_tri_sdp(n::Int, edges::Vector{Tuple{Int,Int,Float64}})
    model = Model(SCS.Optimizer)
    #set_silent(model)

    # Build adjacency list for triangle detection
    adj = [Set{Int}() for _ in 1:n]
    for (i, j, _) in edges
        push!(adj[i], j)
        push!(adj[j], i)
    end

    @variable(model, 0 <= x[1:n] <= 1)

    # Lifted (n+1)×(n+1) PSD matrix M = [1, x'; x, X]
    @variable(model, M[1:n+1, 1:n+1], Symmetric)
    @constraint(model, M in PSDCone())

    @constraint(model, M[1, 1] == 1)
    for i in 1:n
        @constraint(model, M[1, i+1] == x[i])
        @constraint(model, M[i+1, i+1] == x[i])
    end

    # Helper: X_ij = M[i+1, j+1]
    Xij(i, j) = M[i+1, j+1]

    # McCormick inequalities only for edges
    for (i, j, _) in edges
        @constraint(model, Xij(i,j) >= 0)
        @constraint(model, Xij(i,j) >= x[i] + x[j] - 1)
        @constraint(model, Xij(i,j) <= x[i])
        @constraint(model, Xij(i,j) <= x[j])
    end

    # Triangle inequalities only for triangles (i,j,k) where all 3 pairs are edges
    triangles_seen = Set{Tuple{Int,Int,Int}}()
    for (i, j, _) in edges
        for k in intersect(adj[i], adj[j])
            t = Tuple(sort([i, j, k]))
            t in triangles_seen && continue
            push!(triangles_seen, t)
            a, b, c = t
            @constraint(model, Xij(a,b) + Xij(a,c) <= x[a] + Xij(b,c))
            @constraint(model, Xij(a,b) + Xij(b,c) <= x[b] + Xij(a,c))
            @constraint(model, Xij(a,c) + Xij(b,c) <= x[c] + Xij(a,b))
            @constraint(model, x[a] + x[b] + x[c] - Xij(a,b) - Xij(a,c) - Xij(b,c) <= 1)
        end
    end

    @objective(model, Max,
        sum(w * (x[i] + x[j] - 2 * Xij(i,j)) for (i, j, w) in edges))

    optimize!(model)
    return objective_value(model)
end

# =============================================================================
# CLI Entry Point
# Usage: julia maxcut.jl <instance_file> <formulation>
#   formulation: 1 = SDP, 2 = MC+Tri, 3 = MC+Tri+SDP
# =============================================================================

if length(ARGS) != 2
    println("Usage: julia maxcut.jl <instance_file> <formulation>")
    println("  formulation: 1 = SDP | 2 = MC+Tri | 3 = MC+Tri+SDP")
    exit(1)
end

instance_file = ARGS[1]
formulation   = parse(Int, ARGS[2])

if !isfile(instance_file)
    println("Error: file '$instance_file' not found.")
    exit(1)
end

if formulation ∉ (1, 2, 3)
    println("Error: formulation must be 1, 2, or 3.")
    exit(1)
end

n, edges = read_graph(instance_file)
println("Instance: $instance_file  |  Nodes: $n  |  Edges: $(length(edges))")

if formulation == 1
    println("Formulation: SDP")
    val = solve_sdp(n, edges)
    println("SDP bound: $val")
elseif formulation == 2
    println("Formulation: MC+Tri")
    val = solve_mc_tri(n, edges)
    println("MC+Tri bound: $val")
else
    println("Formulation: MC+Tri+SDP")
    val = solve_mc_tri_sdp(n, edges)
    println("MC+Tri+SDP bound: $val")
end
