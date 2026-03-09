using CSV, DataFrames, CairoMakie, DelaunayTriangulation
const DT = DelaunayTriangulation

cd(@__DIR__)

function net_gen()

    data = CSV.read(net_name, DataFrame)
    
    vecdata = [(data[i,1], data[i,2]) for i in 1:size(data,1)]
    tri  = triangulate(vecdata)

    clip_points   = [(0.0, 0.0), (L, 0.0), (L, L), (0.0, L)]
    clip_vertices = (1, 2, 3, 4, 1)
    clip_polygon  = (clip_points, clip_vertices)

    vorn = voronoi(tri; clip=true, clip_polygon=clip_polygon, predicates = FastKernel())

    nodes_df, edges_df = voronoi_nodes_and_edges(vorn; keep_boundary=false)

    idx = build_network_indices(nodes_df, edges_df; L)

    fig1 = Figure()

    ax1 = Axis(
        fig1[1, 1]; 
        aspect = DataAspect(),
    )

    voronoiplot!(
        ax1,
        vorn; 
        show_generators=true
    )

    scatter!(
        ax1,
        data[!,1],
        data[!,2],
        marker = Circle,
        color = :red,
        markersize = 2.0*R,
        markerspace = :data,
        alpha = 0.25,
    )

    triplot!(ax1, tri; strokecolor=:red, strokewidth=1)

    fig2 = Figure()

    ax2 = Axis(
        fig2[1, 1]; 
        aspect = DataAspect(),
        limits = ((-0.05L, 1.05L), (-0.05L, 1.05L)),
        xgridvisible = false,
        ygridvisible = false,
    )

    scatter!(
        ax2,
        data[!,1],
        data[!,2],
        marker = Circle,
        color = :red,
        markersize = 2.0*R,
        markerspace = :data,
        alpha = 0.25,
    )

    used = falses(nrow(nodes_df))

    for r in eachrow(edges_df)
        x = [nodes_df.x[r.node_in], nodes_df.x[r.node_out]]
        y = [nodes_df.y[r.node_in], nodes_df.y[r.node_out]]
        if maximum(x) == L || minimum(x) == 0.0
            lines!(
                ax2,
                x, 
                y,
                color = :green,
                linewidth = 1.5,
            )
        else
            lines!(
                ax2,
                x, 
                y,
                color = :blue,
                linewidth = 1.5,
            )
        end
        used[r.node_in] = true
        used[r.node_out] = true
    end

    ids = findall(used)

    text!(
        ax2,
        string.(ids),
        position = [(nodes_df.x[i], nodes_df.y[i]) for i in ids],
        align = (:left, :bottom),
        fontsize = Int(round(10*10/L, digits=0)),
    )

    lines!(
        ax2,
        [clip_points; clip_points[1]],
        color = :black,
        linewidth = 2.0,
    )

    display(fig1)
    display(fig2)
    save("voronoi.pdf", fig1, pdf_version="1.4")
    save("network.pdf", fig2, pdf_version="1.4")

    save("voronoi.png", fig1)
    save("network.png", fig2)

    return (nodes_df, edges_df, idx, data)

end

function voronoi_nodes_and_edges(vorn; keep_boundary=false)
    
    nodes = DT.get_polygon_points(vorn)

    # Dict{(u,v) => cell_id}, oriented Voronoi edge -> adjacent cell (generator index)
    adj = DT.get_adjacent(DT.get_adjacent(vorn))
    
    # Nodes table (IDs align with indices used in adj keys)
    nodes_df = DataFrame(
        node_id = collect(1:length(nodes)),
        x = first.(nodes),
        y = last.(nodes),
    )

    # Build unique undirected edges
    seen = Set{Tuple{Int,Int}}()
    v1s = Int[]; v2s = Int[]
    c1s = Int[]; c2s = Int[]
    dls = Union{Float64,Missing}[]

    for ((u, v), _) in adj
        a, b = u < v ? (u, v) : (v, u)   # canonical undirected edge key
        (a, b) in seen && continue
        push!(seen, (a, b))

        cell1 = get(adj, (a, b), 0)
        cell2 = get(adj, (b, a), 0)

        # Boundary edges (from clipping) usually have only one adjacent cell, no dual Delaunay edge
        if cell1 == 0 || cell2 == 0 || cell1 == cell2
            keep_boundary || continue
            push!(v1s, a); push!(v2s, b)
            push!(c1s, cell1); push!(c2s, cell2)
            push!(dls, missing)
            continue
        end

        # Dual Delaunay edge length = distance between the two generators (sites)
        p1 = DT.get_generator(vorn, cell1)
        p2 = DT.get_generator(vorn, cell2)
        delaunay_len = hypot(p1[1] - p2[1], p1[2] - p2[2])
        hmin = delaunay_len - 2.0*R 

        push!(v1s, a); push!(v2s, b)
        push!(c1s, cell1); push!(c2s, cell2)
        push!(dls, delaunay_len)
    end

    # Filter out edges that touch the boundary (y=0 or y=L) since they do not represent valid throats between pores
    eids = Int(0)
    v1=Int[]; v2=Int[]; hmin=Float64[]
    for i in 1:length(v1s)
        y = [nodes_df.y[v1s[i]], nodes_df.y[v2s[i]]]
        if maximum(y) == L || minimum(y) == 0.0
            continue
        end
        eids += 1
        push!(v1, v1s[i])
        push!(v2, v2s[i])
        push!(hmin, dls[i]/2.0 - R)
    end

    edges_df = DataFrame(
        edge_ID = 1:length(v1),
        node_in = v1,
        node_out = v2,
        hmin = round.(hmin, digits=4),
    )

    return nodes_df, edges_df
end

function build_network_indices(nodes_df, edges_df; L)

    nnodes = nrow(nodes_df)
    nedges = nrow(edges_df)

    # boundary classification (only left/right are prescribed)
    left  = falses(nnodes)
    right = falses(nnodes)
    for v in 1:nnodes
        left[v]  = nodes_df.x[v] .== 0.0
        right[v] = nodes_df.x[v] .== L
    end

    # only nodes that appear in edges_df matter
    used = unique(vcat(edges_df.node_in, edges_df.node_out))

    left_nodes  = sort!(collect(filter(v -> left[v], used)))
    right_nodes = sort!(collect(filter(v -> right[v], used)))
    interior_nodes = sort!(collect(setdiff(used, union(left_nodes, right_nodes))))

    q_idx = collect(1:nedges)
    p_idx = zeros(Int, nnodes)
    for (k, v) in enumerate(interior_nodes)
        p_idx[v] = nedges + k
    end

    # Incidence list for mass conservation:
    # define Q positive from node_in -> node_out
    inc = [Tuple{Int,Int}[] for _ in 1:nnodes]  # per node: (edge_row, sign)
    for e in 1:nedges
        node_in = edges_df.node_in[e]
        node_out = edges_df.node_out[e]
        push!(inc[node_in], (e, +1))  # leaving node_in
        push!(inc[node_out], (e, -1))  # entering node_out
    end

    is_left  = nodes_df.x .== 0.0
    is_right = nodes_df.x .== L

    return(
        left_nodes = left_nodes,
        right_nodes = right_nodes,
        interior_nodes = interior_nodes,
        q_idx = q_idx,
        p_idx = p_idx,
        inc = inc,
        nnodes = nnodes,
        nedges = nedges,
        is_left = is_left,
        is_right = is_right,
        n_unknown = nedges + length(interior_nodes)
    )
end

(nodes_df, edges_df, idx, data) = net_gen();

