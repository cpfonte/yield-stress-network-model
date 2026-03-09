
include("SlitFlow1D.jl")

# Helper to read pressure at a node (either prescribed or unknown)
@inline function node_pressure(v, z, idx, p_in, p_out)
    if v in idx.left_nodes
        return p_in # left boundary node (p_in)
    elseif v in idx.right_nodes
        return p_out # right boundary node (p_out)
    else
        return z[idx.p_idx[v]]  # interior node
    end
end

# Residuals of nonlinear system F(z)
function F!(res, z, par)

    p_in, p_out, τ0, K, n, α, τS, β, R, idx, edges_df, nodes_df = par

    # 1) Mass conservation at unknown-pressure (interior) nodes
    for (k, v) in enumerate(idx.interior_nodes)
        s = 0.0
        for (e, sgn) in idx.inc[v]
            q = z[idx.q_idx[e]]
            s += sgn * q
        end
        res[k] = s
    end

    # 2) Edge pressure-flow equations
    offset = length(idx.interior_nodes)
    for e in 1:idx.nedges
        node_in = edges_df.node_in[e]
        node_out = edges_df.node_out[e]
        h  = edges_df.hmin[e]

        p1 = node_pressure(node_in, z, idx, p_in, p_out)
        p2 = node_pressure(node_out, z, idx, p_in, p_out)
        Δp_edge = p1 - p2

        q = z[idx.q_idx[e]]

        res[offset + e] = Δp_edge - Δp(q, h, τ0, K, n, α, τS, β, R)
    end

    return nothing
end