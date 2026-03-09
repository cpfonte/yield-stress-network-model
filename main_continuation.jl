using NonlinearSolve, LinearSolve, LinearAlgebra, CSV, DataFrames, CairoMakie, Printf

cd(@__DIR__)

const p_out = 0.0                        # Outlet pressure

# Geometric parameters
const R = 1.0                            # Pillar radius
const L = 10.0                           # Domain length

# CSV file containing the coordinates of the pillar centers
const net_name = "sample_porous_medium.csv"   

# Rheological parameters    
const τ0 = 1.0                           # Yield stress
const K = 1.0                            # Consistency index
const n = 1.0                            # Flow index

# Slip parameters
const α = 0.1*R/K                        # Slip coefficient
const β = 1.0                            # Slip index
const τS = 0.2*τ0                        # Slip yield stress

# Nonlinear solver parameters
const abstol = 1e-8
const reltol = 1e-4
const maxiters = 50

include("NetworkGeneration.jl")
include("Residuals.jl")

# Number of unknowns in the system (pressure at interior nodes + flow rates in edges)
const Neq = idx.n_unknown

# Show iteration progress
show_iter = Val(true)  # Set to Val(true) to show iteration progress, Val(false) to hide

# Show flow plot for each solution
show_flow = Val(true)  # Set to Val(true) to show flow plot for each solution, Val(false) to hide

function main()

    # We use a continuation method to sweep through a range of inlet pressures. This approach helps to ensure convergence of the nonlinear solver, especially near the yield point where the solution can change rapidly.

    # Range of inlet pressures
    lo, hi = 5.0, 2000.0
    # Number of points in the sweep
    N = 50
    # We use a non-uniform spacing for the inlet pressures, with more points clustered near the lower end of the range (close to the yield point). The parameter γ controls how strongly the points are clustered towards the lower end; a larger γ results in more clustering.
    γ = 3.0
    t = range(0.0, 1.0; length = N)
    p_in = lo .+ (hi - lo) .* (1 .- t) .^ γ

    # Run the sweep of inlet pressures and collect the corresponding flow rates at the inlet and outlet
    Q_in, Q_out = run_sweep(p_in)

    # Plotting the results - Pressure gradient vs. Bingham number
    U = Q_out/L
    B = τ0 ./ (K .* (U ./R) .^ n)
    G = p_in ./ L .* R ./ (K .* (U ./R) .^ n)

    fig = Figure()

    ax = Axis(
        fig[1, 1]; 
        limits = ((1e-1, 1e4), (1e1, 1e4)),
        xgridvisible = true,
        ygridvisible = true,
        xscale = log10,
        yscale = log10,
        xlabel = "Bingham number",
        ylabel = "Pressure gradient",
    )

    # Some of the points may be non-finite due to numerical issues, so we filter those out before plotting the line connecting the points.
    good = findall(isfinite.(B) .& isfinite.(G))
    isempty(good) || lines!(ax, B[good], G[good], color=:red)

    scatter!(
        ax,
        B,
        G,
        color = :blue,
        markersize = 8.0,
        marker = Circle,
    )

    display(fig)
    save("G_vs_B.pdf", fig, pdf_version="1.4")
    save("G_vs_B.png", fig)

    return nothing
end

function run_sweep(p_in)

    # The function run_sweep takes a vector of inlet pressures and runs the nonlinear solver for each pressure to compute the corresponding flow rates at the inlet and outlet. It uses a continuation method, where the solution from the previous pressure is used as the initial guess for the next pressure, which helps with convergence.

    par  = (p_in[1], p_out, 0.0, K, n, α, 0.0, β, R, idx, edges_df, nodes_df)

    # We initiate the solution at the first inlet pressure with a simple initial guess (e.g., all pressures set to 1.0 and all flow rates set to 1.0). This seems to be good enough to get convergence for the first point, and then we can use the solution from this point as the initial guess for the next point in the sweep.
    prob = NonlinearProblem(F!, fill(1.0, Neq), par)

    @time sol = solve(
        prob,
        alg = NewtonRaphson( 
            autodiff   = AutoFiniteDiff(fdtype = Val(:central)),
        );
        abstol = abstol,
        reltol = reltol,
        maxiters = maxiters,
        show_trace = show_iter,
    )

    # Initialize arrays to store flow rates at inlet and outlet for each inlet pressure
    Q_in = fill(NaN, length(p_in))
    Q_out = fill(NaN, length(p_in))

    # Use solution from initial solve as initial guess for the pressure sweep
    z_p = sol.u
    z_pp = sol.u
    p_p = p_in[2]
    p_pp = p_in[1]

    for i in eachindex(p_in)

        p  = (p_in[i], p_out, τ0, K, n, α, τS, β, R, idx, edges_df, nodes_df)

        # For the first two points, we use the previous solution as the initial guess. After that, we use a linear extrapolation based on the last two solutions to get a better initial guess, which can help convergence especially near the yield point where the solution changes rapidly.
        if i <= 2
            z0 = z_p
        else
            z0 = z_p + (z_p .- z_pp) * (p_in[i] - p_p) / (p_p - p_pp)
            ω = 1.0  # Relaxation factor for the extrapolation, can be tuned if needed to improve convergence
            z0 = z_p .+ ω .* (z0 .- z_p)
        end

        prob = NonlinearProblem(F!, z0, p)

        println("\nSolving for inlet pressure p_in = ", p_in[i])
            
        @time sol = solve(
            prob,
            # We use a trust-region method for the nonlinear solver, which is generally more robust than Newton-Raphson, especially for problems with strong nonlinearity or near singularities (like near the yield point). The autodiff option allows the solver to compute the Jacobian using finite differences, which is convenient and often efficient for problems of this size.
            alg = TrustRegion(
                        autodiff   = AutoFiniteDiff(fdtype = Val(:central)),
                    );
            abstol = abstol,
            reltol = reltol,
            maxiters = maxiters,
            show_trace = show_iter,
        )

        println()
        println("Final summary")
        println("─────────────")
        println("  Status                : $(sol.retcode)")
        println("  Residual norm ‖F(u*)‖₂: ", @sprintf("%.3e", norm(sol.resid)))

        # Compute total flow rate at inlet and outlet
        if SciMLBase.successful_retcode(sol) == true

            Q_in[i] = net_left(sol.u, idx, edges_df)   # net inflow at x==0 into the domain
            Q_out[i] = net_right(sol.u, idx, edges_df)  # net outflow at x==L out of the domain
            # Update previous solutions for next iteration
            z_pp = z_p
            z_p = sol.u
            p_pp = p_p
            p_p = p_in[i]
            CSV.write("output.csv", DataFrame(p_in = p_in, Q_in = Q_in, Q_out = Q_out))

            if show_flow == Val(true)
                plot_flow(sol, idx, edges_df, nodes_df)
            end

        else

            # Run simulation agaun with different method if it failed
            println("\nSimulation failed for p_in = ", p_in[i], ", trying with different method.")

            @time sol = solve(
                prob,
                alg = RobustMultiNewton(
                        autodiff   = AutoFiniteDiff(fdtype = Val(:central)),
                    );
                abstol = abstol,
                reltol = reltol,
                maxiters = maxiters,
                show_trace = show_iter,
            )

            println()
            println("Final summary")
            println("─────────────")
            println("  Status                : $(sol.retcode)")
            println("  Residual norm ‖F(u*)‖₂: ", @sprintf("%.3e", norm(sol.resid)))

            if SciMLBase.successful_retcode(sol) == true

                Q_in[i] = net_left(sol.u, idx, edges_df)   # net inflow at x==0 into the domain
                Q_out[i] = net_right(sol.u, idx, edges_df)  # net outflow at x==L out of the domain
                # Update previous solutions for next iteration
                z_pp = z_p
                z_p = sol.u
                p_pp = p_p
                p_p = p_in[i]
                CSV.write("output.csv", DataFrame(p_in = p_in, Q_in = Q_in, Q_out = Q_out))

                if show_flow == Val(true)
                    plot_flow(sol, idx, edges_df, nodes_df)
                end

            else
                
                # Run the simulation again with different initial guess if it failed
                println("\nSimulation failed for p_in = ", p_in[i], ", trying with different initial guess.")
                
                # First solve without yield stress and slip yield stress to get a good initial guess
                p  = (p_in[i], p_out, 0.0, K, n, α, 0.0, β, R, idx, edges_df, nodes_df)
                prob = NonlinearProblem(F!, fill(1.0, Neq), p)

                @time sol = solve(
                    prob,
                    alg = NewtonRaphson( 
                        autodiff   = AutoFiniteDiff(fdtype = Val(:central)),
                    );
                    abstol = abstol,
                    reltol = reltol,
                    maxiters = maxiters,
                    show_trace = show_iter,
                )

                z0 = sol.u

                # Now solve the full problem again
                p  = (p_in[i], p_out, τ0, K, n, α, τS, β, R, idx, edges_df, nodes_df)

                prob = NonlinearProblem(F!, z0, p)

                @time sol = solve(
                    prob,
                    alg = TrustRegion(
                        autodiff   = AutoFiniteDiff(fdtype = Val(:central)),
                    );
                    abstol = abstol,
                    reltol = reltol,
                    maxiters = maxiters,
                    show_trace = show_iter,
                )

                println()
                println("Final summary")
                println("─────────────")
                println("  Status                : $(sol.retcode)")
                println("  Residual norm ‖F(u*)‖₂: ", @sprintf("%.3e", norm(sol.resid)))

                if SciMLBase.successful_retcode(sol) == true
                    Q_in[i] = net_left(sol.u, idx, edges_df)   # net inflow at x==0 into the domain
                    Q_out[i] = net_right(sol.u, idx, edges_df)  # net outflow at x==L out of the domain
                    z_pp = z_p
                    z_p = sol.u
                    p_pp = p_p
                    p_p = p_in[i]
                    CSV.write("output.csv", DataFrame(p_in = p_in, Q_in = Q_in, Q_out = Q_out))

                    if show_flow == Val(true)
                        plot_flow(sol, idx, edges_df, nodes_df)
                    end

                else
                    println("\nSimulation failed for p_in = ", p_in[i])
                    println("Skipping this pressure. Continuing the sweep.")
                end
            end
        end
    end

    return Q_in, Q_out

end

# Calculate net flow across left boundary (should equal Qin)
function net_left(z, idx, edges_df)
    Q = 0.0
    for e in 1:idx.nedges
        a, b = edges_df.node_in[e], edges_df.node_out[e]
        la, lb = idx.is_left[a], idx.is_left[b]
        (la ⊻ lb) || continue              # only edges that cross the left boundary
        q = z[idx.q_idx[e]]
        Q += la ? q : -q                   # + if boundary endpoint is node_in, - if node_out
    end
    return Q
end

# Calculate net flow across right boundary (should equal Qout)
function net_right(z, idx, edges_df)
    Q = 0.0
    for e in 1:idx.nedges
        a, b = edges_df.node_in[e], edges_df.node_out[e]
        ra, rb = idx.is_right[a], idx.is_right[b]
        (ra ⊻ rb) || continue              # only edges that cross the right boundary
        q = z[idx.q_idx[e]]
        Q += rb ? q : -q                   # + if boundary endpoint is node_out, - if node_in
    end
    return Q
end

function plot_flow(sol, idx, edges_df, nodes_df)

    Qmag = [abs(sol.u[idx.q_idx[e]]) for e in 1:idx.nedges]
    qmin = 0.0; qmax = extrema(Qmag)[2]
    den = max(qmax - qmin, eps())     # avoid divide-by-zero if all equal

    grad = cgrad(:hot)                # choose any Makie colormap you like

    fig = Figure()

    ax = Axis(
        fig[1, 1]; 
        aspect = DataAspect(),
        limits = ((0, L), (0, L)),
        xgridvisible = false,
        ygridvisible = false,
        backgroundcolor = :black,
    )

    for (e, r) in enumerate(eachrow(edges_df))
        x = [nodes_df.x[r.node_in],  nodes_df.x[r.node_out]]
        y = [nodes_df.y[r.node_in],  nodes_df.y[r.node_out]]

        t = (Qmag[e] - qmin)/den      # normalise to [0,1]
        col = grad[t]                 # pick colour from colormap

        lines!(
            ax,
            x,
            y;
            color = col,
            linewidth = 8.0,
        )
    end

    scatter!(
        ax,
        data[!,1],
        data[!,2],
        marker = Circle,
        color = :grey,
        markersize = 2.0*R,
        markerspace = :data,
    )

    # 3 colourbar with matching limits
    Colorbar(fig[1, 2]; colormap=:hot, limits=(qmin, qmax), label="|Q|")

    display(fig)
    save("flow_rate.pdf", fig, pdf_version="1.4")
    save("flow_rate.png", fig)

    return nothing
end

main();