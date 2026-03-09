using Roots

@inline function flow1D(G, h, τ0, K, n, α, τS, β, t; δ=1e-8, ϵ=1e-8)
    
    absG = smooth_abs(G, δ)          # ≈ |G|
    sgnG = G / absG                  # ≈ sign(G), smooth at 0

    μ = 1000.0

    Δ0 = smooth_pospart(absG*h - τ0, ϵ)
    ΔS = smooth_pospart(absG*h - τS, ϵ)

    out = 2.0 / (absG^2) * n * K * (Δ0 / K)^((n + 1) / n)
    out *= ((n + 1.0) * absG * h + n * τ0)
    out /= ((n + 1.0) * (2.0 * n + 1.0))
    out += 2.0 * h * α * ΔS^β * t
    out += 2.0 / 3.0 * h^3/μ * absG

    return out * sgnG
end

function HBgrad(x, q, hmin, τ0, K, n, α, τS, β, R;
                δ=1e-6, ϵ=1e-6, htol=1e-12, atol=1e-10, rtol=1e-8)

    q == 0 && return 0.0
    s = sign(q)
    qabs = abs(q)

    # You integrate x in [0,R], so clamp there defensively
    x = clamp(x, 0.0, R)

    # Geometry with guards
    rad = max(R^2 - (x - R)^2, 0.0)
    hmin_eff = max(hmin, htol)             # prevent h going <= 0 at x=R
    h = R + hmin_eff - sqrt(rad)

    denom2 = max(x*(2R - x), 0.0)
    hprime = (2.0*(x - R)) / sqrt(denom2 + eps(Float64))
    t = 1.0 / sqrt(hprime^2 + 1.0)

    # Solve for G >= 0
    f(G) = flow1D(G, h, τ0, K, n, α, τS, β, t; δ=δ, ϵ=ϵ) - qabs

    a = 0.0
    fa = f(a)          # = -qabs < 0

    # Start upper bracket from yield/slip/Newtonian estimates
    G_yield = (τ0 > 0) ? τ0 / h : 0.0
    G_slip  = (τS > 0) ? τS / h : 0.0
    G_newt  = 3.0*qabs*K / (2.0*h^3 + eps())  # safe

    b = max(1.01*max(G_yield, G_slip, G_newt), 1e-12)

    fb = f(b)
    k = 0
    while fb <= 0 && k < 80
        b *= 2
        fb = f(b)
        k += 1
    end
    fb <= 0 && error("HBgrad: failed to bracket root (q=$qabs, h=$h, b=$b)")

    Gpos = Roots.find_zero(f, (a, b), Roots.Brent(); atol=atol, rtol=rtol)
    return s * Gpos
end

function Δp(q, hmin, τ0, K, n, α, τS, β, R; N::Int=51)
    sign_q = sign(q)
    q = abs(q)

    @assert N ≥ 3 "Simpson needs at least 3 points"
    @assert isodd(N) "Simpson 1/3 needs odd N (even number of subintervals)"

    h = R/(N-1)

    # f0 + fN
    s = HBgrad(0.0, q, hmin, τ0, K, n, α, τS, β, R) +
        HBgrad(R,   q, hmin, τ0, K, n, α, τS, β, R)

    # 4*odd + 2*even interior points
    @inbounds for j in 2:(N-1)
        x = (j-1)*h
        fj = HBgrad(x, q, hmin, τ0, K, n, α, τS, β, R)
        s += (isodd(j) ? 4.0 : 2.0) * fj
    end

    integral_0R = (h/3.0) * s
    return 2.0 * integral_0R * sign_q
end

@inline smooth_abs(x, δ) = sqrt(x^2 + δ^2)  # smooth approximation to |x|, smooth at 0

@inline smooth_pospart(x, ϵ) = 0.5 * (x + sqrt(x^2 + ϵ^2))  # smooth max(x,0)

