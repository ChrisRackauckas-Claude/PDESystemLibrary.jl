"""
# The Inviscid Burgers Equation in 1D

The Inviscid Burgers equation is a model for the evolution of a fluid.
The fluid is assumed to be incompressible and inviscid, meaning that the fluid is not viscous and does not change in volume.
The fluid is also assumed to be one-dimensional, meaning that the fluid is only moving in one axis.
"""
function inviscid_burgers_monotonic()
    @parameters x t
    @variables u(..)
    Dx = Differential(x)
    Dt = Differential(t)
    x_min = 0.0
    x_max = 1.0
    t_min = 0.0
    t_max = 6.0

    analytic = u(t, x) ~ x / (t + 1)

    analytic_u(t, x) = x / (t + 1)

    eq = Dt(u(t, x)) ~ -u(t, x) * Dx(u(t, x))

    bcs = [
        u(0, x) ~ x,
        u(t, x_min) ~ analytic_u(t, x_min),
        u(t, x_max) ~ analytic_u(t, x_max),
    ]

    domains = [
        t ∈ Interval(t_min, t_max),
        x ∈ Interval(x_min, x_max),
    ]

    dx = 0.05

    tags = ["1D", "Monotonic", "Inviscid", "Burgers", "Advection", "Dirichlet"]

    @named inviscid_burgers_monotonic = PDESystem(
        eq, bcs, domains, [t, x], [u(t, x)];
        analytic = analytic, metadata = tags, eval_module = @__MODULE__
    )

    return inviscid_burgers_monotonic
end

"""
# The Burgers Equation in 2D

The Burgers equation is a model for the evolution of a fluid.
This time the model has a viscosity term, which means that the fluid is viscous. The fluid is also assumed to be two-dimensional.
"""
function burgers_2d()
    @parameters x y t
    @variables u(..) v(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2

    R = 80.0

    x_min = y_min = t_min = 0.0
    x_max = y_max = 1.0
    t_max = 1.0

    #Exact solutions from: https://www.sciencedirect.com/science/article/pii/S0898122110003883

    u_exact(x, y, t) = 3 / 4 - 1 / (4 * (1 + exp(R * (-t - 4x + 4y) / 32)))
    v_exact(x, y, t) = 3 / 4 + 1 / (4 * (1 + exp(R * (-t - 4x + 4y) / 32)))
    analytic = [
        u(t, x, y) ~ u_exact(x, y, t),
        v(t, x, y) ~ v_exact(x, y, t),
    ]

    eq = [
        Dt(u(t, x, y)) + u(t, x, y) * Dx(u(t, x, y)) + v(t, x, y) * Dy(u(t, x, y)) ~
            (1 / R) * (Dxx(u(t, x, y)) + Dyy(u(t, x, y))),
        Dt(v(t, x, y)) + u(t, x, y) * Dx(v(t, x, y)) + v(t, x, y) * Dy(v(t, x, y)) ~
            (1 / R) * (Dxx(v(t, x, y)) + Dyy(v(t, x, y))),
    ]

    domains = [
        t ∈ Interval(t_min, t_max),
        x ∈ Interval(x_min, x_max),
        y ∈ Interval(y_min, y_max),
    ]

    bcs = [
        u(0, x, y) ~ u_exact(x, y, 0),
        u(t, 0, y) ~ u_exact(0, y, t),
        u(t, x, 0) ~ u_exact(x, 0, t),
        u(t, 1, y) ~ u_exact(1, y, t),
        u(t, x, 1) ~ u_exact(x, 1, t),
        v(0, x, y) ~ v_exact(x, y, 0),
        v(t, 0, y) ~ v_exact(0, y, t),
        v(t, x, 0) ~ v_exact(x, 0, t),
        v(t, 1, y) ~ v_exact(1, y, t),
        v(t, x, 1) ~ v_exact(x, 1, t),
    ]

    tags = ["2D", "Non-Monotonic", "Viscous", "Burgers", "Advection", "Dirichlet"]

    @named burgers_2d = PDESystem(
        eq, bcs, domains, [t, x, y], [u(t, x, y), v(t, x, y)],
        analytic = analytic, metadata = tags, eval_module = @__MODULE__
    )

    return burgers_2d
end

all_systems = vcat(all_systems, [inviscid_burgers_monotonic(), burgers_2d()])
