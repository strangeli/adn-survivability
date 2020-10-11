using OrdinaryDiffEq: ODEProblem, Rodas5, solve!, DiscreteCallback, CallbackSet, _initialize_dae!, u_modified!
using DiffEqCallbacks: GeneralDomain
using DiffEqBase: solve
using Setfield

import PowerDynamics: simulate

"""
```Julia
modified_copy(component, key, value)
```

provides a copy of an immutable type and modifies a given field

# Arguments
- `component`: type instance
- `key`: field name that should be modified
- `value`: new value for field `key`

Note that `value` needs to have the appropriate field type.
"""
function modified_copy(component, key::Symbol, value)
    kwargs = Dict{Symbol,Any}()
    # collect values
    for fn in component |> typeof |> fieldnames
        kwargs[fn] = getfield(component, fn)
    end
    # modify field
    @assert value isa fieldtype(typeof(component), key)
    kwargs[key] = value
    return eval(typeof(component))(; kwargs...)
end

using NLsolve: nlsolve, converged
using LinearAlgebra: pinv

struct RootRhs_ic
    rhs
    mpm
end
function (rr::RootRhs_ic)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.0)
    rr.mpm * dx .- dx
end

function RootRhs_ic(of::ODEFunction)
    mm = of.mass_matrix
    @assert mm != nothing
    mpm = pinv(mm) * mm
    RootRhs_ic(of.f, mpm)
end

function find_valid_initial_condition(pg, ic_guess)
    rr = RootRhs_ic(rhs(pg))
    nl_res = nlsolve(
        rr,
        ic_guess,
        autodiff = :forward,
        show_trace = false,
        method=:trust_region,
    )
    if converged(nl_res) == true
        return nl_res.zero #State(pg, nl_res.zero)
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end

"""
```Julia
    NodeShortCircuitAlt(;node_number,Y,tspan_fault)
```
# Keyword Arguments
- `node_number`: number  of the node
- `Y`: admittance of the short circuit
- `tspan_fault`: short circuit timespan
"""
Base.@kwdef struct NodeShortCircuitAlt
    node_number
    Y
    tspan_fault
    shunt_symbol = :Y_n
end

"Error to be thrown if something goes wrong during short circuit"
struct NodeShortCircuitAltError <: PowerDynamicsError
    msg::String
end

function (nsc::NodeShortCircuitAlt)(powergrid)
    nodes = copy(powergrid.nodes)
    sc_node = powergrid.nodes[nsc.node_number]

    # add shunt to the node component at the fault location
    if !(hasproperty(sc_node, nsc.shunt_symbol))
        throw(NodeShortCircuitAltError("Node number: $(nsc.node_number) must have a shunt field called $(nsc.shunt_symbol)."))
    end
    lens = Setfield.compose(Setfield.PropertyLens{nsc.shunt_symbol}())
    faulted_component = Setfield.set(sc_node, lens, nsc.Y)

    nodes[nsc.node_number] = faulted_component
    ControlledPowerGrid(
        powergrid.graph,
        nodes,
        powergrid.lines,
        powergrid.controller,
    )
end

"""
```Julia
simulate(nsc::NodeShortCircuitAlt, powergrid, x1, timespan)
```
Simulates a [`NodeShortCircuitAlt`](@ref)
"""
function simulate(nsc::NodeShortCircuitAlt, powergrid, x1, timespan; additional_tstops=[])
    @assert first(timespan) <= nsc.tspan_fault[1] "fault cannot begin in the past"
    @assert nsc.tspan_fault[2] <= last(timespan) "fault cannot end in the future"

    nsc_powergrid = nsc(powergrid)

    problem = ODEProblem{true}(rhs(powergrid), x1, timespan)

    init_alg =  BrownFullBasicInit()

    function errorState(integrator)
        sol1 = integrator.sol
        println("searching fault state")
        #x2 = find_valid_initial_condition(nsc_powergrid, sol1[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(nsc_powergrid)
        DiffEqBase.initialize_dae!(integrator, init_alg)
        println("found fault state")
        println(State(powergrid, integrator.u)[:, :v])
        #integrator.u = x2
        #u_modified!(integrator,true)
    end

    function regularState(integrator)
        sol2 = integrator.sol
        println("searching post-fault state")
        #x3 = find_valid_initial_condition(powergrid, sol2[end]) # Jump the state to be valid for the new system.
        integrator.f = rhs(powergrid)
        DiffEqBase.initialize_dae!(integrator, init_alg)
        println("found post-fault state")
        println(State(powergrid, integrator.u)[:, :v])
        #integrator.u = x3
        #u_modified!(integrator,true)
    end

    t1 = nsc.tspan_fault[1]
    t2 = nsc.tspan_fault[2]
    dt = 1e-7

    # discrete callback requires tstop
    cb1 = DiscreteCallback(
        ((u, t, integrator) -> t in nsc.tspan_fault[1]),
        errorState,
    )
    cb2 = DiscreteCallback(
        ((u, t, integrator) -> t in nsc.tspan_fault[2]),
        regularState,
    )

    # cb1 = ContinuousCallback(
    #     ((u, t, integrator) -> t - nsc.tspan_fault[1]),
    #     errorState,
    # )
    # cb2 = ContinuousCallback(
    #     ((u, t, integrator) -> t - nsc.tspan_fault[2]),
    #     regularState,
    # )


    # Vlim = 0.1 #pu
    # v_idx = Vector{Int}()
    # offset = 0.
    # for node in cpg.nodes
    #     push!(v_idx, offset + 1)
    #     push!(v_idx, offset + 2)
    #     offset += dimension(node)
    # end
    # function is_valid_voltage(u, resid)
    #     #v = State(powergrid, u)[:, :v]
    #     #resid .= sum(map(x-> x<Vlim ? Vlim - x : 0., v))
    #     for i in 1:length(u)
    #         if i ∈ v_idx
    #             resid[i] =  abs(u[i]) < Vlim ? abs(u[i]) - Vlim : 0.
    #         else
    #             resid[i] = 0.
    #         end
    #     end
    # end
    #
    # cb3 = GeneralDomain(is_valid_voltage; autonomous=true)

    tstops = [t1, t2]
    append!(tstops, additional_tstops)

    sol = solve(
        problem,
        Rodas5(autodiff = true),
        #force_dtmin = true,
        initializealg=BrownFullBasicInit(),
        callback = CallbackSet(cb1, cb2),
        tstops = tstops,
        #tstops = [t1-dt, t1+dt, t2-dt, t2+dt],
    )
    return PowerGridSolution(sol, powergrid)
end

simulate(nsc::NodeShortCircuitAlt, op::State, timespan) =
    simulate(nsc, op.grid, op.vec, timespan)
