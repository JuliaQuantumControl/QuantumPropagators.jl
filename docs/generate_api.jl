#! format: off
import QuantumPropagators
import Documenter
import TimerOutputs


"""Return a list of symbols for the names directly defined in `pkg`.

This filters out re-exported names and sub-modules. By default, for `all=true`,
both public (exported) and private names are included. With `all=false`, the
list is filtered to include only exported names.
"""
function get_local_members(pkg; all=true)
    return [
        m for m in names(pkg, all=all) if !(
            (startswith(String(m), "#")) ||       # compiler-generated names
            (getfield(pkg, m) isa Union{Dict,Array,Set}) ||  # global variable
            (m == Symbol(pkg)) ||                 # the package itself
            (m == :eval) ||                       # compiler-injected "eval"
            (m == :include) ||                    # compiler-injected "include"
            ((getfield(pkg, m)) isa Module) ||    # sub-modules
            (_parentmodule(getfield(pkg, m), pkg)) ≠ pkg  # re-exported
        )
    ]
end

_parentmodule(m, pkg) = parentmodule(m)
_parentmodule(m::TimerOutputs.TimerOutput, pkg) = pkg
_parentmodule(m::Number, pkg) = pkg
_parentmodule(m::VersionNumber, pkg) = pkg


"""Return a list of Module instances for all the sub-modules of `pkg`.
"""
function get_submodules(pkg)
    return [
        getfield(pkg, m) for m in names(pkg, all=true)
        if (getfield(pkg, m) isa Module) && !(m == Symbol(pkg))
    ]
end


"""Return the canonical fully qualified name of an object (function, type).

Given e.g. a function object, it returns a string containing the canonical
fully qualified name of the original function definition.
"""
function canonical_name(obj)
    mod = parentmodule(obj)
    modpath = fullname(mod)
    modname = join((String(sym) for sym in modpath), ".")
    objname = String(nameof(obj))
    return "$modname.$objname"
end


module_api_id(mod) = replace("$mod", "." => "") * "API"


function summarize_submodules(pkg)
    lines = String[
        "The full list of sub-modules and their public members is:",
        ""
    ]
    for mod in get_submodules(pkg)
        push!(lines, "* [`$mod`](@ref $(module_api_id(mod)))")
        for name in get_local_members(mod, all=false)
            push!(lines, "    * [`$name`](@ref $mod.$name)")
        end
    end
    return join(lines, "\n")
end


PREAMBLE = raw"""
The highest-level API of the `QuantumPropagators.jl` package (apart from some [convenience functions](@ref "Convenience functions")) consists of a single function:

* [`propagate`](@ref) — Propagate a quantum state over an entire time grid under a given generator

At a slightly lower level, propagation of quantum states in encapsulated by [The Propagator interface](@ref):

* [`init_prop`](@ref) — Initialize a `propagator` object, which is of some concrete (method-dependent) sub-type of [`AbstractPropagator`](@ref QuantumPropagators.AbstractPropagator)
* [`reinit_prop!`](@ref) — Re-initialize the `propagator`
* [`prop_step!`](@ref) — Advance the `propagator`  by a single time step forward or backward

In some cases, the ability to mutate the propagator after each propagation step can be useful. This can be achieved with the following private methods:

* [`set_state!`](@ref QuantumPropagators.set_state!) — Mutate the current quantum `state` of the `propagator` (**not exported**)
* [`set_t!`](@ref QuantumPropagators.set_t!) — Mutate the current time of the `propagator` (**not exported**)

The dynamics of a quantum state are determined by a time-dependent dynamical generator (a Hamiltonian or Liouvillian).
The `QuantumPropagators` package re-exports the two main initialization routines for generators from [`QuantumPropagators.Generators`](@ref QuantumPropagatorsGeneratorsAPI):

* [`hamiltonian`](@ref) — Construct a time-dependent generator for a propagation in Hilbert space under the Schrödinger equation
* [`liouvillian`](@ref) — Construct a time-dependent generator for a propagation in Liouville space under the master equation in Lindblad form

To set up the time-dependent control fields in a Hamiltonian, methods from the submodules [`QuantumPropagators.Controls`](@ref QuantumPropagatorsControlsAPI), [`QuantumPropagators.Shapes`](@ref QuantumPropagatorsShapesAPI), and [`QuantumPropagators.Amplitudes`](@ref QuantumPropagatorsAmplitudesAPI) can be used.

The above constitutes the main interface of `QuantumPropagators`. At the lowest level, further functionality is provided by sub-modules like [`QuantumPropagators.Cheby`](@ref QuantumPropagatorsChebyAPI), which defines a standalone API specifically for the Chebychev propagation method.
""" * summarize_submodules(QuantumPropagators) * """

### Convenience functions

There are some "convenience functions" that wrap around [`propagate`](@ref) for common tasks:

* [`propagate_sequence`](@ref)

"""


function write_module_api(out, mod, description; reference_title="Reference")

    members = [
        m for m in names(mod)
        if !(
            (String(Symbol(mod)) |> endswith(".$m")) ||
            m == Symbol(mod)
        )
    ]

    public_members = get_local_members(mod, all=false)

    all_local_members = get_local_members(mod, all=true)

    documented_members = [
        k.var for k in keys(Documenter.DocSystem.getmeta(mod))
    ]
    documented_private_members = [
        name for name in documented_members
        if (name ∉ public_members) && (name ∈ all_local_members)
    ]

    reexported_members = [
        m for m in members
        if m ∉ public_members
    ]

    write(out, "\n\n## [`$mod`](@id $(module_api_id(mod)))\n\n")
    if length(description) > 0
        write(out, "\n\n")
        write(out, description)
        write(out, "\n\n")
    end
    write(out, "\n\n### $reference_title\n\n")
    if length(public_members) > 0
        write(out, "\nPublic Members:\n\n")
        for name ∈ public_members
            println(out, "* [`$name`](@ref $mod.$name)")
        end
        write(out, "\n")
    end
    if length(reexported_members) > 0
        write(out, "\nRe-exported Members:\n\n")
        for name ∈ reexported_members
            obj = getfield(mod, name)
            ref = canonical_name(obj)
            println(out, "* [`$name`](@ref $ref)")
        end
        write(out, "\n")
    end
    if length(documented_private_members) > 0
        write(out, "\nPrivate Members:\n")
        for name ∈ documented_private_members
            println(out, "* [`$name`](@ref $mod.$name)")
        end
        write(out, "\n")
    end
    if length(public_members) > 0
        write(out, "\n\n#### Public members\n\n")
        println(out, "```@docs")
        for name ∈ public_members
            println(out, "$mod.$name")
        end
        println(out, "```")
    end
    if length(documented_private_members) > 0
        write(out, "\n\n#### Private members\n\n")
        println(out, "```@docs")
        for name ∈ documented_private_members
            println(out, "$mod.$name")
        end
        println(out, "```")
    end

end


outfile = joinpath(@__DIR__, "src", "api", "quantumpropagators.md")
println("Generating API for QuantumPropagators in $outfile")
open(outfile, "w") do out
    write(out, "```@meta\n")
    write(out, "EditURL = \"../../generate_api.jl\"\n")
    write(out, "```\n\n")
    write(out, "# API")
    write_module_api(
        out,
        QuantumPropagators,
        PREAMBLE;
        reference_title="Reference for top-level `QuantumPropagators` module"
    )
    write_module_api(
        out,
        QuantumPropagators.Generators,
        """The following routines define dynamical generators (Hamiltonians,
        Liouvillians). This includes the initialization of generators and the
        methods that define how these generators contain controls and control
        amplitudes. These methods must be defined for any custom generator or
        control amplitude.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Shapes,
        """The following routines define useful function ``S(t)`` that can be
        used for control functions or amplitudes.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Controls,
        """The following routines define methods that must be defined for any
        control function, or for generators with respect to control functions.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.ParameterizedFunctions,
        """The following types implement useful instances of
        [`QuantumPropagators.Controls.ParameterizedFunction`](@ref).
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Amplitudes,
        raw"""The following types define control amplitudes
        ``a(\{ϵ_l(t)\}, t)`` that depend on one or more control functions
        ``ϵ_l(t)``.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Storage,
        """The following routines allow to manage and extend storage arrays
        `storage` parameter in [`propagate`](@ref). See the discussion of
        [Expectation Values](@ref) for more details.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Cheby,
        """The following routines implement time propagation via expansion in
        Chebychev polynomials.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Newton,
        """The following routines implement time propagation via expansion in
        Newton polynomials using a restarted Arnoldi scheme to determine
        evaluation points.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.ExpProp,
        """The following routines implement time propagation via explicit
        exponentiation of the dynamical generator.
        """
    )
    write_module_api(
        out,
        QuantumPropagators.SpectralRange,
        """[Chebychev propagation](@ref QuantumPropagatorsChebyAPI) relies on
        estimating the spectral range of the Hamiltonian, which in turn may be
        done via [Arnoldi iteration](@ref QuantumPropagatorsArnoldiAPI).
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Arnoldi,
        """[Arnoldi iteration](https://en.wikipedia.org/wiki/Arnoldi_iteration)
        is an approximate method to find extremal eigenvalues of a dynamical
        generator by projecting it into a
        [Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace). It is
        used to estimate spectral ranges for
        [Chebychev Propagation](@ref QuantumPropagatorsChebyAPI) and to find
        appropriate Leja points to support
        [Newton Propagation](@ref QuantumPropagatorsNewtonAPI)
        """
    )
    write_module_api(
        out,
        QuantumPropagators.Interfaces,
        """The following routines allow to test whether custom data structures
        match the interface requirements of `QuantumPropagators`.
        """
    )
    write(out, """

    ## Index

    ```@index
    ```
    """)
end
