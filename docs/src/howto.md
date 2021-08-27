# Howtos

## Howto extend QuantumPropagators with a new propagation method

* Define a new workspace type that is unique to the propagation method, e.g. `MyNewMethodWrk`
* Specialize the method [`initpropwrk`](@ref) for `method::Val{:mynewmethod}`
* Specialize the method [`propstep!`](@ref) for `wrk::MyNewMethodWrk`

By defining only the above two methods, it becomes possible to call

~~~julia
propagate(state, genfunc, tlist, :mynewmethod; kwargs...)
~~~
