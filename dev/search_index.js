var documenterSearchIndex = {"docs":
[{"location":"benchmarks/#Benchmarks","page":"Benchmarks","title":"Benchmarks","text":"","category":"section"},{"location":"benchmarks/","page":"Benchmarks","title":"Benchmarks","text":"TODO","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#High-level-routines","page":"API","title":"High level routines","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"The QuantumPropagators.jl package provides the following high-level routines:","category":"page"},{"location":"api/","page":"API","title":"API","text":"initpropwrk — Initialize a work space for propagation\npropagate — Propagate a state over an entire time grid.\npropstep! — Perform a single propagation step in-place.","category":"page"},{"location":"api/","page":"API","title":"API","text":"These delegate the routines implementing the various propagation methods, detailed in the reference sections below.","category":"page"},{"location":"api/","page":"API","title":"API","text":"See the Index for the full list of routines.","category":"page"},{"location":"api/#Reference","page":"API","title":"Reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumPropagators]\nPages = [\"propagate.jl\"]","category":"page"},{"location":"api/#QuantumPropagators.initpropwrk-Tuple{Any, Any, Val{:auto}, Vararg{Any, N} where N}","page":"API","title":"QuantumPropagators.initpropwrk","text":"Initialize a workspace for propagation.\n\nwrk = initpropwrk(state, tlist, method=:auto, generator...; kwargs...)\n\nThe resulting wrk can be passed to propagate or propstep!.\n\nArguments\n\nstate: An exemplary state for the propagation (e.g., the initial state)\ntlist: The time grid over which propagate will be called. Must include at least to points in order to determine the propagation time step to prepare. If the propagation will be over a tlist with a variable dt, the full tlist must be passed here.\ngenerator: An exemplary (non-time-dependent) dynamical generator. For full generality (if method=:cheby), the given generator should have a spectral range sufficiently large to encompass the entire propagation. If given multiple times, a spectral envelope enclosing all the generators will be determined automatically. In this case, you should pass the generators with the extremal values of all the controls.\nmethod: The propagation method to use. The default value of :auto attempts to choose the best method available, based on the properties of the given state, tlist, and generator. Alternative values are :cheby and :newton, and :expprop.\n\nAll other kwargs are filtered and passed to the contructor for returned workspace, e.g. limit for method=:cheby or m_max for method=:newton.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.propagate","page":"API","title":"QuantumPropagators.propagate","text":"Propagate a state over an entire time grid.\n\npropagate(state, genfunc, tlist, method=:auto;\n          backwards=false; storage=nothing,\n          observables=(<store state>, ), hook=nothing)\n\npropagates state over the time grid in tlist, using piecewise-constant dynamical generators (Hamiltonians or Liouvillians) determined by genfunc, and returns the resulting propagated state. The propagation is performed by calling propstep! for every interval in tlist.\n\nFor the i'th time interval, genfunc(tlist, i) must return the generator for that time interval. Generally, when approximating a time-continuous dynamical generator as piecewise-constant on the time grid, it should be evaluated at the midpoint of the interval. A possible exception is the first and last interval, which may be better evaluated at tlist[1] and tlist[end] to ensure exact boundary conditions like control fields that are exactly zero.\n\nIn addition to the two positional parameters indicating the time interval, genfunc will also receive the state (the input state for the propagation step), backwards, storage, observables, and init as keyword arguments. These additional parameters may be used for unusual equations of motion beyond the standard Schrödinger or Liouville-von-Neumann equation, e.g. state would enter the genfunc for a Gross–Pitaevskii equation. For standard equations of motion that do not use the additional parameters, it is best to capture the keyword arguments to genfunc with a definition like\n\ngenfunc(tlist, i; kwargs...) = ...\n\nFor valid propagation methods, see initpropwrk.\n\nIn general, there is no requirement that tlist has a constant time step, although some propagation methods (most notably cheby!) only support a uniform time grid.\n\nIf storage is given as an Array, it will be filled with data determined by the observables. The default \"observable\" results in the propagated states at every point in time being stored. The storage array should be created with init_storage. See its documentation for details.\n\nThe storage parameter may also be given as true, and a new storage array will be created internally with init_storage and returned instead of the propagated state.\n\nIf backwards is true, the input state is assumed to be at time tlist[end], and the propagation progresses backwards in time (with a negative time step dt). If storage is given, it will be filled back-to-front during the backwards propagation.\n\nIf hook is given as a callable, it will be called after each propagation step, as hook(state, generator, tlist, i, wrk, observables) where i is the index of the time interval on tlist covered by the propagation step (0 for the initial state, respectives lastindex(tlist) for the backward propagation).  The hook is called before calculating any observables. Example usage includes writing data to file, or modyfing state, e.g. removing amplitude from the lowest and highest level to mitigate \"truncation error\".\n\nThe propagate routine returns the propagated state at tlist[end], respectively tlist[1] if backwards=true, or a storage array with the stored states / observable data if storage=true.\n\n\n\n\n\n","category":"function"},{"location":"api/#QuantumPropagators.propstep!-Tuple{Any, Any, Any, ChebyWrk}","page":"API","title":"QuantumPropagators.propstep!","text":"Perform a single propagation step in-place.\n\npropstep!(state, generator, dt, wrk;, kwargs...)\n\nThe propagation method is determined by wrk, see initpropwrk. The kwargs are forwarded to the underlying method\n\n\n\n\n\n","category":"method"},{"location":"api/#Storage-Reference","page":"API","title":"Storage Reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"The following routine allow manage and extend storage arrays storage parameter in propagate. See Storage of states or expectation values for more details.","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumPropagators]\nPages = [\"storage.jl\"]","category":"page"},{"location":"api/#QuantumPropagators.get_from_storage!-Tuple{Any, AbstractVector{T} where T, Any}","page":"API","title":"QuantumPropagators.get_from_storage!","text":"Obtain data from storage\n\nget_from_storage!(state, storage, i)\n\nextracts data from the storage for the i'th time slot. Invese of write_to_storage!\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.init_storage-Tuple{Any, AbstractVector{T} where T}","page":"API","title":"QuantumPropagators.init_storage","text":"Create a storage array for propagate.\n\nstorage = init_storage(state, tlist)\n\ncreates a storage array suitable for storing a state for each point in tlist.\n\nstorage = init_storage(state, tlist, observables))\n\ncreates a storage array suitable for the data generated by the observables applied to state, see map_observables, for each point in tlist.\n\nstorage = init_storage(data, nt))\n\ncreates a storage arrays suitable for storing data nt times, where nt=length(tlist). By default, this will be a vector of typeof(data) and length nt, or a n × nt Matrix with the same eltype as data if data is a Vector of length n.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.map_observable-Tuple{Any, Any}","page":"API","title":"QuantumPropagators.map_observable","text":"Apply a single observable to state.\n\ndata = map_observable(observable, state)\n\nBy default, observable is assumed to be callable, and the above is equivalent to data = observable(state).\n\nIf observable is a matrix and state is a vector evaluate the expectation value of the observable as dot(state, observable, state).\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.map_observables-Tuple{Any, Any}","page":"API","title":"QuantumPropagators.map_observables","text":"Obtain \"observable\" data from state.\n\ndata = map_observables(observables, state)\n\ncalculates the data for a tuple of observables applied to state. For a single observable (tuple of length 1), simply return the result of map_observable.\n\nFor multiple observables, return the tuple resulting from applying map_observable for each observable. If the tuple is \"uniform\" (all elements are of the same type, e.g. if each observable calculates the expectation value of a Hermitian operator), it is converted to a Vector. This allows for compact storage in a storage array, see init_storage.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.write_to_storage!-Tuple{Any, Integer, Any, Any}","page":"API","title":"QuantumPropagators.write_to_storage!","text":"Place data into storage for time slot i.\n\n    write_to_storage!(storage, i, state, observables)\n\nFor a storage array created by init_storage, store the data obtains from map_observables into the storage for time slot i. This delegates to the more general\n\n    write_to_storage!(storage, i, data)\n\nConceptually, this corresponds roughly to storage[i] = data, but storage may have its own idea on how to store data for a specific time slot. For example, with the default init_storage Vector data will be stored in a matrix, and write_to_storage! will in this case write data to the i'th column of the matrix.\n\nFor a given type of storage and data, it is the developer's responsibility that init_storage and write_to_storage! are compatible.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chebychev-reference","page":"API","title":"Chebychev reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"The following routines implement time propagation via expansion in Chebychev polynomials.","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumPropagators]\nPages = [\"cheby.jl\"]","category":"page"},{"location":"api/#QuantumPropagators.ChebyWrk","page":"API","title":"QuantumPropagators.ChebyWrk","text":"Workspace for the Chebychev propagation routine.\n\n    ChebyWrk(Ψ, Δ, E_min, dt; limit=1e-12)\n\ninitializes the workspace for the propagation of a state similar to Ψ under a Hamiltonian with eigenvalues between E_min and E_min + Δ, and a time step dt. Chebychev coefficients smaller than the given limit are discarded.\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumPropagators.cheby!-NTuple{4, Any}","page":"API","title":"QuantumPropagators.cheby!","text":"Evaluate Ψ = exp(-i H dt) Ψ in-place.\n\nArgs:\n\nΨ: on input, initial vector. Will be overwritten with result.\nH: Hermitian operator\ndt: time step\nE_min: minimum eigenvalue of H, to be used instead of the E_min from the  initialization of wrk. The same wrk may be used for different values  E_min, as long as the spectra radius Δ and the time step dt are the  same as those used for the initialization of wrk.\n\nThe routine will not allocate any internal storage. This implementation requires copyto! lmul!, and axpy! to be implemented for Ψ, and the three-argument mul! for Ψ and H.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.cheby_coeffs!","page":"API","title":"QuantumPropagators.cheby_coeffs!","text":"Calculate Chebychev coefficients in-place.\n\nn = cheby_coeffs!(coeffs, Δ, dt, limit=1e-12)\n\noverwrites the first n values in coeffs with new coefficients larger than limit for the given new spectral radius Δ and time step dt. The coeffs array will be resized if necessary, and may length > n on exit.\n\n\n\n\n\n","category":"function"},{"location":"api/#QuantumPropagators.cheby_coeffs-Tuple{Any, Any}","page":"API","title":"QuantumPropagators.cheby_coeffs","text":"Calculate Chebychev coefficients.\n\nReturn an array of coefficiencts larger than limit.\n\nArguments\n\nΔ: the spectral radius of the underlying operator\ndt: the time step\n\n\n\n\n\n","category":"method"},{"location":"api/#Newton-reference","page":"API","title":"Newton reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"The following routines implement time propagation via expansion in Newton polynomials using a restarted Arnoldi scheme to determine evaluation points.","category":"page"},{"location":"api/#Public-members","page":"API","title":"Public members","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumPropagators]\nPages = [\"newton.jl\"]\nPublic = true\nPrivate = false","category":"page"},{"location":"api/#QuantumPropagators.NewtonWrk","page":"API","title":"QuantumPropagators.NewtonWrk","text":"    NewtonWrk(v0, m_max=10)\n\nWorkspace for the Newton-with-restarted-Arnoldi propagation routine.\n\nInitializes the workspace for the propagation of a vector v0, using a maximum Krylov dimension of m_max in each restart iteration. Note that m_max should be smaller than the length of v0.\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumPropagators.newton!","page":"API","title":"QuantumPropagators.newton!","text":"newton!(Ψ, H, dt, wrk, func=(z -> exp(-1im*z)); norm_min=1e-14, relerr=1e-12,\n        max_restarts=50)\n\nEvaluate Ψ = func(H*dt) Ψ using a Newton-with-restarted-Arnoldi scheme.\n\nArguments\n\nΨ: The state to propagate, will be overwritten in-place with the propagated state\nH: Operator acting on Ψ. Together with dt, this is the argument to func\ndt: Implicit time step. Together with H, this is the argument to func\nwkr: Work array, initialized with NewtonWrk\nfunc: The function to apply to H dt, taking a single (scalar) complex-valued argument z in place of H dt. The default func is to evaluate the time evoluation operator for the Schrödinger equation\nnorm_min: the minium norm at which to consider a state similar to Ψ as zero\nrelerr: The relative error defining the convergence condition for the restart iteration. Propagation stops when the norm of the accumulated Ψ is stable up to the given relative error\nmax_restart: The maximum number of restart iterations. Exceeding max_restart will throw an AssertionError.\n\n\n\n\n\n","category":"function"},{"location":"api/#Private-members","page":"API","title":"Private members","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumPropagators]\nPages = [\"newton.jl\"]\nPublic = false\nPrivate = true","category":"page"},{"location":"api/#QuantumPropagators.arnoldi!-Union{Tuple{T}, Tuple{Matrix{ComplexF64}, Array{T, N} where N, Int64, T, Any}, Tuple{Matrix{ComplexF64}, Array{T, N} where N, Int64, T, Any, Float64}} where T","page":"API","title":"QuantumPropagators.arnoldi!","text":"arnoldi!(Hess, q, m, Ψ, H, dt=1.0; extended=true, norm_min=1e-15)\n\nCalculate the Hessenberg matrix and Arnoldi vectors of H dt, from Ψ.\n\nFor a given order m, the m×m Hessemberg matrix is calculated and stored in in the pre-allocated Hess. Further  an array of m normalized Arnoldi vectors is stored in in the pre-allocated q, plus one additional unnormalized Arnoldi vector.  The unnormalized m+1st vector could be used to easily extend a given m×m Hessenberg matrix to a (m+1)×(m+1) matrix.\n\nIf the extended Hessenberg matrix is requested (extended=true, default), the m+1st Arnoldi vector is also normalized, and it's norm will be stored in m+1, m entry of the (extended) Hessenberg matrix, which is an (m+1)×(m+1) matrix.\n\nReturn the size m of the calculated Hessenberg matrix. This will usually be the input m, except when the Krylov dimension of H starting from Ψ is less then m. E.g., if Ψ is an eigenstate of H, the returned m will be 1.\n\nSee http://en.wikipedia.org/wiki/Arnoldi_iteration for a description of the algorithm.\n\nArguments\n\nHess::Matrix{ComplexF64}: Pre-allocated storage for the Hessemberg matrix.  Can be uninitialized on input. The matrix must be at least of size m×m, or  (m+1)×(m+1) if extended=true. On output, the m×m sub-matrix of Hess  (with the returned output m) will contain the Hessenberg matrix, and all  other elements of Hess be be set to zero.\nq: Pre-allocated array of states similar to Ψ, as storage for the calculated Arnoldi vectors. These may be un-initialized on input. Must be at least of length m+1\nm: The requested dimensions of the output Hessenberg matrix.\nΨ: The starting vector for the Arnoldi procedure. This can be of any type,  as long as Φ = H * Ψ results in a vector similar to Ψ, there is an inner  products of Φ and Ψ (Ψ⋅Φ is defined), and norm(Ψ) is defined.\nH: The operator (up to dt) for which to calculate the Arnoldi procedure. Can be of any type, as long as H * Ψ is defined.\ndt: The implicit time step; the total operator for which to calculate the Arnoldi procedure is H * dt\nextended: If true (default), calculate the extended Hessenberg matrix, and normalized the final Arnoldi vector\nnorm_min: the minimum value of the norm of Ψ at which Ψ should be  considered the zero vector\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.diagonalize_hessenberg_matrix-Tuple{Any, Any}","page":"API","title":"QuantumPropagators.diagonalize_hessenberg_matrix","text":"diagonalize_hessenberg_matrix(Hess, m; accumulate=false)\n\nDiagonalize the m × m top left submatrix of the given Hessenberg matrix.\n\nIf accumulate is true, return the concatenated eigenvalues for Hess[1:1,1:1] to Hess[1:m,1:m], that is, all sumatrices of size 1 through m.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.extend_leja!-Tuple{OffsetArrays.OffsetVector{ComplexF64, AA} where AA<:AbstractVector{ComplexF64}, Any, OffsetArrays.OffsetVector{ComplexF64, AA} where AA<:AbstractVector{ComplexF64}, Any}","page":"API","title":"QuantumPropagators.extend_leja!","text":"extend_leja!(leja, n, newpoints, n_use)\n\nGiven an array of n (ordered) Leja points, extract n_use points from newpoints, and append them to the existing Leja points. The array leja should be sufficiently large to hold the new Leja points, which are appended after index n_old. It will be re-allocated if necessary and may have a size of up to 2*(n+n_use).\n\nArguments\n\nleja: Array of leja values. Must contain the \"old\" leja values to be kept  in leja(0:n-1). On output, n_use new leja points will be in  leja(n+:n+n_use-1), for the original value of n.  The leja array must  use zero-based indexing.\nn: On input, number of \"old\" leja points in leja. On output, total number of leja points (i.e. n=n+n_use)\nnewpoints: On input, candidate points for new leja points.  The n_use best values will be chosen and added to leja. On output, the values of new_points are undefined.\nn_use: Number of points that should be added to leja\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumPropagators.extend_newton_coeffs!-Tuple{OffsetArrays.OffsetVector{ComplexF64, AA} where AA<:AbstractVector{ComplexF64}, Int64, OffsetArrays.OffsetVector{ComplexF64, AA} where AA<:AbstractVector{ComplexF64}, Any, Int64, Float64}","page":"API","title":"QuantumPropagators.extend_newton_coeffs!","text":"extend_newton_coeffs!(a, n_a, leja, func, n_leja, radius)\n\nExtend the array a of existing Newton coefficients for the expansion of the func from n_a coefficients to n_leja coefficients. Return a new value n_a=n_a+n_leja with the total number of Newton coefficients in the updated a.\n\nArguments\n\na: On input, a zero-based array of length n_a or greater, containing Newton coefficients. On output, array containing a total n_leja coefficients. The array a will be resized if necessary, and may have a length greater than n_leja on output\nn_a:  The number of Newton coefficients in a, on input. Elements of a  beyond the first n_a elements will be overwritten.\nleja: Array of normalized Leja points, containing at least n_leja elements.\nfunc: Function for which to calcluate Newton coeffiecients\nn_leja: The number of elements in leja to use for calculating new coefficients, and the total number of Newton coefficients on output\nradius: Normalization radius for divided differences\n\n\n\n\n\n","category":"method"},{"location":"api/#Matrix-exponentiation-reference","page":"API","title":"Matrix exponentiation reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"The following routines implement time propagation via explicit exponentiation of the dynamical generator.","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumPropagators]\nPages = [\"expprop.jl\"]\nPublic = true\nPrivate = false","category":"page"},{"location":"api/#QuantumPropagators.ExpPropWrk","page":"API","title":"QuantumPropagators.ExpPropWrk","text":"    ExpPropWrk(v0)\n\nWorkspace for propagation via direct matrix exponentiation.\n\nInitializes the workspace for the propagation of a vector v0\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumPropagators.expprop!","page":"API","title":"QuantumPropagators.expprop!","text":"expprop!(Ψ, H, dt, wrk, func=(z -> exp(-1im*z)))\n\nEvaluate Ψ = func(H*dt) Ψ by directly evaluating U = func(H*dt), i.e. by matrix exponentiation for the default func, and then multiplying U and Ψ in-place with mul!.\n\nThe workspace wrk must be initialized with ExpPropWrk to provide storage for a temporary state.\n\n\n\n\n\n","category":"function"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"howto/#Howtos","page":"Howtos","title":"Howtos","text":"","category":"section"},{"location":"howto/#Howto-extend-QuantumPropagators-with-a-new-propagation-method","page":"Howtos","title":"Howto extend QuantumPropagators with a new propagation method","text":"","category":"section"},{"location":"howto/","page":"Howtos","title":"Howtos","text":"Define a new workspace type that is unique to the propagation method, e.g. MyNewMethodWrk\nSpecialize the method initpropwrk for method::Val{:mynewmethod}\nSpecialize the method propstep! for wrk::MyNewMethodWrk","category":"page"},{"location":"howto/","page":"Howtos","title":"Howtos","text":"By defining only the above two methods, it becomes possible to call","category":"page"},{"location":"howto/","page":"Howtos","title":"Howtos","text":"propagate(state, genfunc, tlist, :mynewmethod; kwargs...)","category":"page"},{"location":"overview/#Overview","page":"Overview","title":"Overview","text":"","category":"section"},{"location":"overview/#Storage-of-states-or-expectation-values","page":"Overview","title":"Storage of states or expectation values","text":"","category":"section"},{"location":"overview/","page":"Overview","title":"Overview","text":"The propagate routine allows the storage of data for every point of the time grid.  This is done by passing it a storage object created with init_storage, or simply storage=true in order to create the appropriate storage automatically.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"By default, the storage will be used to store the propagated states at each point in time. More generally, what goes into storage can be customized via the observables parameter of propagate.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"After each propagation step, with a propagate state at time slot i,","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"map_observables generates data from the propagated state \nwrite_to_storage! places that data into storage for time slot i","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"After propagate returns, the get_from_storage! routine can be used to extract data from any time slot. This interface hides the internal memory organization of storage, which is set up by init_storage based on the type of state and the given observables. This system can can extended with multiple dispatch, allowing to optimize the storage for custom data types. Obviously, init_storage, map_observables, write_to_storage!, and get_from_storage! must all be consistent.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"The default implementation of these routine uses either a standard Vector or a Matrix as storage.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"Roughly speaking, when storing states, if the state of some arbitrary type, the storage will be a Vector where the i'th entry points to a copy of the propagated state at the i'th time slot. If the state is a Vector, the storage will be a Matrix containing the state for the i'th time slot in the i'th column.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"When a tuple of observables is passed to propagate, if map_observable returns data of the same type for each observable, storage will be a Matrix containing the values from the different observables for the i'th time slot in the i'th column. This would be typical for the storage of expectation values, e.g. with","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"observables=(state->dot(state, Ô₁, state), state->dot(state, Ô₂, state))","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"where Ô₁, Ô₂ are two Hermitian operators. In this case, storage would be a 2 × nt Float64 array. Calling get_from_storage!(data, storage, i) would be equivalent to copyto!(data, storage[:,i]) and extract the i'th column of storage, i.e. the vector [⟨Ô₁⟩, ⟨Ô₂A]⟩] at time slot i. Alternatively, storage[1,:] would return the values of ⟨Ô₁⟩ over time. This would be useful for plotting, and illustrates the benefits of using a Matrix as storage.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"Usually, the observables should be functions acting on the state, but map_observable can be extended to other types of observables as well. For example, the situation were state is a vector and the observables are matrices is also supported; if  Ô₁, Ô₂ are matrices, ","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"observables=(Ô₁, Ô)","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"would have the same result of storing the expectation values of those two operators.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"The observables are not required to yield real values: the term \"observables\" is used very loosely here. We could directly calculate, e.g., the complex amplitude α of a coherent state in quantum optics, the number of levels with non-zero population (as an integer), or the propagated state transformed from a moving frame to a lab frame.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"If there are multiple observables that return data of different types, by default storage will be a Vector that contains tuples with the result for each observable. For example, with","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"observables=(state->dot(state, Ô₁, state), state->count_poplevels(state))","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"where count_poplevels is a function that counts the number of levels with non-zero population, the resulting storage would be a Vector{Tuple{Float64, Int64}}.","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"If there is a single observable that yields a vector, that vector is stored in the i'th column of a storage matrix. This is in fact what happens when storing the propagated states (observables=(Ψ->copy(Ψ), )) if Ψ is a Vector, but there are other use cases, such as calculating the population in all levels in one go, with observables=(Ψ -> abs.(Ψ).^2, ).","category":"page"},{"location":"overview/","page":"Overview","title":"Overview","text":"If there is a single variable that yields a non-vector object, storage will be a Vector where the i'th entry points to the object. This is in fact what happens by default when storing states  that are e.g. instances of QuantumOptics.Ket. In such a case, it might be advisable to add new methods for init_storage and write_to_storage! that implement a more efficient in-place storage.","category":"page"},{"location":"background/#Background","page":"Background","title":"Background","text":"","category":"section"},{"location":"background/#Chebychev-Propagation","page":"Background","title":"Chebychev Propagation","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"TODO","category":"page"},{"location":"background/#Newton-Propagation","page":"Background","title":"Newton Propagation","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"TODO","category":"page"},{"location":"background/#Propagation-with-Explicit-Matrix-Exponentiation","page":"Background","title":"Propagation with Explicit Matrix Exponentiation","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"TODO","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"TODO","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Pages=[\n    \"1.md\",\n    \"2.md\",\n]\nDepth = 1","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = QuantumPropagators","category":"page"},{"location":"#QuantumPropagators","page":"Home","title":"QuantumPropagators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for QuantumPropagators.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages=[\n \"overview.md\",\n \"background.md\",\n \"examples/index.md\",\n \"howto.md\",\n \"benchmarks.md\",\n \"api.md\",\n]\nDepth = 2","category":"page"}]
}
