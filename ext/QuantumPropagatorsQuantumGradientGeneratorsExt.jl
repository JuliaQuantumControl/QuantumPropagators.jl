module QuantumPropagatorsQuantumGradientGeneratorsExt

using QuantumGradientGenerators: GradgenOperator, GradVector

Base.size(G::GradgenOperator, dim::Integer) = size(G.G, dim)
Base.similar(::GradVector, ::Type{T}, dims::Tuple{Int,Int}) where {T} =
    Matrix{T}(undef, dims...)
Base.similar(::GradVector, ::Type{T}, dims::Tuple{Int}) where {T} =
    Vector{T}(undef, dims[1])

end
