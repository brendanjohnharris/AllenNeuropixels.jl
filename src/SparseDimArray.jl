using SparseArrays
using DimensionalData
export AbstractSparseDimArray, SparseDimArray
import DimensionalData: dims, refdims, data, name, metadata, parent, rebuild

abstract type AbstractSparseDimArray{T,N,D,A} <: AbstractDimArray{T,N,D,A} end

struct SparseDimArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractSparseArray{T,Int64,N},Na,Me} <: AbstractSparseDimArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
end

DimensionalData.format(dims, x::AbstractSparseDimArray) = DimensionalData.format(dims, axes(x))
# DimensionalData.format(::Base.KeySet, ::Type, ::AbstractRange) = DimensionalData.format(DimensionalData.Dimensions.LookupArrays.Aligned, ::Type, ::AbstractRange)
# format(val::AbstractSparseArray, D::Type, axis::AbstractRange) = format(AutoLookup(), D, val, axis)

function SparseDimArray(data::A, dims::D, refdims::R=(), name::Na=DimensionalData.NoName(), metadata=DimensionalData.NoMetadata()) where {D,R,A,Na}
    @info "Constructing a sparsedimarray properly"
    SparseDimArray(data, DimensionalData.format(dims, axes(data)), refdims, name, metadata)
end

function goSparseDimArray(data::A, dims::D, refdims::R=(), name::Na=DimensionalData.NoName(), metadata=DimensionalData.NoMetadata()) where {D,R,A,Na}
    @info "Constructing a sparsedimarray properly"
    SparseDimArray(data, DimensionalData.format(dims, axes(data)), refdims, name, metadata)
end



# function FeatureArray(data::AbstractArray, features::Union{Tuple{Symbol}, Vector{Symbol}}, args...)
#     FeatureArray(data, (Dim{:feature}(features), fill(AnonDim(), ndims(data)-1)...), args...)
# end

# FeatureArray(D::DimArray) = FeatureArray(D.data, D.dims, D.refdims, D.name, D.metadata)
# # DimensionalData.DimArray(D::FeatureArray) = DimArray(D.data, D.dims, D.refdims, D.name, D.metadata)

dims(A::AbstractSparseDimArray) = A.dims
refdims(A::AbstractSparseDimArray) = A.refdims
data(A::AbstractSparseDimArray) = A.data
name(A::AbstractSparseDimArray) = A.name
metadata(A::AbstractSparseDimArray) = A.metadata
parent(A::AbstractSparseDimArray) = data(A)
Base.Array(A::AbstractSparseDimArray) = parent(A)

function DimensionalData.rebuild(A::AbstractSparseDimArray, data::AbstractSparseArray, dims::Tuple, refdims::Tuple, name, metadata)
    SparseDimArray(data, dims, refdims, name, metadata)
end


AbstractSparseDimVector = AbstractSparseDimArray{T,1,D,A} where {T, D, A}
