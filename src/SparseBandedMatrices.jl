module SparseBandedMatrices

import LinearAlgebra: mul!
using .Threads: @threads
using PrecompileTools: @setup_workload, @compile_workload

"""
    SparseBandedMatrix{T} <: AbstractMatrix{T}
    SparseBandedMatrix{T}(undef, N, M)
    SparseBandedMatrix{T}(indices, diagonals, N, M)

Construct an `N × M` matrix that stores only its nonzero diagonals. The type implements
the `AbstractMatrix` indexing, conversion, and multiplication interfaces while avoiding
storage for diagonals that are identically zero.

The `undef` form starts with no stored diagonals; assigning an element or calling
[`setdiagonal!`](@ref) creates the corresponding diagonal. The second form initializes the
stored diagonals directly.

# Arguments

- `T`: The matrix element type.
- `undef`: Selects construction without any stored diagonals.
- `indices::Vector{Int}`: Internal diagonal identifiers. For an `N × M` matrix, the
  identifier for element `(i, j)` is `N - i + j`.
- `diagonals::Vector{Vector{T}}`: Values for each identifier in `indices`. Each vector must
  have the length of its represented diagonal.
- `N::Integer`: The number of rows.
- `M::Integer`: The number of columns.

# Returns

- `SparseBandedMatrix{T}`: A sparse-diagonal matrix with dimensions `(N, M)`.

# Fields

- `size::Tuple{Int, Int}`: The logical matrix dimensions.
- `indices::Vector{Int}`: Sorted identifiers for the stored diagonals. For an
  `N x M` matrix, the identifier of entry `(i, j)` is `N - i + j`.
- `diags::Vector{Vector{T}}`: Values for the corresponding stored diagonals. The
  `k`th vector in `diags` belongs to the diagonal identified by `indices[k]`.

# Interface

`SparseBandedMatrix` follows the standard `AbstractMatrix` interface. Generic
matrix code may use `size`, `axes`, `getindex`, `setindex!`, iteration, and
conversion to `Matrix`; it should not depend on the storage fields above. An
entry on a diagonal that has not been stored reads as `zero(T)`. Assigning an
entry creates its diagonal when necessary, and assigning a complete diagonal is
more efficient with [`setdiagonal!`](@ref).

# Examples

```julia
julia> using SparseBandedMatrices

julia> matrix = SparseBandedMatrix{Float64}(undef, 4, 4);

julia> matrix[1, 1] = 5.0;

julia> matrix[3, 2] = -1.0;

julia> Matrix(matrix)
4×4 Matrix{Float64}:
 5.0   0.0  0.0  0.0
 0.0   0.0  0.0  0.0
 0.0  -1.0  0.0  0.0
 0.0   0.0  0.0  0.0
```

# Storage

Use the standard matrix interface and [`setdiagonal!`](@ref). The concrete storage fields
and diagonal identifier encoding are implementation details.
"""
struct SparseBandedMatrix{T} <: AbstractMatrix{T}
    size::Tuple{Int, Int}
    indices::Vector{Int}
    diags::Vector{Vector{T}}
    function SparseBandedMatrix{T}(::UndefInitializer, N, M) where {T}
        size = (N, M)
        indices = Int[]
        diags = Vector{T}[]
        return new(size, indices, diags)
    end
    function SparseBandedMatrix{T}(ind_vals, diag_vals, N, M) where {T}
        size = (N, M)
        perm = sortperm(ind_vals)
        indices = ind_vals[perm]
        for i in 1:(length(indices) - 1)
            @assert indices[i] != indices[i + 1]
        end
        diags = diag_vals[perm]
        return new(size, indices, diags)
    end
end

function Base.size(M::SparseBandedMatrix)
    return M.size
end

function Base.getindex(M::SparseBandedMatrix{T}, i::Int, j::Int, I::Int...) where {T}
    @boundscheck checkbounds(M, i, j, I...)
    rows, cols = size(M)
    wanted_ind = rows - i + j
    ind = searchsortedfirst(M.indices, wanted_ind)
    if (ind <= length(M.indices) && M.indices[ind] == wanted_ind)
        if (i > j)
            return M.diags[ind][j]
        else
            return M.diags[ind][i]
        end
    end
    return zero(T)
end

function Base.setindex!(M::SparseBandedMatrix{T}, val, i::Int, j::Int, I::Int...) where {T}
    @boundscheck checkbounds(M, i, j, I...)
    rows = size(M, 1)
    wanted_ind = rows - i + j
    ind = searchsortedfirst(M.indices, wanted_ind)
    if (ind > length(M.indices) || M.indices[ind] != wanted_ind)
        insert!(M.indices, ind, wanted_ind)
        insert!(M.diags, ind, zeros(T, rows - abs(wanted_ind - rows)))
    end
    if (i > j)
        M.diags[ind][j] = val isa T ? val : convert(T, val)::T
    else
        M.diags[ind][i] = val isa T ? val : convert(T, val)::T
    end
    return val
end

"""
    setdiagonal!(matrix::SparseBandedMatrix, values, lower::Bool)

Set one complete stored diagonal of `matrix`. The length of `values` identifies the
diagonal: `lower = true` places it against the bottom-left edge, while `lower = false`
places it against the top-right edge. Existing values on that diagonal are overwritten.

This operation is more efficient than assigning every element separately because it
inserts or updates the diagonal storage in one call.

# Arguments

- `matrix::SparseBandedMatrix{T}`: The matrix to modify.
- `values`: Values for the diagonal. They are converted to `T` when necessary, and their
  length must not exceed the number of rows.
- `lower::Bool`: Select the lower (`true`) or upper (`false`) diagonal of the given length.

# Returns

- `values`: The original input collection.

# Throws

- `ErrorException`: If `length(values)` exceeds the number of rows.

# Examples

```julia
julia> using SparseBandedMatrices

julia> matrix = SparseBandedMatrix{Float64}(undef, 5, 5);

julia> setdiagonal!(matrix, [3.0, 4.0, 5.0], true);

julia> setdiagonal!(matrix, [1.0, 2.0], false);

julia> findall(!iszero, Matrix(matrix))
5-element Vector{CartesianIndex{2}}:
 CartesianIndex(3, 1)
 CartesianIndex(4, 2)
 CartesianIndex(5, 3)
 CartesianIndex(1, 4)
 CartesianIndex(2, 5)
```
"""
function setdiagonal!(M::SparseBandedMatrix{T}, diagvals, lower::Bool) where {T}
    rows, cols = size(M)
    if length(diagvals) > rows
        error("size of diagonal is too big for the matrix")
    end
    if lower
        wanted_ind = length(diagvals)
    else
        wanted_ind = 2 * rows - length(diagvals)
    end

    ind = searchsortedfirst(M.indices, wanted_ind)
    if (ind > length(M.indices) || M.indices[ind] != wanted_ind)
        insert!(M.indices, ind, wanted_ind)
        insert!(M.diags, ind, diagvals isa Vector{T} ? diagvals : convert(Vector{T}, diagvals)::Vector{T})
    else
        for i in eachindex(diagvals)
            M.diags[ind][i] = diagvals[i] isa T ? diagvals[i] : convert(T, diagvals[i])::T
        end
    end
    return diagvals
end

# C = Cb + aAB
function mul!(C::Matrix{T}, A::SparseBandedMatrix{T}, B::Matrix{T}, a::Number, b::Number) where {T}
    @assert size(A, 2) == size(B, 1)
    @assert size(A, 1) == size(C, 1)
    @assert size(B, 2) == size(C, 2)
    if iszero(b)
        fill!(C, zero(T))
    else
        C .*= b
    end

    rows, cols = size(A)
    @inbounds for (ind, location) in enumerate(A.indices)
        @threads for i in 1:length(A.diags[ind])
            # value: diag[i]
            # index in array:
            #       if ind < rows(A), then index = (rows - loc + i, i)
            #       else index = (i, loc - cols + i)
            val = A.diags[ind][i] * a
            if location < rows
                index_i = rows - location + i
                index_j = i
            else
                index_i = i
                index_j = location - cols + i
            end
            #A[index_i, index_j] * B[index_j, j] = C[index_i, j]
            for j in 1:size(B, 2)
                C[index_i, j] = fma(val, B[index_j, j], C[index_i, j])
            end
        end
    end
    return C
end

# C = C*b + a*B*A
function mul!(C::Matrix{T}, A::Matrix{T}, B::SparseBandedMatrix{T}, a::Number, b::Number) where {T}
    @assert size(A, 2) == size(B, 1)
    @assert size(A, 1) == size(C, 1)
    @assert size(B, 2) == size(C, 2)

    if iszero(b)
        fill!(C, zero(T))
    else
        C .*= b
    end

    rows, cols = size(B)
    @inbounds for (ind, location) in enumerate(B.indices)
        @threads for i in eachindex(B.diags[ind])
            val = B.diags[ind][i] * a
            if location < rows
                index_i = rows - location + i
                index_j = i
            else
                index_i = i
                index_j = location - cols + i
            end
            @simd for j in 1:size(A, 1)
                C[j, index_j] = fma(val, A[j, index_i], C[j, index_j])
            end
        end
    end
    return C
end

function mul!(C::SparseBandedMatrix{T}, A::SparseBandedMatrix{T}, B::SparseBandedMatrix{T}, a::Number, b::Number) where {T}
    @assert size(A, 2) == size(B, 1)
    @assert size(A, 1) == size(C, 1)
    @assert size(B, 2) == size(C, 2)

    if iszero(b)
        empty!(C.indices)
        empty!(C.diags)
    else
        C .*= b
    end

    rows_a, cols_a = size(A)
    rows_b, cols_b = size(B)
    @inbounds for (ind_a, location_a) in enumerate(A.indices)
        @threads for i in eachindex(A.diags[ind_a])
            val_a = A.diags[ind_a][i] * a
            if location_a < rows_a
                index_ia = rows_a - location_a + i
                index_ja = i
            else
                index_ia = i
                index_ja = location_a - cols_a + i
            end
            min_loc = rows_b - index_ja + 1
            max_loc = 2 * rows_b - index_ja
            for (ind_b, location_b) in enumerate(B.indices)
                #index_ib = index_ja
                #       if ind < rows(A), then index = (rows - loc + i, i)
                #rows - loc + j = index_ja, j = index_ja - rows + loc
                #       else index = (i, loc - cols + i)
                # if location < rows(B), then
                if location_b <= rows_b && location_b >= min_loc
                    j = index_ja - rows_b + location_b
                    index_jb = j
                    val_b = B.diags[ind_b][j]
                    C[index_ia, index_jb] = muladd(val_a, val_b, C[index_ia, index_jb])
                elseif location_b > rows_b && location_b <= max_loc
                    j = index_ja
                    index_jb = location_b - cols_b + j
                    val_b = B.diags[ind_b][j]
                    C[index_ia, index_jb] = muladd(val_a, val_b, C[index_ia, index_jb])
                end
            end
        end
    end
    return C
end

function mul!(C::Matrix{T}, A::SparseBandedMatrix{T}, B::SparseBandedMatrix{T}, a::Number, b::Number) where {T}
    @assert size(A, 2) == size(B, 1)
    @assert size(A, 1) == size(C, 1)
    @assert size(B, 2) == size(C, 2)

    if iszero(b)
        fill!(C, zero(T))
    else
        C .*= b
    end

    rows_a, cols_a = size(A)
    rows_b, cols_b = size(B)
    @inbounds for (ind_a, location_a) in enumerate(A.indices)
        @threads for i in eachindex(A.diags[ind_a])
            val_a = A.diags[ind_a][i] * a
            if location_a < rows_a
                index_ia = rows_a - location_a + i
                index_ja = i
            else
                index_ia = i
                index_ja = location_a - cols_a + i
            end
            min_loc = rows_b - index_ja + 1
            max_loc = 2 * rows_b - index_ja
            for (ind_b, location_b) in enumerate(B.indices)
                #index_ib = index_ja
                #       if ind < rows(A), then index = (rows - loc + i, i)
                #rows - loc + j = index_ja, j = index_ja - rows + loc
                #       else index = (i, loc - cols + i)
                # if location < rows(B), then
                if location_b <= rows_b && location_b >= min_loc
                    j = index_ja - rows_b + location_b
                    index_jb = j
                    val_b = B.diags[ind_b][j]
                    C[index_ia, index_jb] = muladd(val_a, val_b, C[index_ia, index_jb])
                elseif location_b > rows_b && location_b <= max_loc
                    j = index_ja
                    index_jb = location_b - cols_b + j
                    val_b = B.diags[ind_b][j]
                    C[index_ia, index_jb] = muladd(val_a, val_b, C[index_ia, index_jb])
                end
            end
        end
    end
    return C
end

export SparseBandedMatrix, setdiagonal!

@setup_workload begin
    # Minimal setup - create small test arrays
    @compile_workload begin
        # Precompile Float64 operations (most common)
        A = SparseBandedMatrix{Float64}(undef, 10, 10)
        A[1, 1] = 1.0
        A[5, 5] = 2.0
        _ = A[1, 1]
        setdiagonal!(A, [1.0, 2.0, 3.0], true)
        setdiagonal!(A, [4.0, 5.0], false)

        # Precompile mul! with Matrix (SparseBandedMatrix * Matrix)
        B = ones(10, 2)
        C = zeros(10, 2)
        mul!(C, A, B, 1.0, 0.0)

        # Precompile mul! from right (Matrix * SparseBandedMatrix)
        B2 = ones(2, 10)
        C2 = zeros(2, 10)
        mul!(C2, B2, A, 1.0, 0.0)

        # Precompile SparseBandedMatrix * SparseBandedMatrix -> Matrix
        A2 = SparseBandedMatrix{Float64}(undef, 10, 10)
        A2[1, 1] = 1.0
        setdiagonal!(A2, [1.0, 2.0], true)
        C3 = zeros(10, 10)
        mul!(C3, A, A2, 1.0, 0.0)

        # Precompile SparseBandedMatrix * SparseBandedMatrix -> SparseBandedMatrix
        C4 = SparseBandedMatrix{Float64}(undef, 10, 10)
        mul!(C4, A, A2, 1.0, 0.0)
    end
end

end
