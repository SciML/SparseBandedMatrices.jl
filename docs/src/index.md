# SparseBandedMatrices.jl Documentation

## Overview

`SparseBandedMatrices.jl` provides a fast implementation of sparse banded matrices in Julia. The package is optimized for matrices where only certain diagonals contain non-zero elements, which is a common pattern in numerical linear algebra.

This package was primarily developed for use in butterfly LU factorizations implemented in [RecursiveFactorization.jl](https://github.com/JuliaLinearAlgebra/RecursiveFactorization.jl) and [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl), but can be useful for any application involving sparse banded matrices.

### Key Features

- **Efficient storage**: Only non-zero diagonals are stored, reducing memory usage
- **Fast multiplication**: Optimized matrix-matrix and matrix-vector multiplication
- **Simple interface**: Standard Julia array indexing and operations
- **Thread-safe**: Parallel operations using Julia's threading capabilities

## Installation

```julia
using Pkg
Pkg.add("SparseBandedMatrices")
```

## Quick Start

```julia
using SparseBandedMatrices

# Create an empty 5×5 sparse banded matrix
A = SparseBandedMatrix{Float64}(undef, 5, 5)

# Set individual elements
A[1,1] = 5.0
A[2,1] = 3.0
A[1,2] = 4.0

# Set an entire diagonal at once (more efficient)
setdiagonal!(A, [3.0, 4.0, 5.0], true)
```

## Usage Examples

### Creating Matrices

There are two main ways to create a `SparseBandedMatrix`:

#### Empty Matrix

Create an uninitialized matrix and populate it element-by-element or diagonal-by-diagonal:

```julia
A = SparseBandedMatrix{Float64}(undef, 5, 5)
A[1,1] = 5.0
setdiagonal!(A, [3.0, 4.0, 5.0], true)
```

#### Pre-specified Diagonals

Create a matrix with known diagonals upfront:

```julia
# Create a 6×6 matrix with two diagonals
# The constructor takes:
#   - A vector of diagonal indices [1, 8]
#   - A vector of vectors containing the diagonal values
#   - Matrix dimensions (rows, columns)
B = SparseBandedMatrix{Float64}([1, 8], [[3.0], [-2.0, 5.0, 1.0, 3.0]], 6, 6)
```

**Note:** Diagonal indices are stored internally using a specific indexing scheme where for an N×M matrix, diagonal indices range from 1 to N+M-1. The length of each diagonal vector must match the actual length of that diagonal in the matrix.

### Matrix Operations

`SparseBandedMatrix` supports standard Julia array operations:

```julia
# Matrix-vector multiplication
v = rand(5)
result = A * v

# Matrix-matrix multiplication
B = rand(5, 3)
C = A * B

# Element access
value = A[2, 3]
A[2, 3] = 7.0

# Size queries
rows, cols = size(A)
```

### Using `setdiagonal!`

The `setdiagonal!` function efficiently sets an entire diagonal:

```julia
A = SparseBandedMatrix{Float64}(undef, 5, 5)

# Set a lower diagonal
# The `true` parameter indicates this is a lower diagonal
setdiagonal!(A, [3.0, 4.0, 5.0], true)

# Set an upper diagonal
# The `false` parameter indicates this is an upper diagonal
setdiagonal!(A, [1.0, 2.0], false)
```

## Understanding Sparse Banded Matrices

A sparse banded matrix is a matrix where only certain diagonals contain non-zero elements. Unlike dense banded matrices (which store all elements within a band), sparse banded matrices only store the diagonals that actually contain non-zero values, regardless of how far apart they are.

This makes `SparseBandedMatrix` particularly efficient for:
- Butterfly factorizations where non-zero elements appear on specific, potentially distant diagonals
- Structured problems where the pattern of non-zero diagonals is known in advance
- Cases where traditional banded storage would waste memory on zero diagonals

## Performance Considerations

- Use `setdiagonal!` when possible instead of setting individual elements, as it's more efficient for initializing diagonals
- Pre-specify diagonals during construction when you know them in advance to avoid repeated memory reallocations
- The package is optimized for multiplication operations, making it ideal for iterative solvers and matrix factorizations
- Thread-parallel operations are automatically used for matrix multiplications when Julia threading is enabled

## API Reference

```@docs
SparseBandedMatrix
setdiagonal!
```