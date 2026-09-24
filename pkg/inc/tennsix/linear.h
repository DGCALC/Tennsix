//
// Created by Lucius Schoenbaum on 11/6/21.
//

// Convenience header for grabbing
// types for basic linear algebra types: vector, matrix, column matrix,
// over your choice of scalar T.

#pragma once

#include <tennsix/super.h>
#include <scalar/754.h>

// scalar-independent array Vector<T, n>
template<typename T, size_t size, typename PP = T>
using Vector = Super<T, T, size, 1>;

// linear maps over Vector --
// - index-references using subscript operator:
//   A[{i,j}] is the (ij)th element, counting from zero
// - multiply by a vector or matrix using *, as in:
//   A * x = b, A*B*x = b, etc.
template<typename T, size_t size, size_t cosize = size>
using Matrix = Super<T, T, size, cosize>;

// Be careful with this, but --
// - the object is a Vector of Vectors (as you see here),
//    and indexes work accordingly
// - multiply a column matrix with a vector using %, as in:
//    A % x = b
//   ....this is the ordinary matrix multiplication of A and x
template<typename T, size_t size, size_t cosize>
using ColumnMatrix = Super<Vector<T, size>, T, cosize, 1>;

// three-index matrix (3d matrix)
// with index pattern A[i][{j,k}], which looks quite fine (to me at least).
// layout identical to A[i][j][k], row-major.
template<typename T, size_t size, size_t cosize, size_t trisize>
using TriMatrix = Vector<Matrix<T, cosize, trisize>, size>;

// four-index matrix (4d matrix)
// with index pattern A[{i,j}][{k,l}], which looks quite fine (to me at least).
// layout identical to A[i][j][k][l], row-major.
template<typename T, size_t size, size_t cosize, size_t trisize, size_t quadsize>
using QuadMatrix = Matrix<Matrix<T, trisize, quadsize>, size, cosize>;

// five-index matrix (5d matrix)
// with index pattern A[i][{j,k}][{l,m}], which is getting kind of ridiculous, perhaps.
// layout identical to A[i][j][k][l][m], row-major.
template<typename T, size_t size, size_t cosize, size_t trisize, size_t quadsize, size_t pentasize>
using PentaMatrix = Vector<QuadMatrix<T, cosize, trisize, quadsize, pentasize>, size>;

// six-index matrix (6d matrix)
// with index pattern A[{i,j}][{k,l}][{m,n}], which is getting kind of ridiculous, perhaps.
// layout identical to A[i][j][k][l][m][n], row-major.
template<typename T, size_t size, size_t cosize, size_t trisize, size_t quadsize, size_t pentasize, size_t hexasize>
using HexaMatrix = Matrix<QuadMatrix<T, trisize, quadsize, pentasize, hexasize>, size, cosize>;


// 2d and 3d scalar vectors

using PaperVector = Vector<double754, 2>;
using PaperMatrix = Matrix<double754, 2, 2>;

using GeoVector = Vector<double754, 3>;
using GeoMatrix = Matrix<double754, 3, 3>;




