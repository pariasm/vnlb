/* Modified work: Copyright (c) 2019, Pablo Arias <pariasm@gmail.com>
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * 
 * This program is free software: you can use, modify and/or redistribute it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version. You should have received a copy of this license
 * along this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 */
#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

#include <vector>
#include <string>

/* Invert a symmetric real matrix in place.
 *
 * mat: Symmetric input matrix. Output will be stored here.
 * N  : Matrix size
 *
 * returns 0 if normal exit, 1 if input matrix not positive definite
 */
int inverseMatrix(
	std::vector<float> &mat,
	const unsigned N);

/* Transpose a real square matrix in place.
 *
 * mat: Input matrix. Output will be stored here.
 * N  : Matrix size NxN.
 */
void transposeMatrix(
	std::vector<float> &mat,
	const unsigned N);

/* Compute empirical covariance matrix.
 *
 * data: n x d data matrix (each data point stored as a row)
 * cov : output covariance matrix
 * n   : number of data points
 * d   : dimension of data points
 */
void covarianceMatrix(
	std::vector<float> const& data,
	std::vector<float> &cov,
	const unsigned n,
	const unsigned d);

/* Matrix multiplication (matrices stored in row-major layout).
 *
 * AB: n x l product matrix (output)
 * A : n x m matrix
 * B : m x l matrix
 * n, m, l: size of matrices
 */
void productMatrix(
	std::vector<float> &AB,
	std::vector<float> const& A,
	std::vector<float> const& B,
	const unsigned n,
	const unsigned m,
	const unsigned l);

/* Matrix multiplication using BLAS SGEMM.
 *
 * AB: n x l product matrix (output)
 * A : n x m matrix
 * B : m x l matrix
 * n, m, l: size of matrices
 * transA: true for transposing A
 * transA: true for transposing B
 * colMajor: true if matrices are stored in row-major layout
 */
void productMatrix(
	std::vector<float> &AB,
	std::vector<float> const& A,
	std::vector<float> const& B,
	const unsigned n,
	const unsigned m,
	const unsigned l,
	const bool transA,
	const bool transB,
	const bool colMajor = true,
	unsigned lda = 0,
	unsigned ldb = 0);

/* Compute the r leading eigenvectors and eigenvalues of a
 * symmetric matrix.
 *
 * NOTES:
 * - uses LAPACKE_ssyevx
 * - matrices are stored in column-major ordering
 * - columns of input matrices are contiguous in memory
 * - only the upper triangular triangular part of mat is used
 * - WARNING: the upper triangular part of mat is destroyed
 * - the output U contains the eigenvectors as columns, and is
 *   stored in column-major ordering
 *
 * mat: nxn input matrix
 * n  : size of the matrix
 * r  : number of eigenvectors and eigenvalues
 * S  : vector with the r eigenvalues
 * U  : nxr matrix with eigenvectors
 *
 * returns LAPACKE_ssyevx error code.
 */
int matrixEigs(
	std::vector<float> &mat,
	const unsigned n,
	const unsigned r,
	std::vector<float> &S,
	std::vector<float> &U);

/* Write a matrix as an ascii file.
 *
 * mat: mxn input matrix
 * m,n: matrix size
 * filename: filename
 */
void printMatrix(
	std::vector<float> &mat,
	unsigned m,
	unsigned n,
	std::string filename);

#endif // LIB_MATRIX_H_INCLUDED
