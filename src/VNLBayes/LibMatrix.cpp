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

/* LibMatrix.cpp
 * Tools for matrix manipulation.
 *
 * @author Pablo Arias <pariasm@gmail.com>
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 */

#include "LibMatrix.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

extern "C" {
#include <cblas.h>
#include <lapacke.h>
}

using namespace std;

int inverseMatrix(
	vector<float> &mat,
	const unsigned N)
{
	// Initializations
	float z;
	unsigned p, q, r, s, t, j, k;

	// Implementation based on ccmath functions by Daniel A. Atkinson.

	for (j = 0, p = 0; j < N ; j++, p += N + 1) {
		for (q = j * N; q < p ; q++)
			mat[p] -= mat[q] * mat[q];

		if (mat[p] <= 0.f) return 1;

		mat[p] = sqrtf(mat[p]);

		for (k = j + 1, q = p + N; k < N ; k++, q += N) {
			for (r = j * N, s = k * N, z = 0.f; r < p; r++, s++) {
				z += mat[r] * mat[s];
			}

			mat[q] -= z;
			mat[q] /= mat[p];
		}
	}

	transposeMatrix(mat, N);

	for (j = 0, p = 0; j < N; j++, p += N + 1) {
		mat[p] = 1.f / mat[p];

		for (q = j, t = 0; q < p; t += N + 1, q += N) {
			for (s = q, r = t, z = 0.f; s < p; s += N, r++) {
				z -= mat[s] * mat[r];
			}

			mat[q] = z * mat[p];
		}
	}

	for (j = 0, p = 0; j < N; j++, p += N + 1) {
		for (q = j, t = p - j; q <= p; q += N, t++) {
			for (k = j, r = p, s = q, z = 0.f; k < N; k++, r++, s++)
				z += mat[r] * mat[s];

			mat[t] = mat[q] = z;
		}
	}

	return 0;
}

void transposeMatrix(
	vector<float> &mat,
	const unsigned N)
{
	// Implementation based on ccmath functions by Daniel A. Atkinson.

	for (unsigned i = 0; i < N - 1; i++) {
		unsigned p = i * (N + 1) + 1;
		unsigned q = i * (N + 1) + N;

		for (unsigned j = 0; j < N - 1 - i; j++, p++, q += N) {
			const float s = mat[p];
			mat[p] = mat[q];
			mat[q] = s;
		}
	}
}

void covarianceMatrix(
	vector<float> const& data,
	vector<float> &cov,
	const unsigned n,
	const unsigned d)
{
	const float coefNorm = 1.f / (float) (n);

	for (unsigned i = 0; i < d; i++)
	for (unsigned j = 0; j < i + 1; j++) {
		float val = 0.f;
		for (unsigned k = 0; k < n; k++)
			val += data[i * n + k] * data[j * n + k];

		cov[i * d + j] = val * coefNorm;
		cov[i + j * d] = val * coefNorm;
	}
}

void productMatrix(
	vector<float> &AB,
	vector<float> const& A,
	vector<float> const& B,
	const unsigned n,
	const unsigned m,
	const unsigned l)
{
	// Implementation based on ccmath functions by Daniel A. Atkinson.

	vector<float> q0(m, 0.f);

	for (unsigned i = 0; i < l; i++) {
		for (unsigned k = 0; k < m; k++) {
			q0[k] = B[i + k * l];
		}

		for (unsigned j = 0; j < n; j++) {
			float z = 0.f;

			for (unsigned k = 0; k < m; k++) {
				z += A[j * m + k] * q0[k];
			}

			AB[i + j * l] = z;
		}
	}
}

void productMatrix(
	vector<float> &AB,
	vector<float> const& A,
	vector<float> const& B,
	const unsigned n,
	const unsigned l,
	const unsigned m,
	const bool transA,
	const bool transB,
	const bool colMajor,
	unsigned lda,
	unsigned ldb)
{
	if (!lda) lda = colMajor ? (transA ? m : n) : (transA ? n : m);
	if (!ldb) ldb = colMajor ? (transB ? l : m) : (transB ? m : l);

	cblas_sgemm(colMajor ? CblasColMajor : CblasRowMajor,  // matrix storage mode
	            transA   ? CblasTrans    : CblasNoTrans,   // op(A)
	            transB   ? CblasTrans    : CblasNoTrans,   // op(B)
	            n,                                         // rows(op(A)) [= rows(AB)   ]
	            l,                                         // cols(op(B)) [= cols(AB)   ]
	            m,                                         // cols(op(A)) [= rows(op(B))]
	            1.f,                                       // alpha
	            A.data(), lda,                             // A, lda
	            B.data(), ldb,                             // B, ldb
	            0.f,                                       // beta
	            AB.data(), colMajor ? n : l);              // AB, ldab
}

int matrixEigs(
	vector<float> &mat,
	const unsigned n,
	const unsigned r,
	vector<float> &S,
	vector<float> &U)
{
	// set parameters for LAPACKE SSYEV function
	// SSYEVX: Single SYmmetric EigenValues and eigenvectors eXpert
	lapack_int m;        // total values of eigenvalues found

	S.resize(n);         // array of dimension n. The first m entries
	                     // contain the eigenvalues found in ascending order.

	U.resize(n*r);       // the first m columns contain the output eigenvectors

	lapack_int lda = n;  // leading dimension for input matrix
	lapack_int ldu = n;  // leading dimension for output eigenvectors matrix

	lapack_int ifail[n]; // if jobz == 'V' and info > 0, then ifail contains
	                     // the indices of the eigenvectors that didn't converge

	lapack_int info;     // info =  0 : successful exit
	                     // info = -i : ith argument is wrong
	                     // info =  i : i eigenvectors failed to converge

	info = LAPACKE_ssyevx(LAPACK_COL_MAJOR,
			'V',          // compute eigenvalues and eigenvectors
			'I',          // range 'I': eigenvals/vecs between IL-th and IU-th
			'U',          // use upper triangular part of A
			n,            // order of matrix A
			mat.data(),   // matrix A
			lda,          // stride of matrix
			-1, -1,       // not used (used only when range is 'V'
			n - r + 1, n, // IL and IU indices for eigenvals/vecs range
			0,            // abstol for stopping criterion
			&m,           // total values of eigenvectors found
			S.data(),     // eigenvalues output
			U.data(),     // eigenvectors matrix output
			ldu,          // eigenvectors matrix stride
			ifail         // eigenvectors that did not converge
			);

	return(info);
}

void printMatrix(
	std::vector<float> &matrix,
	unsigned rows,
	unsigned cols,
	std::string filename)
{
	FILE *file = fopen(filename.c_str(),"w");

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			fprintf(file, "%.16g ", matrix[i*cols + j]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}
