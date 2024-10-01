#include "nr.h"
#include "nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

//
// *** DCT Type I ***
//

void create_f_ext_dct1(float ***data, float ***f_ext, int n) {

	int i, j;
	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			f_ext[1][i][j] = data[1][i][j];
			if (j == 1 || j == n) f_ext[1][i][j] *= sqrt(2);
			if (i == 1 || i == n) f_ext[1][i][j] *= sqrt(2);
		}
	}

	for (i = 1; i <= n; i++) {
		for (j = n+1; j <= 2*n-2; j++) {
			f_ext[1][i][j] = data[1][i][2*n-j];
			if (j == 1 || j == n) f_ext[1][i][j] *= sqrt(2);
			if (i == 1 || i == n) f_ext[1][i][j] *= sqrt(2);
		}
	}

	for (i = n+1; i <= 2*n-2; i++) {
		for (j = 1; j <= n; j++) {
			f_ext[1][i][j] = data[1][2*n-i][j];
			if (j == 1 || j == n) f_ext[1][i][j] *= sqrt(2);
			if (i == 1 || i == n) f_ext[1][i][j] *= sqrt(2);
		}
	}

	for (i = n+1; i <= 2*n-2; i++) {
		for (j = n+1; j <= 2*n-2; j++) {
			f_ext[1][i][j] = data[1][2*n-i][2*n-j];
		}
	}
	
}

void dct1_2d(float ***data, int n) {

	float ***f_ext = f3tensor(1, 1, 1, 2*n-2, 1, 2*n-2);

	create_f_ext_dct1(data, f_ext, n);	

	float **speq = matrix(1, 1, 1, 2*(2*n-2));

	rlft3(f_ext, speq, 1, 2*n-2, 2*n-2, 1);

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			data[1][i][j] = f_ext[1][i][j*2-1] * 0.5 / (n-1);
			if (j == 1 || j == n) data[1][i][j] /= sqrt(2);
			if (i == 1 || i == n) data[1][i][j] /= sqrt(2);
		}
		data[1][i][n] = speq[1][2*i-1] * 0.5 / (n-1) / sqrt(2);
		if (i == 1 || i == n) data[1][i][n] /= sqrt(2);
	}

	free_matrix(speq, 1, 1, 1, 2*(2*n-2));
	
}

//
// *** DCT Type II ***
//

void create_f_ext_dct2(float ***data, float ***f_ext, int n) {
	
	for (int i = 2; i <= 2*n; i+=2) {
		for (int j = 2; j <= 2*n; j+=2) {
			f_ext[1][i][j] = data[1][i/2][j/2];
		}
	}

	for (int i = 2; i <= 2*n; i+=2) {
		for (int j = 2*n+2; j <= 4*n; j+=2) {
			f_ext[1][i][j] = data[1][i/2][n-(j-2*n)/2+1];
		}
	}

	for (int i = 2*n+2; i <= 4*n; i+=2) {
		for (int j = 2; j <= 2*n; j+=2) {
			f_ext[1][i][j] = data[1][n-(i-2*n)/2+1][j/2];
		}
	}
	
	for (int i = 2*n+2; i <= 4*n; i+=2) {
		for (int j = 2*n+2; j <= 4*n; j+=2) {
			f_ext[1][i][j] = data[1][n-(i-2*n)/2+1][n-(j-2*n)/2+1];
		}
	}

}

void dct2_2d(float ***data, int n) {

	float ***f_ext = f3tensor(1, 1, 1, n*4, 1, n*4);

	create_f_ext_dct2(data, f_ext, n);

	float **speq = matrix(1, 1, 1, 2*4*n);

	rlft3(f_ext, speq, 1, 4*n, 4*n, 1);

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= 2*n; j+=2) {
			data[1][i][(j-1)/2+1] = f_ext[1][i][j] / (2 * n);
			// Orthogonalize by dividing the first row + column by sqrt(2)
			if (i == 1) data[1][i][(j-1)/2+1] /= sqrt(2);
			if (j == 1) data[1][i][(j-1)/2+1] /= sqrt(2);
		}
	}

	free_matrix(speq, 1, 1, 1, 2*4*n);
	free_f3tensor(f_ext, 1, 1, 1, n*4, 1, n*4);
	
}

//
// *** DCT Type III ***
//

void create_f_ext_dct3_1d(float *data, float *f_ext, int n) {

	float complex temp;
	
	for (int i = 1; i <= n; i++) {

		temp = data[i] * 2.0 * csqrt(n/2.0) * cexp((i-1) * M_PI / (2.0 * n) * I);
		// Orthogonalize by multiplying the first element by sqrt(2)
		if (i == 1) temp *= csqrtf(2.0);
		f_ext[2*i-1] = crealf(temp);
		f_ext[2*i] = cimagf(temp);

		
	}

}

void dct3_1d(float *data, int n) {

	float *f_ext = vector(1, 2*n);

	create_f_ext_dct3_1d(data, f_ext, n);

	realft(f_ext, 2*n, -1);

	// Numerical recipes stores the result of realfft horribly:
	// the first 2 elements are the first and last real elements of the result, 
	// and the n-2 elements afterwards are elements 1 through n-1
	data[1] = f_ext[1] / n;
	data[n] = f_ext[2] / n;
	for (int i = 2; i <= n-1; i++) {
		data[i] = f_ext[i+1] / n; // Normalize and cut out the real part
	}

	free_vector(f_ext, 1, 2*n);

}

void dct3_2d(float ***data, int n) {

	int i, j;

	for (i = 1; i <= n; i++) {
		dct3_1d(data[1][i], n);
	}

	float *data_temp = vector(1, n);

	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			data_temp[j] = data[1][j][i];
		}

		dct3_1d(data_temp, n);

		for (j = 1; j <= n; j++) {
			data[1][j][i] = data_temp[j];
		}

	}

	free_vector(data_temp, 1, n);
		
}

//
// *** DST Type I ***
//

void create_f_ext_dst1(float ***data, float ***f_ext, int n) {

	int i,j;

	for (i = 2; i <= n+1; i++) {
		for (j = 2; j <= n+1; j++) {
			f_ext[1][i][j] = data[1][i-1][j-1];
		}
	}	

	for (i = n+3; i <= 2*n+2; i++) {
		for (j = 2; j <= n+1; j++) {
			f_ext[1][i][j] = -data[1][n-(i-n-2)+1][j-1];
		}
	}	

	for (i = 2; i <= n+1; i++) {
		for (j = n+3; j <= 2*n+2; j++) {
			f_ext[1][i][j] = -data[1][i-1][n-(j-n-2)+1];
		}
	}

	for (i = n+3; i <= 2*n+2; i++) {
		for (j = n+3; j <= 2*n+2; j++) {
			f_ext[1][i][j] = data[1][n-(i-n-2)+1][n-(j-n-2)+1];
		}
	}	

}

// Note: Requires no additional orthogonalization since the DST is already
// orthogonal by default
void dst1_2d(float ***data, int n) {

	float **speq = matrix(1, 1, 1, 4*n+4);
	float ***f_ext = f3tensor(1, 1, 1, 2*n+2, 1, 2*n+2);

	create_f_ext_dst1(data, f_ext, n);

	rlft3(f_ext, speq, 1, 2*n+2, 2*n+2, 1);

	for (int i = 2; i <= n+1; i++) {
		for (int j = 3; j <= 2*n+1; j+=2) {
			data[1][i-1][(j-1)/2] = -0.5 / (n+1) * f_ext[1][i][j];
		}	
	}

	free_matrix(speq, 1, 1, 1, 4*n+4);
	free_f3tensor(f_ext, 1, 1, 1, 2*n+2, 1, 2*n+2);

}

//
// *** DST Type II ***
//

void create_f_ext_dst2(float ***data, float ***f_ext, int n) {

	int i, j;

	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			if ((i % 2 == 0) != (j % 2 == 0)) f_ext[1][i][j] = -data[1][i][j];
			else f_ext[1][i][j] = data[1][i][j];
		}
	}	

}

void dst2_2d(float ***data, int n) {

	float ***f_ext = f3tensor(1, 1, 1, n, 1, n);

	create_f_ext_dst2(data, f_ext, n);
	
	// Since the DCT is already orthogonalized, the resulting DST is too
	dct2_2d(f_ext, n);

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			data[1][i][j] = f_ext[1][n-i+1][n-j+1];
		}	
	}

	free_f3tensor(f_ext, 1, 1, 1, n, 1, n);

}

//
// *** DST Type III ***
//

void create_f_ext_dst3(float ***data, float ***f_ext, int n) {

	int i, j;

	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			f_ext[1][i][j] = data[1][i][j];
		}
	}	

}

void dst3_2d(float ***data, int n) {

	float ***f_ext = f3tensor(1, 1, 1, n, 1, n);

	create_f_ext_dst3(data, f_ext, n);

	// Since the DCT is already orthogonalized, the resulting DST is too
	dct3_2d(f_ext, n);

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			if ((i % 2 == 0) != (j % 2 == 0)) data[1][i][j] = -f_ext[1][i][j];
			else data[1][i][j] = f_ext[1][i][j];
		}	
	}

	free_f3tensor(f_ext, 1, 1, 1, n, 1, n);

}
