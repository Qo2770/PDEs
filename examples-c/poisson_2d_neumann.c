#include "nr.h"     // Include Numerical Recipes FFT functions
#include "nrutil.h" // Include Numerical Recipes utility functions
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4 // Size of the grid (must be a power of 2)
#define L 1.0 // Physical size of the domain
#define SIG 0.2 // Width << 1

void initialize_rhs(float ***rhs, int n, float l) {
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            float x = (i - 1) * l / n;
            float y = (j - 1) * l / n;
			float rsq = (x - 0.5 * L ) * (x - 0.5 * L) + (y - 0.5 * L ) * (y - 0.5 * L);
			rhs[1][i][j] = exp(-rsq / (2 * SIG * SIG)) * (rsq - 2 * SIG * SIG) / (SIG * SIG * SIG * SIG);
        }
    }
}

void create_f_ext(float ***data, float ***f_ext, int n) {
	
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

void cosine_fft_2d(float ***data, int n) {

	float ***f_ext = f3tensor(1, 1, 1, n*4, 1, n*4);

	create_f_ext(data, f_ext, n);

	float **speq = matrix(1, 1, 1, 2*4*n);

	for (int i = 1; i <= N*4; i++) {
		for (int j = 1; j <= N*4; j++) {
			printf("%f ", f_ext[1][i][j]);
		}
		printf("\n");
	}
	printf("f_ext\n");
	printf("\n");

	rlft3(f_ext, speq, 1, 4*N, 4*N, 1);

	for (int i = 1; i <= N*4; i++) {
		for (int j = 1; j <= N*4; j++) {
			printf("%f ", f_ext[1][i][j]);
		}
		printf("\n");
	}
	printf("f_ext2\n");

	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= 2*n; j+=2) {
			data[1][i][(j-1)/2+1] = f_ext[1][i][j] * 2/n * 1/4;
			if (i == 1) data[1][i][(j-1)/2+1] = data[1][i][(j-1)/2+1] * sqrt(0.5);
			if (j == 1) data[1][i][(j-1)/2+1] = data[1][i][(j-1)/2+1] * sqrt(0.5);
		}
	}

	free_matrix(speq, 1, 1, 1, 2*4*n);
	
}

int main() {

	int i, j;

	float ***data, **speq;
	data = f3tensor(1, 1, 1, N, 1, N);
	
	float **phi = matrix(1, N, 1, N);   // Solution in real space

	// Initialize the data
	initialize_rhs(data, N, L);

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			printf("%f ", data[1][i][j]);
		}
		printf("\n");
	}
	printf("\n");


	// Perform forward sine FFT
	cosine_fft_2d(data, N);

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			printf("%f ", data[1][i][j]);
		}
		printf("\n");
	}
	printf("\n");

	// Solve Poisson equation in Fourier space
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {

			float kx = (i-1 <= N / 2) ? i-1 : (i-1) - N;
			float ky = (j-1 <= N / 2) ? j-1 : (j-1) - N;

			// Square and multiply with the square of the prefix 2*pi/L
			float k_squared = -(4 * (M_PI * M_PI) / (L * L)) * ( kx * kx + ky * ky );

			// At 0,0 don't normalize
			if (k_squared != 0) {
				data[1][i][j] /= k_squared;
			}
		}
	}
	
	// Perform inverse sine FFT
	cosine_fft_2d(data, N);

	// Copy data array back to phi
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			phi[i][j] = data[1][i][j] * 4 / (N * N); // Normalize the solution
		}
	}

	// Print the solution values
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			printf("%f ", phi[i][j] - phi[1][1]);
		}
		printf("\n");
	}

	// Clean up
	free_matrix(phi, 1, N, 1, N);
	free_f3tensor(data, 1, 1, 1, N, 1, N);

	return 0;
}
