#include "nr.h"     // Include Numerical Recipes FFT functions
#include "nrutil.h" // Include Numerical Recipes utility functions
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

#define N 15 // Size of the grid (2*N+2 must be a power of 2)
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

	for (int i = 2; i <= n+1; i++) {
		for (int j = 2; j <= n+1; j++) {
			f_ext[1][i][j] = data[1][i-1][j-1];
		}
	}	

	for (int i = n+3; i <= 2*n+2; i++) {
		for (int j = 2; j <= n+1; j++) {
			f_ext[1][i][j] = -data[1][n-(i-n-2)+1][j-1];
		}
	}	

	for (int i = 2; i <= n+1; i++) {
		for (int j = n+3; j <= 2*n+2; j++) {
			f_ext[1][i][j] = -data[1][i-1][n-(j-n-2)+1];
		}
	}

	for (int i = n+3; i <= 2*n+2; i++) {
		for (int j = n+3; j <= 2*n+2; j++) {
			f_ext[1][i][j] = data[1][n-(i-n-2)+1][n-(j-n-2)+1];
		}
	}	

}

void sine_fft_2d(float ***data, int n, int isign) {

	float **speq = matrix(1, 1, 1, 4*n+4);
	float ***f_ext = f3tensor(1, 1, 1, 2*n+2, 1, 2*n+2);

	create_f_ext(data, f_ext, n);

	rlft3(f_ext, speq, 1, 2*n+2, 2*n+2, 1);

	for (int i = 2; i <= 1*n+1; i++) {
		for (int j = 3; j <= 2*n+2; j+=2) {
			if (isign > 0) data[1][i-1][(j-1)/2] = -f_ext[1][i][j] / 4;
			else data[1][i-1][(j-1)/2] = -f_ext[1][i][j] * 4 / pow(N+1, 2.0);
		}	
	}

	free_matrix(speq, 1, 1, 1, 4*n+4);
	free_f3tensor(f_ext, 1, 1, 1, 2*n+2, 1, 2*n+2);

}

int main() {

	int i, j;

	float ***data;
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
	sine_fft_2d(data, N, 1);

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

	data[1][1][1] = 0.0;

	// Perform inverse sine FFT
	sine_fft_2d(data, N, -1);

	// Copy data array back to phi
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			phi[i][j] = data[1][i][j] * 4 / ((N+1)*(N+1)); // Normalize the solution
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
