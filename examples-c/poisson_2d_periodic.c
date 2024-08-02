#include "nr.h"     // Include Numerical Recipes FFT functions
#include "nrutil.h" // Include Numerical Recipes utility functions
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 8 // Size of the grid (must be a power of 2)
#define L 1.0 // Physical size of the domain
#define SIG 0.1 // Width << 1

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

int main() {

	int i, j;

	float ***data, **speq;
	data = f3tensor(1, 1, 1, N, 1, N);
	speq = matrix(1, 1, 1, 2*N);
	
	float **phi = matrix(1, N, 1, N);   // Solution in real space

	// Initialize the data
	initialize_rhs(data, N, L);

	// Perform forward FFT
	rlft3(data, speq, 1, N, N, 1);

	// Solve Poisson equation in Fourier space
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N/2; j++) {

			float kx = (i-1 <= N / 2) ? i-1 : (i-1) - N;
			float ky = (j-1 <= N / 2) ? j-1 : (j-1) - N;

			// Square and multiply with the square of the prefix 2*pi/L
			float k_squared = -(4 * (M_PI * M_PI) / (L * L)) * ( kx * kx + ky * ky );

			// At 0,0 don't normalize
			if (k_squared != 0) {
				data[1][i][2*(j-1)+1] /= k_squared;
				data[1][i][2*(j-1)] /= k_squared;
			}
		}
	}
	
	// Compute the last row of the FFT, since that is cut out
	// to make room for the complex coeff and is in a seperate 
	// array
	for (i = 1; i <= N; i++) {
		float kx = (i-1 <= N / 2) ? i-1 : (i-1) - N;
		float ky = (j-1 <= N / 2) ? j-1 : (j-1) - N;
		float k_squared = -(4 * (M_PI * M_PI) / (L * L)) * ( kx * kx + ky * ky );
		speq[1][2*(i-1)+1] /= k_squared;
		speq[1][2*(i-1)] /= k_squared;	
	}

	// Perform inverse FFT
	rlft3(data, speq, 1, N, N, -1);

	// Copy data array back to phi
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			phi[i][j] = data[1][i][j] * 2 / (N * N); // Normalize the solution
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
	free_matrix(speq, 1, 1, 1, 2 * N);
	free_f3tensor(data, 1, 1, 1, N, 1, N);

	return 0;
}
