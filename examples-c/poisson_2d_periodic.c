#include "nr.h"     // Include Numerical Recipes FFT functions
#include "nrutil.h" // Include Numerical Recipes utility functions
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4 // Size of the grid (must be a power of 2)
#define L 1.0 // Physical size of the domain
#define SIG 0.2 // Width << 1

void initialize_rhs(float *rhs, int n, float l) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2*n; j+=2) {
            float x = i * l / n;
            float y = (j/2.0) * l / n;
			float rsq = (x - 0.5 * L) * (x - 0.5 * L) + (y - 0.5 * L) * (y - 0.5 * L);
			int idx = i*2*n+j+1;
			rhs[idx] = exp(-rsq / (2 * SIG * SIG)) * (rsq - 2 * SIG * SIG) / (SIG * SIG * SIG * SIG);
        }
    }
}

int main() {

	int i, j;

	float *data;
	
	data = vector(1, 2*(N*N)+1); // Flattened data array, 2*N^2 to hold complex values, 1-indexed
	unsigned long nn[] = {0, N, N}; // Size of all dims (uses 1-indexing so first element is padding)
	
	float **phi = matrix(1, N, 1, N);   // Solution in real space

	// Initialize the data
	initialize_rhs(data, N, L);

	// Perform forward FFT
	fourn(data, nn, 2, 1);

	// Solve Poisson equation in Fourier space
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {

			float kx = (i <= N / 2) ? i : i - N;
			float ky = (j <= N / 2) ? j : j - N;

			// Square and multiply with the square of the prefix 2*pi/L
			float k_squared = -(4 * (M_PI * M_PI) / (L * L)) * ( kx * kx + ky * ky );

			// At 0,0 don't normalize
			if (k_squared != 0) {
				int idx = i*2*N+2*(j+1);
				data[idx-1] /= k_squared;
				data[idx] /= k_squared;
			}
		}
	}
	// Zero out first component
	data[1] = 0.0;

	// Perform inverse FFT
	fourn(data, nn, 2, -1);

	// Copy data array back to phi
	for (i = 0; i < N; i++) {
		for (j = 0; j < 2*N; j+=2) {
			phi[i+1][j/2+1] = data[i*2*N+j+1] * 1 / (N * N); // Normalize the solution
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
	free_vector(data, 1, 2*(N*N)+1);

	return 0;
}
