#include "nr.h"     // Include Numerical Recipes FFT functions
#include "nrutil.h" // Include Numerical Recipes utility functions
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "transforms.h"

#define N 8 // Size of the grid (2*N+2 must be a power of 2)
#define L 1.0 // Physical size of the domain
#define SIG 0.2 // Width << 1

void initialize_rhs(float ***rhs, int n, float l) {
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            float x = (i - 0.5) * l / n;
            float y = (j - 0.5) * l / n;
			float rsq = (x - 0.5 * L ) * (x - 0.5 * L) + (y - 0.5 * L ) * (y - 0.5 * L);
			rhs[1][i][j] = exp(-rsq / (2 * SIG * SIG)) * (rsq - 2 * SIG * SIG) / (SIG * SIG * SIG * SIG);
        }
    }
}

int main() {

	int i, j;

	float ***data;
	data = f3tensor(1, 1, 1, N, 1, N);
	
	float **phi = matrix(1, N, 1, N);   // Solution in real space

	// Initialize the data
	initialize_rhs(data, N, L);

	// Perform forward sine FFT
	dst2_2d(data, N);

	// Solve Poisson equation in Fourier space
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {

			float kx = i;
			float ky = j;

			// Square and multiply with the square of the prefix 2*pi/L
			float k_squared = -((M_PI * M_PI) / (L * L)) * ( kx * kx + ky * ky );

			// At 0,0 don't normalize
			if (k_squared != 0) {
				data[1][i][j] /= k_squared;
			} 
		}
	}
	// Set 0,0 to zero
	data[1][1][1] = 0.0;

	// Perform inverse sine FFT
	dst2_2d(data, N);

	// Copy data array back to phi
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			phi[i][j] = data[1][i][j]; 
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
