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

void cosine_fft_2d(float ***data, int n) {

	for (int i = 1; i <= n+1; i++) {

		cosft1(data[1][i], n);

	}
	
	float *temp_data = vector(1, n+1);
	for (int i = 1; i <= n+1; i++) {

		for (int j = 1; j <= n+1; j++) {
			
			temp_data[j] = data[1][j][i];

		}

		cosft1(temp_data, n);

		for (int j = 1; j <= n+1; j++) {
			
			data[1][j][i] = temp_data[j];

		}

	}

	free_vector(temp_data, 1, n+1);

}

int main() {

	int i, j;

	float ***data, **speq;
	data = f3tensor(1, 1, 1, N+1, 1, N+1);
	
	float **phi = matrix(1, N, 1, N);   // Solution in real space

	// Initialize the data
	initialize_rhs(data, N, L);

	// Perform forward sine FFT
	cosine_fft_2d(data, N);

	// Solve Poisson equation in Fourier space
	for (i = 1; i <= N+1; i++) {
		for (j = 1; j <= N+1; j++) {

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
