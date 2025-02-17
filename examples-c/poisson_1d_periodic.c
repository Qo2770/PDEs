#include "nr.h"     // Include Numerical Recipes FFT functions
#include "nrutil.h" // Include Numerical Recipes utility functions
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 16 // Size of the grid (must be a power of 2)
#define L 1.0 // Physical size of the domain

//
// Adapted from Brian's Python code
//

void initialize_rhs(float *rhs) {
    for (int i = 0; i < N; i++) {
		float x = ((float)i)/N;
		float pi_sq = M_PI * M_PI;
		rhs[2*i+1] = -16 * pi_sq * sin(4 * M_PI * x) - 8 * pi_sq * pow(sin(2 * M_PI * x), 2.0) + 8 * pi_sq * pow(cos(2 * M_PI * x), 2.0);
    }
}

void construct_ksq(float *ksq) {

	int i;
	float k;

	for (i = 0; i <= N; i++) {
		
		k = (i <= N / 2) ? i : i - N;
		
		// Square and multiply with the square of the prefix 2*pi/L
		ksq[i] = -(4 * (M_PI * M_PI) / (L * L)) * ( k * k );

	}

} 

void apply_ksq(float *rhs, float *ksq) {
	
	int i;
	
	for (i = 0; i <= N; i++) {
		
		// At 0,0 don't normalize
		if (ksq[i] != 0) {
			rhs[2*i+1] /= ksq[i];
			rhs[2*i+2] /= ksq[i]; 
		} 
	}

	// Zero out first component
	rhs[1] = 0.0;
	rhs[2] = 0.0;

}

int main() {

	int i, j;

	float *data;
	float *ksq;
	
	data = vector(1, 2*N); // Flattened data array, 2*N^2 to hold complex values, 1-indexed
	ksq = vector(0, N-1); // Flattened data array, 2*N^2 to hold complex values, 1-indexed
	unsigned long nn = N; // Size of all dims (uses 1-indexing so first element is padding)
	
	float *phi = vector(1, N);   // Solution in real space

	// Initialize the data
	initialize_rhs(data);

	for (i = 1; i <= 2*N; i+=2) {
		printf("%f ", data[i]);
	}
	printf("\n");

	// Perform forward FFT
	four1(data, nn, 1);
	
	// Construct solver for the Poisson equation in Fourier space
	construct_ksq(ksq);

	// Apply the constructed solver to the input
	apply_ksq(data, ksq);

	// Perform inverse FFT
	four1(data, nn, -1);

	// Copy data array back to phi
	for (i = 1; i <= N; i++) {
		phi[i] = data[2*(i-1)+1] * 1 / N; // Normalize the solution
	}

	// Print the solution values
	for (i = 1; i <= N; i++) {
		printf("%f ", phi[i] - phi[1]);
	}

	// Clean up
	free_vector(phi, 1, N);
	free_vector(ksq, 0, N-1);
	free_vector(data, 1, 2*N);

	return 0;
}
