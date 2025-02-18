#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//
// Modified numerical recipes transform
//

#define SWAP(a, b)      tempr = a; a = b; b = tempr

void four1(double* data, int n, int isign) {
	int nn, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
	nn = n << 1;
	j = 1;
	for (i = 1; i < nn; i += 2) {
		if (j > i) {
			SWAP(data[j - 1], data[i - 1]);
			SWAP(data[j], data[i]);
		}
		m = n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (nn > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= nn; i += istep) {
				j = i + mmax;
				tempr = wr * data[j - 1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j - 1];
				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

// 
// End NR code
//

#define N 16 // Size of the grid (must be a power of 2)
#define L 1.0 // Physical size of the domain

//
// Adapted from Brian's Python code
//

void initialize_rhs(double *rhs) {
    for (int i = 0; i < N; i++) {
		double x = ((double)i)/N;
		double pi_sq = M_PI * M_PI;
		rhs[2*i] = -16 * pi_sq * sin(4 * M_PI * x) - 8 * pi_sq * pow(sin(2 * M_PI * x), 2.0) + 8 * pi_sq * pow(cos(2 * M_PI * x), 2.0);
    }
}

void construct_ksq(double *ksq) {

	int i;
	double k;

	for (i = 0; i < N; i++) {
		
		k = (i <= N / 2) ? i : i - N;
		
		// Square and multiply with the square of the prefix 2*pi/L
		ksq[i] = -(4 * (M_PI * M_PI) / (L * L)) * ( k * k );

	}

} 

void apply_ksq(double *rhs, double *ksq) {
	
	int i;
	
	for (i = 0; i < N; i++) {
		
		// At 0,0 don't normalize
		if (ksq[i] != 0) {
			rhs[2*i] /= ksq[i];
			rhs[2*i+1] /= ksq[i]; 
		} 
	}

	// Zero out first component
	rhs[0] = 0.0;
	rhs[1] = 0.0;

}

int main() {

	int i, j;

	double *data;
	double *ksq;
	
	data = calloc(2*N, sizeof(double)); // Flattened data array, 2*N^2 to hold complex values, 1-indexed
	ksq = malloc(N * sizeof(double)); // Flattened data array, 2*N^2 to hold complex values, 1-indexed
	unsigned long nn = N; // Size of all dims (uses 1-indexing so first element is padding)
	
	double *phi = malloc(N * sizeof(double));   // Solution in real space

	// Initialize the data
	initialize_rhs(data);

	for (i = 0; i < 2*N; i+=2) {
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
	for (i = 0; i < N; i++) {
		phi[i] = data[2*i] * 1 / N; // Normalize the solution
	}

	// Print the solution values
	for (i = 0; i < N; i++) {
		printf("%f ", phi[i] - phi[0]);
	}

	// Clean up
	free(data);
	free(ksq);
	free(phi);

	return 0;
}
