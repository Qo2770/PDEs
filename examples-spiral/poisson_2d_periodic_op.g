ImportAll(realdft);
ImportAll(filtering);
ImportAll(dct_dst);

# === PREAMBLE ===

Class(_FList, FList, rec(tolist := self >> self.list));
SkewCirculant := l -> Toeplitz(Reversed(l)::(-DropLast(Reversed(l), 1)));
spiral.code.FList.tolist := self >> List(self.list, V);
_mat := m -> When(IsSPL(m), MatSPL(m), m);
AssertEqualMat := (a, b) -> let(am := _mat(a), bm := _mat(b), When(am=bm, true, Error(InfinityNormMat(am-bm))));
AssertEqualVec := (a, b) -> When(a=b, true, Error(InfinityNormMat([a-b])));
RCDiagTaps := (t, taps) -> RCDiag(_FList(TReal, MatSPL(t) * taps));
RotDiag := lst -> RCDiag(_FList(TReal, List(Zip2(lst, Flat(Replicate(Length(lst)/2, [1, -1]))), Product)));
_Diag := l -> Diag(_FList(TReal, l));

# ================

n := 8;
l := 1;

#f := Mat([ [1..4], [5..8], [9..12], [13..16] ]);
#f := Mat([
#	[0, 0.0005, 0.0086, 0.0005],
#	[0.0005, 2.027, 18.6732, 2.027],
#	[0.0086, 18.6732, -200, 18.6732],
#	[0.0005, 2.027, 18.6732, 2.027]
#]);
f := Mat([ 
	[ 6.666213e-08, 1.220733e-05, 4.789328e-04, 4.190820e-03, 8.571302e-03, 4.190820e-03, 4.789328e-04, 1.220733e-05 ], 
	[ 1.220733e-05, 2.040752e-03, 7.111223e-02, 5.513290e-01, 1.066115e+00, 5.513290e-01, 7.111223e-02, 2.040752e-03 ], 
	[ 4.789328e-04, 7.111223e-02, 2.026977e+00, 1.169231e+01, 1.867320e+01, 1.169231e+01, 2.026977e+00, 7.111223e-02 ], 
	[ 4.190820e-03, 5.513290e-01, 1.169231e+01, 2.358128e+01, -2.003021e+01, 2.358128e+01, 1.169231e+01, 5.513290e-01 ], 
	[ 8.571302e-03, 1.066115e+00, 1.867320e+01, -2.003021e+01, -2.000000e+02, -2.003021e+01, 1.867320e+01, 1.066115e+00 ], 
	[ 4.190820e-03, 5.513290e-01, 1.169231e+01, 2.358128e+01, -2.003021e+01, 2.358128e+01, 1.169231e+01, 5.513290e-01 ], 
	[ 4.789328e-04, 7.111223e-02, 2.026977e+00, 1.169231e+01, 1.867320e+01, 1.169231e+01, 2.026977e+00, 7.111223e-02 ], 
	[ 1.220733e-05, 2.040752e-03, 7.111223e-02, 5.513290e-01, 1.066115e+00, 5.513290e-01, 7.111223e-02, 2.040752e-03 ] 
]);

2d_dft := MatSPL(MDDFT([n, n], 1));
2d_idft := (1/n^2) * MatSPL(MDDFT([n, n], -1));

factor := 2*d_PI/l;

k := factor * ([0..n/2]::-Reversed([1..(n/2-1)]));

kx := Mat(Replicate(n, k));
ky := kx.transpose();

f_flat := Flat(f.element);
kx_flat := Flat(kx.element);
ky_flat := Flat(ky.element);

kx_sq := List(kx_flat, v->v*v);
ky_sq := List(ky_flat, v->v*v);

# Setup frequency constants
delsq := -(kx_sq + ky_sq);
delsq_inv := List(delsq, v->1/v);

# NOTE: Work-around since inf does not work properly in SPIRAL
# (dividing by it gives nan instead of 0)
delsq_inv[1] := 0;

delsq_inv_mat := MatSPL(_Diag(delsq_inv));

# Construct the 2D Poisson operator and apply to f
# FFT the input, multiply with delsq elements, then IFFT
op := 2d_idft * delsq_inv_mat * 2d_dft;
phi := op * f_flat;

# Convert to reals
phi_real := List(phi, v->Value(TReal, Re(v)));

# Zero-mean the result by convention
phi_real_zero_mean := phi_real - phi_real[1];

