ImportAll(realdft);
ImportAll(filtering);
ImportAll(dct_dst);

n := 4;
l := 1;

#f := Mat([ [1..4], [5..8], [9..12], [13..16] ]);
f := Mat([
	[0, 0.0005, 0.0086, 0.0005],
	[0.0005, 2.027, 18.6732, 2.027],
	[0.0086, 18.6732, -200, 18.6732],
	[0.0005, 2.027, 18.6732, 2.027]
]);

factor := 2*d_PI/l;

k := factor * ([0..n/2]::-Reversed([1..(n/2-1)]));

kx := Mat(Replicate(n, k));
ky := kx.transpose();

# FIXME: Generalize for any n
HProd_4 := function(mat_a, mat_b)
	local a, b;
	a := HStack(Diag(MatSPL(mat_a) * [ 1, 0, 0, 0 ]), 
		Diag(MatSPL(mat_a) * [ 0, 1, 0, 0 ]), 
		Diag(MatSPL(mat_a) * [ 0, 0, 1, 0 ]), 
		Diag(MatSPL(mat_a) * [ 0, 0, 0, 1 ])
		);
	b := VStack((mat_b * Mat([ [ 1, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ] ])), 
		(mat_b * Mat([ [ 0, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ] ])),
		(mat_b * Mat([ [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 0 ] ])),
		(mat_b * Mat([ [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 0 ], [ 0, 0, 0, 1 ] ]))
		);
	return a * b;
end;

# TODO: Need better Hadamard Product
delsq := -(HProd_4(kx, kx) + HProd_4(ky, ky));
delsq_inv := Mat(List(MatSPL(delsq), r->List(r, v->1/v)));

# FIXME: Why does the inf mess everything up? (It turns the entire results matrix into nan's)
# FIXME: How to properly check for -inf?
delsq_inv_fix := Mat(List(MatSPL(delsq_inv), r->List(r, v->When(v < -9999999, 1, v))));

fhat := HProd_4(DFT(n) * f * DFT(n), delsq_inv_fix);

# TODO: This is a result of the earlier mess with inf in delsq
fhat_fix := Mat(List(MatSPL(fhat), r->List(r, v->When(v < -50, 0, v))));
phi := (1/n^2) * ( IDFT(n) * fhat_fix * IDFT(n) );

phi_real := Mat(List(MatSPL(phi), r->List(r, v->Value(TReal, Re(v)))));
phi_real_zero_mean := MatSPL(phi_real) - MatSPL(phi_real)[1][1];

pm(phi_real_zero_mean);

