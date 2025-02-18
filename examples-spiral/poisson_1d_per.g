ImportAll(realdft);
ImportAll(filtering);
ImportAll(dct_dst);
ImportAll(fftx.nonterms);

Class(_FList, FList, rec(tolist := self >> self.list));
_Diag := l -> Diag(_FList(TReal, l));
_mat := m -> When(IsSPL(m), MatSPL(m), m);
AssertEqualMat := (a, b) -> let(am := _mat(a), bm := _mat(b), When(am=bm, true, Error(InfinityNormMat(am-bm))));

# === === ===

# Setup args
n := 8;
l := 1;

# Data
f := [ 78.956835209999994, -157.91367041999999, -78.956835209999994, 
	157.91367041999999, 78.956835209999994, -157.91367041999999, 
	-78.956835209999994, 157.91367041999999 ];

# Poisson factor
factor := (2*d_PI)/l;

# Setup DFTs
dft := DFT(n, -1);
idft := (1/n) * DFT(n, 1);

# Zero-mean projection operators
ones_large := Mat(Replicate(n, Replicate(n, 1)));
ones_small := Mat(Replicate(n-1, Replicate(n-1, 1)));
bt_op_iso := I(n-1) + (-1/n) * ones_small;
bt_inv_iso := ones_small + I(n-1);
bt_op := HStack(O(n, 1), VStack(Mat([Replicate(n-1, -1/n)]), bt_op_iso));
bt_inv := DirectSum(O(1,1), bt_inv_iso);

# Poisson operators
k := ([1..n/2]::-Reversed([1..(n/2-1)]));
ksq := List(k, v->-v*v);
ksq_inv := List(ksq, v->1/v);
ones_large := Mat(Replicate(n, Replicate(n, 1)));

ksq_diag := DirectSum(O(1,1), _Diag(ksq));
ksq_diag_inv := DirectSum(O(1,1), _Diag(ksq_inv));

# Construct operator
zero_mean_op := bt_inv * ksq_diag_inv;
zero_mean_inv := ksq_diag * bt_op;

AssertEqualMat(DirectSum(O(1,1), I(n-1)), zero_mean_op * zero_mean_inv);

op_1d := dft * ksq_diag_inv * idft;
op_1d_inv := idft * ksq_diag * dft;
op := bt_op * dft * zero_mean_op * idft;

# Construct inverse filter
filt := Filt(n, [-1,2,-1]);

sct := Scat(fAdd(n+2, n, 1));
wrap := (Scat(fBase(n+2,0)) * Gath(fBase(n,n-1)) + Scat(fBase(n+2,n+1)) * Gath(fBase(n,0)));
ext := sct+wrap;

filtperi := filt * ext;
fpd := DFT(n, -1) * filtperi * DFT(n, 1);
fpdm := MatSPL(fpd);
fptaps := List([1..n], i->fpdm[i][i]);
fptapsi := List(fptaps, i->When(i=0,0,1/i));

ifilt := DFT(n, 1) * _Diag(fptapsi) * DFT(n, -1);

# Construct base changes
op1d_ifilt := MatSPL(HStack(op_1d, ifilt)); # TriangulizeMat is in-place and requires MatSPL
ifilt_op1d := MatSPL(HStack(ifilt, op_1d)); # TriangulizeMat is in-place and requires MatSPL

TriangulizeMat(op1d_ifilt);
TriangulizeMat(ifilt_op1d);

cb_op1d_ifilt := op1d_ifilt{[1..n]}{[(n+1)..(2*n)]};
cb_ifilt_op1d := ifilt_op1d{[1..n]}{[(n+1)..(2*n)]};

# Confirm 
AssertEqualMat(MatSPL(op_1d) * cb_op1d_ifilt, ifilt);
AssertEqualMat(MatSPL(ifilt) * cb_ifilt_op1d, MatSPL(op_1d));

pm(ifilt);
pm(op_1d);

# Apply operator to input
#res_raw := f * MatSPL(op);

#res_raw / (factor * factor);

#IdSubtractMean := I(n) + Mat(- 1/n * Replicate(n, Replicate(n, 1)));

#AssertEqualMat(IdSubtractMean, op_1d * op_1d_inv);
#AssertEqualMat(IdSubtractMean, bt * bt * op_1d * op_1d_inv * bt * bt);
#AssertEqualMat(IdSubtractMean, bt * bt * op_1d * bt * bt * op_1d_inv * bt * bt);
#AssertEqualMat(IdSubtractMean, bt * bt * op_1d_inv  * op_1d * bt * bt);


