function [r4_, r5_, t4_, t5_] = get_interferometer(ra, ta, rb, tb, la, lb, BET_La, BET_Lb)
	nf = length(BET_La);
	M1 = repmat([[rb, tb]; [tb, -conj(rb)]], [1,1,nf]);
	zer = repmat(0,[1,1,nf]);
	expb = zeros(1,1,nf); expb(1,1,:) = exp(-1i*BET_Lb*lb);
	expa = zeros(1,1,nf); expa(1,1,:) = exp(-1i*BET_La*la);
	%M2 = [[exp(-1i*BET*lb), 0]; [0, exp(-1i*BET*la)]];
	M2 = [[expb, zer]; [zer, expa]];
	M3 = repmat([[ra, ta]; [ta, -conj(ra)]], [1,1,nf]);
	%CoupMat = M1*M2*M3;
	CoupMat = pagemtimes(M1, pagemtimes(M2, M3)); 
	r4_ = CoupMat(1,1,:); r4_ = transpose(r4_(:));
	t4_ = CoupMat(1,2,:); t4_ = transpose(t4_(:));
	t5_ = CoupMat(2,1,:); t5_ = transpose(t5_(:));
	r5_ = CoupMat(2,2,:); r5_ = transpose(r5_(:));
end