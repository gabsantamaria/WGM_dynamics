function [BET, f, Qi, f0, mp, lda0] = get_beta(lambda0, fspan, neff, loss, R0)
	%loss in dB/m
	c = 299792458;

	mp = round(2*pi*R0/(lambda0/neff));
	f0 = c/(2*pi*R0*neff/mp);
	lda0 = c/f0;

	bet_pp = loss/(20*log10(exp(1)));
	%f0 = c/lambda0;
	Qi = pi*neff/(bet_pp*lambda0);
	df = 0.05*f0/Qi;
	f = f0 + [-fspan/2:df:fspan/2];
	w = 2*pi*f;
	BET = (w/c)*neff -1i*bet_pp;
end