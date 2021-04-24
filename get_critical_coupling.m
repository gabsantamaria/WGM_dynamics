function [r] = get_critical_coupling(Qi, mp)
	r = exp(-2*pi*mp/(2*Qi));
end