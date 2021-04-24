function [t] = get_lossless_t(r)
	t = sqrt(1- abs(r)^2);
end