function [Qloaded, Dnu_loaded, Wstored, Ploss_abs, Ploss_leak, P_leak0, P_leak5] = calculate_Ql(a1234_num, b012345_num, L1_, L2_, L3_, L4_, neff_res, neff_coup, r0_, f)
	c = 299792458;
	t0_ = sqrt(1-abs(r0_^2));
	Ploss_abs = (abs(b012345_num(2,:)).^2 - abs(a1234_num(2,:)).^2) + (abs(b012345_num(3,:)).^2 - abs(a1234_num(1,:)).^2) + (abs(b012345_num(5,:)).^2 - abs(a1234_num(3,:)).^2) + (abs(b012345_num(4,:)).^2 - abs(a1234_num(4,:)).^2);
	tau1 = L1_/(c/neff_res); tau2 = L2_/(c/neff_res); tau3 = L3_/(c/neff_coup); tau4 = L4_/(c/neff_coup); 
	Wstored = (abs(a1234_num(1,:)).^2)*tau2 + (abs(a1234_num(2,:)).^2)*tau1 + (abs(a1234_num(4,:)).^2)*tau3 + (abs(a1234_num(3,:)).^2)*tau4; 

	P_leak0 = (abs(a1234_num(1,:)).^2)*(abs(t0_).^2);
	P_leak5 = abs(b012345_num(6,:).^2);
	Ploss_leak =  P_leak0 + P_leak5;

	Qloaded = 2*pi*f.*Wstored./(Ploss_abs+Ploss_leak);
	Dnu_loaded = f./Qloaded;
end