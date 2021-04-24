close all; pause(0.1); c = 299792458;

lambda0 = 1550e-9; fspan = 2e12; 
%--------DESIGN PARAMETERS--------------------
neff_res = 1.8305; 			%effective index of resonator
neff_coup = 1.8305; 		%effective index of coupler waveguides
loss_res = 2.7; 			%loss of resonator (dB/m)
loss_coup = 2.7; 			%loss of coupler waveguides (dB/m)
R0 = 52.13160512687181e-6;	%radius of resonator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finding resonance frequency, asigning frequency vector and finding complex betas
[BET_res_, f, Qi, f0, mp, lda0]  = get_beta(lambda0, fspan, neff_res, loss_res, R0);
bet_pp_coup = loss_coup/(20*log10(exp(1))); BET_coup_ = (2*pi*f/c)*neff_coup -1i*bet_pp_coup;



%--------DESIGN PARAMETERS--------------------
L1_ = 2*pi*R0/2; 	%L1 in my drawing
L2_ = 2*pi*R0/2; 	%L2 in my drawing 
L3_ = 10e-6; 	 	%L3 in my drawing
L4_ = 10e-6;	 	%L4 in my drawing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%att_roundtrip is the single-roundtrip field attenuation factor (1 means no absorptio/radiation loss)
att_roundtrip = get_critical_coupling(Qi, mp);

%--------DESIGN PARAMETERS--------------------
r0_ = -(att_roundtrip^1.02);		%Input waveguide - resonator coupling strength (make it have 180º phase not to induce intra-cavity phase shift)
r2_ = 0+1*(att_roundtrip.^50);		%Resonator - coupler coupling strength
ra = sqrt(0.6);						%Coupling strength of first 50/50 splitter of the M-Z interferometer
rb = -sqrt(0.4); 					%Coupling strength of second 50/50 splitter of the M-Z interferometer
lb = 10e-6; 						%Length of "b" arm of the M-Z interferometer
la = lb + 0.5*c/((500e9)*neff_res);	%Length of "a" arm of the M-Z interferometer
ta = get_lossless_t(ra);			%Transmission coupling strength of first 50/50 splitter of the M-Z interferometer (use get_lossless_t(ra) for lossless coupler)
tb = get_lossless_t(rb);			%Transmission coupling strength of second 50/50 splitter of the M-Z interferometer (use get_lossless_t(rb) for lossless coupler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Building M-Z asymmetric interferometer
[r4_, r5_, t4_, t5_] = get_interferometer(ra, ta, rb, tb, la, lb, BET);

%Solving all modes symbolically and then numerically
[a1234_num, b012345_num] = solveall_simpler(r0_, r2_, r4_, r5_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%				Plotting routines				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; plot((f-f0)*1e-9, abs(r4_));  
ylabel('MZ amplitude response', 'interpreter', 'latex');
yyaxis right;
plot((f-f0)*1e-9, angle(r4_), '--')
xlabel('Laser detuning: $$\nu-\nu_0$$ [GHz]', 'interpreter', 'latex');
ylabel('MZ phase response', 'interpreter', 'latex');
legend('$$|r_4|$$', '$$\angle r_4$$', 'interpreter', 'latex')

figure; semilogy((f-f0)*1e-9, abs(b012345_num(1,:)), 'linewidth', 2); hold on;
%plot((f-f0)*1e-9, abs(b012345_num(2,:)));
semilogy((f-f0)*1e-9, abs(b012345_num(6,:)));
semilogy((f-f0)*1e-9, abs(b012345_num(2,:)), '--', 'linewidth', 2);
semilogy((f-f0)*1e-9, abs(a1234_num(3,:)));
xlabel('Laser detuning: $$\nu-\nu_0$$ [GHz]', 'interpreter', 'latex'); ylabel('Magnitude of the mode [$$\sqrt{\mathrm{W}}$$]', 'interpreter', 'latex');
legend('b0', 'b5', 'b1', 'a3');
%plot((f-f0)*1e-9, abs(b012345_num(6,:)));
[Qloaded, Dnu_loaded, Wstored, Ploss_abs, Ploss_leak, P_leak0, P_leak5] = calculate_Ql(a1234_num, b012345_num, L1_, L2_, L3_, L4_, neff_res, neff_coup, r0_, f);
figure;
semilogy((f-f0)*1e-9,Qloaded); ylabel('Loaded Q', 'interpreter', 'latex');  
%yl = ylim; ylim([yl(1), 1.1*Qi]); hold on;
%yline(Qi/2,'b--'); yline(Qi,'b-.'); hold off;
yyaxis right;
semilogy((f-f0)*1e-9,Dnu_loaded*1e-6); ylabel('$$\Delta\nu_\mathrm{loaded}$$ [MHz]', 'interpreter', 'latex');
xlabel('Laser detuning: $$\nu-\nu_0$$ [GHz]', 'interpreter', 'latex');

legend('Loaded Q', 'Loaded bandwidth', 'critical coupling loaded Q', 'intrinsic Q');


%figure; semilogy((f-f0)*1e-9,Wstored, 'linewidth', 4); 
%yyaxis right;  semilogy((f-f0)*1e-9,Ploss_abs,'linewidth', 2);  hold on;  semilogy((f-f0)*1e-9,Ploss_leak, 'k--');
%legend('Wstored', 'Ploss abs', 'Ploss leak')

figure; semilogy((f-f0)*1e-9,(2*pi*f.*Wstored)./Ploss_abs, 'linewidth', 4); hold on;
semilogy((f-f0)*1e-9,(2*pi*f.*Wstored)./P_leak0, 'linewidth', 2);
semilogy((f-f0)*1e-9,(2*pi*f.*Wstored)./P_leak5, 'linewidth', 2);
legend('Qi', 'Qc0', 'Qc5')
xlabel('Laser detuning: $$\nu-\nu_0$$ [GHz]', 'interpreter', 'latex');
ylabel('Intrinsic/Coupling Q', 'interpreter', 'latex');  
%energy conservation
%(33.4744^2)*(1-get_critical_coupling(Qi, mp)^2)+0.39029^2 + 0.00880521^2