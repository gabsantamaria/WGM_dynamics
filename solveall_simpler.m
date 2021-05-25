%% solveall: function description
function [a1234_num, b012345_num] = solveall_simpler(r0_, r2_, r4_, r5_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_L3_, BET_coup_L4_)
%function [a1234_num, b1234_num, b0_num, b5_num] = solveall_simpler(r0_, r2_, r4_, r5_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_)
%BET_coup_L3_
%BET_coup_L4_

	a0 = 1; a5 = 0;

	t0_ = sqrt(1- abs(r0_)^2);
	t2_ = sqrt(1- abs(r2_)^2);

	syms r0; syms t0; 
	syms r2; syms t2; 
	syms r4; syms t4; 
	syms r5; syms t5;

	M = [[-conj(r0), 0, 0, 0]; [0, r2, t2, 0]; [0, t2, -conj(r2), 0];  [0, 0, 0, r4]];
	syms BET_res;
	%syms BET_coup; %[DEPRECATED]
		syms BET_coup_L3;
		syms BET_coup_L4;
	syms L1; syms L2;
	syms L3; syms L4;

	P = [[0, exp(-1i*BET_res*L2), 0, 0]; [exp(-1i*BET_res*L1), 0, 0, 0]; [0, 0, 0, exp(-1i*BET_coup_L4*L4)]; [0, 0, exp(-1i*BET_coup_L3*L3), 0] ];

	%a1234 = inv(eye(4) - P*M)*P*[1;0;0;0]*a0*t0;
    %a1234 = P*inv(eye(4) - M*P)*[1;0;0;0]*a0*t0;
    b1234 = inv(eye(4) - M*P)*[1;0;0;0]*a0*t0;
%a1234 = P*(   (eye(4) - M*P)\[1;0;0;0]   )*a0*t0;
	%b1234 = M*a1234;
	a1234 = P*b1234;
    %b1234 = inv(P)*a1234;
%b1234 = (P)\a1234;
	b0 = r0*a0       + t0*a1234(1);
	b5 = t5*a1234(4) + r5*a5;

	f_b0   (r0, r2, r4, r5, t0, t2, t4, t5, L1, L2, L3, L4, BET_res, BET_coup_L3, BET_coup_L4) = b0;
	f_b5   (r0, r2, r4, r5, t0, t2, t4, t5, L1, L2, L3, L4, BET_res, BET_coup_L3, BET_coup_L4) = b5;
	f_b1234(r0, r2, r4, r5, t0, t2, t4, t5, L1, L2, L3, L4, BET_res, BET_coup_L3, BET_coup_L4) = b1234;
	f_a1234(r0, r2, r4, r5, t0, t2, t4, t5, L1, L2, L3, L4, BET_res, BET_coup_L3, BET_coup_L4) = a1234;

	h_b0 = matlabFunction(f_b0);
	h_b5 = matlabFunction(f_b5);
	h_b1234 = matlabFunction(f_b1234);
	h_a1234 = matlabFunction(f_a1234);

	a1234_num = 				h_a1234 (r0_, r2_, r4_, r5_, t0_, t2_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_L3_, BET_coup_L4_);
	b012345_num = 				h_b0    (r0_, r2_, r4_, r5_, t0_, t2_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_L3_, BET_coup_L4_);
	b012345_num = [b012345_num;    h_b1234(r0_, r2_, r4_, r5_, t0_, t2_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_L3_, BET_coup_L4_)];
	b012345_num = [b012345_num;    h_b5(r0_, r2_, r4_, r5_, t0_, t2_, t4_, t5_, L1_, L2_, L3_, L4_, BET_res_, BET_coup_L3_, BET_coup_L4_)];


	%r0 = r0_;
	%r2 = r2_;
	%r4 = r4_;
	%r5 = r5_;
	%t0 = t0_;
	%t2 = t2_;
	%t4 = t4_;
	%t5 = t5_;
	%L1 = L1_;
	%L2 = L2_;
	%L3 = L3_;
	%L4 = L4_;
	%BET_res = BET_res_;
	%BET_coup = BET_coup_;
%
%
	%a1234_num = double(subs(a1234));
	%b1234_num = double(subs(b1234));
	%b0_num    = double(subs(b0));
	%b5_num    = double(subs(b5));

end