function Xdot = Dewasme(t,in2,in3)
%DEWASME
%    XDOT = DEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 19:00:27

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th6 = in3(5,:);
th8 = in3(6,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th1+x2;
t3 = th5+x3;
t4 = 1.0./th3;
t5 = 1.0./x5;
t8 = x4+1.0e-4;
t10 = t.*(1.97e+2./1.0e+3);
t11 = x3+5.1e-3;
t6 = 1.0./t2;
t7 = 1.0./t3;
t12 = exp(t10);
t13 = 1.0./t8;
t15 = 1.0./t11;
t9 = t6.*th2.*x2;
t20 = t7.*t13.*th5.*x4.*1.017161140197789;
t22 = t7.*t13.*th5.*x4.*1.75025025025025e-1;
t24 = t7.*t13.*th5.*x4.*8.137289121582315e+1;
t14 = t9.*8.0e+1;
t18 = t9.*1.720720720720721e-1;
t21 = -t20;
t23 = -t22;
t25 = -t24;
t16 = -t14;
t26 = exp(t25);
t27 = t9+t21;
t28 = t14+t25;
t32 = t18+t23;
t17 = exp(t16);
t29 = exp(t28);
t35 = t20.*t26;
t36 = t15.*t32.*x3.*8.0e+1;
t19 = t9.*t17;
t30 = t29+1.0;
t31 = t17+t26;
t37 = -t36;
t33 = 1.0./t30;
t34 = 1.0./t31;
t38 = exp(t37);
t41 = t19+t35;
t39 = t38+1.0;
t40 = 1.0./t39;
mt1 = [-x1.*(t4.*t5.*t12.*5.4175e-4-t34.*t41.*th3-t27.*t29.*t33.*th8+t15.*t32.*t38.*t40.*x3.*1.7531),-x1.*(t34.*t41-t27.*t29.*t33)-t4.*t5.*t12.*(x2-4.5e+2).*5.4175e-4,x1.*(t27.*t29.*t33.*th6+t15.*t32.*t38.*t40.*x3)-t4.*t5.*t12.*x3.*5.4175e-4,x4.*-1.8e+4-x4.*(t34.*t41.*3.438e-1+(t27.*t29.*t33)./1.0e+4-t15.*t32.*t38.*t40.*x3.*(9.99e+2./5.0e+2))-t4.*t5.*t12.*x4.*5.4175e-4+6.3e+2];
mt2 = [t4.*t12.*5.4175e-4];
Xdot = reshape([mt1,mt2],5,1);
