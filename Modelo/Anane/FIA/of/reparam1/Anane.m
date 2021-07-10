function Xdot = Anane(t,in2,in3)
%ANANE
%    XDOT = ANANE(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 20:45:22

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th6 = in3(5,:);
th7 = in3(6,:);
th8 = in3(7,:);
th9 = in3(8,:);
th10 = in3(9,:);
th11 = in3(10,:);
th12 = in3(11,:);
th13 = in3(12,:);
th14 = in3(13,:);
th15 = in3(14,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th5+x2;
t3 = th6+x3;
t4 = 1.0./th8;
t5 = 1.0./th9;
t6 = 1.0./th10;
t7 = 1.0./x5;
t12 = x4+1.0e-4;
t13 = t.*(1.97e+2./1.0e+3);
t8 = t4.*x3;
t9 = 1.0./t2;
t10 = 1.0./t3;
t14 = exp(t13);
t16 = 1.0./t12;
t11 = t8+1.0;
t15 = 1.0./t11;
t17 = t9.*t15.*th1.*x2;
t18 = t17+th7;
t19 = t5.*t17;
t20 = 1.0./t18;
t21 = t19+1.0;
t22 = 1.0./t21;
t23 = t17.*t20.*th3;
t24 = -t23;
t25 = t17+t24;
t26 = t16.*t25.*x4;
t27 = t26-6.7e-3;
Xdot = [x1.*(t23.*th11+t27.*th10-t6.*t7.*t14.*5.4175e-4+t10.*t22.*th2.*th14.*x3);-t17.*x1-t6.*t7.*t14.*(x2-4.5e+2).*5.4175e-4;x1.*(t23.*th15-t10.*t22.*th2.*x3)-t6.*t7.*t14.*x3.*5.4175e-4;x4.*-1.8e+4-x1.*(t27.*th12+t10.*t22.*th2.*th13.*x3)-t6.*t7.*t14.*x4.*1.0835e-3+6.3e+2;t6.*t14.*5.4175e-4];
