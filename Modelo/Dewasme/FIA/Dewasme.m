function Xdot = Dewasme(t,in2,in3)
%DEWASME
%    XDOT = DEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    06-Jul-2021 15:57:26

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th9 = in3(9,:);
th10 = in3(10,:);
th11 = in3(11,:);
th12 = in3(12,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th1+x2;
t3 = th5+x3;
t4 = th7+x3;
t5 = 1.0./th3;
t6 = 1.0./th4;
t7 = 1.0./th9;
t8 = 1.0./x5;
t12 = t.*(1.1e+1./1.0e+2);
t14 = x4+1.0e-4;
t9 = 1.0./t2;
t10 = 1.0./t3;
t11 = 1.0./t4;
t13 = exp(t12);
t16 = 1.0./t14;
t15 = t9.*th2.*x2;
t21 = t6.*t10.*t16.*th5.*th11.*x4;
t17 = t15.*1.0e+2;
t22 = -t21;
t23 = t21.*1.0e+2;
t18 = -t17;
t24 = -t23;
t26 = t15+t22;
t19 = exp(t18);
t25 = exp(t24);
t27 = t17+t24;
t33 = t7.*t11.*t26.*th4.*x3.*1.0e+2;
t20 = t15.*t19;
t28 = exp(t27);
t30 = t19+t25;
t34 = -t33;
t36 = t21.*t25;
t29 = t28+1.0;
t32 = 1.0./t30;
t35 = exp(t34);
t39 = t20+t36;
t31 = 1.0./t29;
t37 = t35+1.0;
t38 = 1.0./t37;
Xdot = [-x1.*(t5.*t8.*t13.*3.025e-4-t32.*t39.*th3-t26.*t28.*t31.*th8+t7.*t11.*t26.*t35.*t38.*th4.*th10.*x3);-x1.*(t32.*t39-t26.*t28.*t31)-t5.*t8.*t13.*(x2-4.5e+2).*3.025e-4;x1.*(t26.*t28.*t31.*th6+t7.*t11.*t26.*t35.*t38.*th4.*x3)-t5.*t8.*t13.*x3.*3.025e-4;x4.*-1.8e+4-x4.*(t32.*t39.*th4+t26.*t28.*t31.*th12-t11.*t26.*t35.*t38.*th4.*x3)-t5.*t8.*t13.*x4.*3.025e-4+6.3e+2;t5.*t13.*3.025e-4];
