function Xdot = Xu(t,in2,in3)
%XU
%    XDOT = XU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 21:52:05

th1 = in3(1,:);
th3 = in3(2,:);
th9 = in3(3,:);
th10 = in3(4,:);
th11 = in3(5,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th9-1.0;
t3 = 1.0./th9;
t4 = 1.0./th11;
t5 = 1.0./x5;
t6 = x3./4.0;
t9 = th9.*(2.3e+1./4.0e+2);
t10 = x4+1.0e-4;
t12 = t.*(1.97e+2./1.0e+3);
t14 = x2+1.034e-1;
t16 = th3.*9.002366863905325e-1;
t17 = th3.*7.20189349112426e+1;
t18 = x3.*1.344483583855441e-1;
t19 = x3+5.019e-1;
t7 = 1.0./t2;
t8 = t6+1.0;
t13 = exp(t12);
t15 = 1.0./t10;
t20 = -t17;
t21 = t18+1.0;
t22 = 1.0./t14;
t23 = 1.0./t19;
t11 = 1.0./t8;
t24 = 1.0./t21;
t25 = t23.*x3.*9.792;
t26 = -t25;
t28 = t4.*t11.*t15.*th3.*x4;
t30 = t22.*t24.*th1.*x2;
t31 = t11.*t15.*t16.*x4;
t32 = t11.*t15.*t17.*x4;
t33 = t11.*t15.*th3.*x4.*(-9.002366863905325e-1);
t27 = exp(t26);
t29 = -t28;
t34 = -t30;
t35 = t30.*8.0e+1;
t39 = t16+t33;
t40 = t20+t32;
t36 = -t35;
t38 = t9+t29;
t41 = exp(t40);
t43 = t23.*t27.*x3.*1.224e-1;
t37 = exp(t36);
t42 = t7.*t38.*8.0e+1;
t47 = t27+t41;
t51 = t39.*t41;
t44 = -t42;
t46 = t30.*t37;
t48 = 1.0./t47;
t53 = t43+t51;
t45 = exp(t44);
t49 = t37+t45;
t52 = t7.*t38.*t45;
t50 = 1.0./t49;
t54 = t46+t52;
t55 = t50.*t54;
t56 = t55.*8.0e+1;
t58 = t34+t55;
t57 = -t56;
t59 = t35+t57;
t60 = exp(t59);
t61 = t60+1.0;
t62 = 1.0./t61;
mt1 = [x1.*(t48.*t53.*(3.1e+1./2.0e+2)+th9.*(t55-2.3e+1./4.0e+2)-t3.*t5.*t13.*5.4175e-4+t60.*t62.*th10.*(t30-t55)),t34.*x1-t3.*t5.*t13.*(x2-4.5e+2).*5.4175e-4,-x1.*(t48.*t53-t60.*t62.*(t30-t55).*3.652e-1+t60.*t62.*th10.*(t30-t55).*3.652e-1)-t3.*t5.*t13.*x3.*5.4175e-4];
mt2 = [x4.*-1.8e+4-x1.*(t48.*t53.*6.427915e-1+t11.*t15.*th3.*x4)-t3.*t5.*t13.*x4.*1.0835e-3+6.3e+2,t3.*t13.*5.4175e-4];
Xdot = reshape([mt1,mt2],5,1);
