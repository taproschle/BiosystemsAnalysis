function Xdot = Xu(t,in2,in3)
%XU
%    XDOT = XU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 21:53:21

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
t9 = t.*(1.1e+1./1.0e+2);
t10 = th9.*(2.3e+1./4.0e+2);
t12 = x4+1.0e-4;
t14 = th3.*9.002366863905325e+1;
t15 = x2+1.034e-1;
t18 = th3.*9.002366863905325e-1;
t19 = x3.*1.344483583855441e-1;
t20 = x3+5.019e-1;
t7 = 1.0./t2;
t8 = t6+1.0;
t11 = exp(t9);
t16 = 1.0./t12;
t17 = -t14;
t21 = t19+1.0;
t22 = 1.0./t15;
t23 = 1.0./t20;
t13 = 1.0./t8;
t24 = 1.0./t21;
t25 = t23.*x3.*(3.06e+2./2.5e+1);
t26 = -t25;
t28 = t4.*t13.*t16.*th3.*x4;
t30 = t13.*t14.*t16.*x4;
t31 = t22.*t24.*th1.*x2;
t32 = t13.*t16.*t18.*x4;
t33 = t13.*t16.*th3.*x4.*(-9.002366863905325e-1);
t27 = exp(t26);
t29 = -t28;
t34 = -t31;
t35 = t31.*1.0e+2;
t39 = t17+t30;
t41 = t18+t33;
t36 = -t35;
t38 = t10+t29;
t40 = exp(t39);
t42 = t23.*t27.*x3.*1.224e-1;
t37 = exp(t36);
t43 = t7.*t38.*1.0e+2;
t47 = t27+t40;
t51 = t40.*t41;
t44 = -t43;
t46 = t31.*t37;
t48 = 1.0./t47;
t53 = t42+t51;
t45 = exp(t44);
t49 = t37+t45;
t52 = t7.*t38.*t45;
t50 = 1.0./t49;
t54 = t46+t52;
t55 = t50.*t54;
t56 = t55.*1.0e+2;
t58 = t34+t55;
t57 = -t56;
t59 = t35+t57;
t60 = exp(t59);
t61 = t60+1.0;
t62 = 1.0./t61;
mt1 = [x1.*(t48.*t53.*(3.1e+1./2.0e+2)+th9.*(t55-2.3e+1./4.0e+2)-t3.*t5.*t11.*3.025e-4+t60.*t62.*th10.*(t31-t55)),t34.*x1-t3.*t5.*t11.*(x2-4.5e+2).*3.025e-4,-x1.*(t48.*t53-t60.*t62.*(t31-t55).*3.652e-1+t60.*t62.*th10.*(t31-t55).*3.652e-1)-t3.*t5.*t11.*x3.*3.025e-4,x4.*-1.8e+4-x1.*(t48.*t53.*6.427915e-1+t13.*t16.*th3.*x4)-t3.*t5.*t11.*x4.*6.05e-4+6.3e+2];
mt2 = [t3.*t11.*3.025e-4];
Xdot = reshape([mt1,mt2],5,1);
