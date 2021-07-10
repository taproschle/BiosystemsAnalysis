function Xdot = Xu(t,in2,in3)
%XU
%    XDOT = XU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 19:51:51

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th7 = in3(4,:);
th9 = in3(5,:);
th10 = in3(6,:);
th11 = in3(7,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th9-1.0;
t3 = 1.0./th7;
t4 = 1.0./th9;
t5 = 1.0./th11;
t6 = 1.0./x5;
t10 = t.*(1.1e+1./1.0e+2);
t11 = th9.*(2.3e+1./4.0e+2);
t13 = x4+1.0e-4;
t15 = th3.*9.002366863905325e+1;
t16 = x2+1.034e-1;
t19 = th3.*9.002366863905325e-1;
t20 = x3.*1.344483583855441e-1;
t21 = x3+5.019e-1;
t7 = t3.*x3;
t8 = 1.0./t2;
t12 = exp(t10);
t17 = 1.0./t13;
t18 = -t15;
t22 = t20+1.0;
t23 = 1.0./t16;
t24 = 1.0./t21;
t9 = t7+1.0;
t25 = 1.0./t22;
t26 = t24.*th2.*x3.*1.0e+2;
t14 = 1.0./t9;
t27 = -t26;
t31 = t23.*t25.*th1.*x2;
t28 = exp(t27);
t29 = t5.*t14.*t17.*th3.*x4;
t32 = t14.*t15.*t17.*x4;
t33 = -t31;
t34 = t14.*t17.*t19.*x4;
t35 = t31.*1.0e+2;
t36 = t14.*t17.*th3.*x4.*(-9.002366863905325e-1);
t30 = -t29;
t37 = -t35;
t39 = t24.*t28.*th2.*x3;
t41 = t18+t32;
t43 = t19+t36;
t38 = exp(t37);
t40 = t11+t30;
t42 = exp(t41);
t44 = t8.*t40.*1.0e+2;
t47 = t31.*t38;
t48 = t28+t42;
t52 = t42.*t43;
t45 = -t44;
t49 = 1.0./t48;
t54 = t39+t52;
t46 = exp(t45);
t50 = t38+t46;
t53 = t8.*t40.*t46;
t51 = 1.0./t50;
t55 = t47+t53;
t56 = t51.*t55;
t57 = t56.*1.0e+2;
t59 = t33+t56;
t58 = -t57;
t60 = t35+t58;
t61 = exp(t60);
t62 = t61+1.0;
t63 = 1.0./t62;
mt1 = [x1.*(t49.*t54.*(3.1e+1./2.0e+2)+th9.*(t56-2.3e+1./4.0e+2)-t4.*t6.*t12.*3.025e-4+t61.*t63.*th10.*(t31-t56)),t33.*x1-t4.*t6.*t12.*(x2-4.5e+2).*3.025e-4,-x1.*(t49.*t54-t61.*t63.*(t31-t56).*3.652e-1+t61.*t63.*th10.*(t31-t56).*3.652e-1)-t4.*t6.*t12.*x3.*3.025e-4,x4.*-1.8e+4-x1.*(t49.*t54.*6.427915e-1+t14.*t17.*th3.*x4)-t4.*t6.*t12.*x4.*6.05e-4+6.3e+2];
mt2 = [t4.*t12.*3.025e-4];
Xdot = reshape([mt1,mt2],5,1);