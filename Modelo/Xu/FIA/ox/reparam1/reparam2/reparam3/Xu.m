function Xdot = Xu(t,in2,in3)
%XU
%    XDOT = XU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 16:26:07

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th7 = in3(5,:);
th8 = in3(6,:);
th9 = in3(7,:);
th10 = in3(8,:);
th11 = in3(9,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th5+x2;
t3 = th9-1.0;
t4 = 1.0./th7;
t5 = 1.0./th8;
t6 = 1.0./th9;
t7 = 1.0./th11;
t8 = 1.0./x5;
t15 = t.*(1.1e+1./1.0e+2);
t16 = th9.*(2.3e+1./4.0e+2);
t18 = x4+1.0e-4;
t21 = th3.*9.002366863905325e+1;
t24 = th3.*9.002366863905325e-1;
t25 = x3+5.019e-1;
t9 = t4.*x3;
t10 = t5.*x3;
t11 = 1.0./t2;
t12 = 1.0./t3;
t17 = exp(t15);
t22 = 1.0./t18;
t23 = -t21;
t26 = 1.0./t25;
t13 = t9+1.0;
t14 = t10+1.0;
t28 = t26.*th2.*x3.*1.0e+2;
t19 = 1.0./t13;
t20 = 1.0./t14;
t30 = -t28;
t27 = t11.*t20.*th1.*x2;
t32 = exp(t30);
t35 = t7.*t19.*t22.*th3.*x4;
t37 = t19.*t21.*t22.*x4;
t38 = t19.*t22.*t24.*x4;
t39 = t19.*t22.*th3.*x4.*(-9.002366863905325e-1);
t29 = -t27;
t31 = t27.*1.0e+2;
t36 = -t35;
t40 = t26.*t32.*th2.*x3;
t43 = t23+t37;
t45 = t24+t39;
t33 = -t31;
t41 = t16+t36;
t44 = exp(t43);
t34 = exp(t33);
t46 = t12.*t41.*1.0e+2;
t49 = t32+t44;
t53 = t44.*t45;
t42 = t27.*t34;
t47 = -t46;
t50 = 1.0./t49;
t55 = t40+t53;
t48 = exp(t47);
t51 = t34+t48;
t54 = t12.*t41.*t48;
t52 = 1.0./t51;
t56 = t42+t54;
t57 = t52.*t56;
t58 = t57.*1.0e+2;
t60 = t29+t57;
t59 = -t58;
t61 = t31+t59;
t62 = exp(t61);
t63 = t62+1.0;
t64 = 1.0./t63;
mt1 = [x1.*(t50.*t55.*(3.1e+1./2.0e+2)+th9.*(t57-2.3e+1./4.0e+2)-t6.*t8.*t17.*3.025e-4+t62.*t64.*th10.*(t27-t57)),t29.*x1-t6.*t8.*t17.*(x2-4.5e+2).*3.025e-4,-x1.*(t50.*t55-t62.*t64.*(t27-t57).*3.652e-1+t62.*t64.*th10.*(t27-t57).*3.652e-1)-t6.*t8.*t17.*x3.*3.025e-4,x4.*-1.8e+4-x1.*(t50.*t55.*6.427915e-1+t19.*t22.*th3.*x4)-t6.*t8.*t17.*x4.*6.05e-4+6.3e+2];
mt2 = [t6.*t17.*3.025e-4];
Xdot = reshape([mt1,mt2],5,1);
