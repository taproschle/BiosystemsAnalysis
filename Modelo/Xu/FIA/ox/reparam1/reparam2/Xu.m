function Xdot = Xu(t,in2,in3)
%XU
%    XDOT = XU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 16:23:44

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
th14 = in3(11,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th5+x2;
t3 = th6+x3;
t4 = th9-1.0;
t5 = 1.0./th7;
t6 = 1.0./th8;
t7 = 1.0./th9;
t8 = 1.0./th11;
t9 = 1.0./x5;
t17 = t.*(1.1e+1./1.0e+2);
t18 = th9.*(2.3e+1./4.0e+2);
t20 = x4+1.0e-4;
t23 = th3.*9.002366863905325e+1;
t26 = th3.*9.002366863905325e-1;
t10 = t5.*x3;
t11 = t6.*x3;
t12 = 1.0./t2;
t13 = 1.0./t3;
t14 = 1.0./t4;
t19 = exp(t17);
t24 = 1.0./t20;
t25 = -t23;
t15 = t10+1.0;
t16 = t11+1.0;
t27 = t13.*th2.*x3.*1.0e+2;
t21 = 1.0./t15;
t22 = 1.0./t16;
t28 = -t27;
t29 = exp(t28);
t30 = t12.*t22.*th1.*x2;
t36 = t8.*t21.*t24.*th3.*x4;
t38 = t21.*t23.*t24.*x4;
t39 = t21.*t24.*t26.*x4;
t40 = t21.*t24.*th3.*x4.*(-9.002366863905325e-1);
t31 = -t30;
t32 = t13.*t29.*th2.*x3;
t33 = t30.*1.0e+2;
t37 = -t36;
t43 = t25+t38;
t45 = t26+t40;
t34 = -t33;
t41 = t18+t37;
t44 = exp(t43);
t35 = exp(t34);
t46 = t14.*t41.*1.0e+2;
t49 = t29+t44;
t53 = t44.*t45;
t42 = t30.*t35;
t47 = -t46;
t50 = 1.0./t49;
t55 = t32+t53;
t48 = exp(t47);
t51 = t35+t48;
t54 = t14.*t41.*t48;
t52 = 1.0./t51;
t56 = t42+t54;
t57 = t52.*t56;
t58 = t57.*1.0e+2;
t60 = t31+t57;
t59 = -t58;
t61 = t33+t59;
t62 = exp(t61);
t63 = t62+1.0;
t64 = 1.0./t63;
t65 = -t62.*t64.*th10.*(t30-t57);
t66 = t62.*t64.*th10.*(t30-t57);
Xdot = [x1.*(t66+t50.*t55.*(3.1e+1./2.0e+2)+th9.*(t57-2.3e+1./4.0e+2)-t7.*t9.*t19.*3.025e-4);t31.*x1-t7.*t9.*t19.*(x2-4.5e+2).*3.025e-4;-x1.*(t50.*t55-th14.*(t65+t62.*t64.*(t30-t57)))-t7.*t9.*t19.*x3.*3.025e-4;x4.*-1.8e+4-x1.*(t50.*t55.*6.427915e-1+t21.*t24.*th3.*x4)-t7.*t9.*t19.*x4.*6.05e-4+6.3e+2;t7.*t19.*3.025e-4];
