function Xdot = Xu(t,in2,in3)
%XU
%    XDOT = XU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 13:58:41

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
th13 = in3(13,:);
th14 = in3(14,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th4.*th9;
t3 = th5+x2;
t4 = th6+x3;
t5 = th9-1.0;
t6 = th13-1.0;
t7 = 1.0./th7;
t8 = 1.0./th8;
t9 = 1.0./th9;
t10 = 1.0./th11;
t11 = 1.0./x5;
t20 = t.*(1.1e+1./1.0e+2);
t22 = x4+1.0e-4;
t12 = t7.*x3;
t13 = t8.*x3;
t14 = 1.0./t3;
t15 = 1.0./t4;
t16 = 1.0./t5;
t17 = 1.0./t6;
t21 = exp(t20);
t25 = 1.0./t22;
t18 = t12+1.0;
t19 = t13+1.0;
t26 = t15.*th2.*x3.*1.0e+2;
t23 = 1.0./t18;
t24 = 1.0./t19;
t27 = -t26;
t28 = exp(t27);
t29 = t14.*t24.*th1.*x2;
t35 = t23.*t25.*th3.*x4;
t30 = -t29;
t31 = t15.*t28.*th2.*x3;
t32 = t29.*1.0e+2;
t36 = -t35;
t38 = t10.*t35;
t41 = t17.*th12.*(t35-th3).*-1.0e+2;
t33 = -t32;
t37 = t36+th3;
t39 = t10.*t36;
t42 = exp(t41);
t34 = exp(t33);
t40 = t2+t39;
t47 = t28+t42;
t51 = -t17.*t42.*th12.*(t35-th3);
t52 = t17.*t42.*th12.*(t35-th3);
t43 = t29.*t34;
t44 = t16.*t40.*1.0e+2;
t48 = 1.0./t47;
t54 = t31+t52;
t45 = -t44;
t56 = t48.*t54;
t46 = exp(t45);
t57 = t56.*th13;
t49 = t34+t46;
t53 = t16.*t40.*t46;
t58 = -t57;
t50 = 1.0./t49;
t55 = t43+t53;
t59 = t50.*t55;
t60 = t59.*1.0e+2;
t62 = t30+t59;
t61 = -t60;
t63 = t32+t61;
t64 = exp(t63);
t65 = t64+1.0;
t66 = 1.0./t65;
t67 = -t64.*t66.*th10.*(t29-t59);
Xdot = [x1.*(t57+th9.*(t59-th4)-t9.*t11.*t21.*3.025e-4+t64.*t66.*th10.*(t29-t59));t30.*x1-t9.*t11.*t21.*(x2-4.5e+2).*3.025e-4;-x1.*(t56-th14.*(t67+t64.*t66.*(t29-t59)))-t9.*t11.*t21.*x3.*3.025e-4;x4.*-1.8e+4-x1.*(t35+th12.*(t56+t58))-t9.*t11.*t21.*x4.*6.05e-4+1.44e+2;t9.*t21.*3.025e-4];