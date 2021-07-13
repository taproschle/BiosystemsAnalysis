function dfdthSym = dfdthDewasme(t,in2,in3)
%DFDTHDEWASME
%    DFDTHSYM = DFDTHDEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    10-Jul-2021 19:56:35

th2 = in3(1,:);
th3 = in3(2,:);
th7 = in3(3,:);
th10 = in3(4,:);
th11 = in3(5,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th7+x3;
t3 = x2.^2;
t4 = x3.^2;
t5 = x4.^2;
t6 = 1.0./th3.^2;
t7 = 1.0./x5;
t11 = t.*(1.1e+1./1.0e+2);
t13 = x4+1.0e-4;
t14 = x2+1.414e-1;
t17 = x3+4.7323;
t8 = 1.0./t2;
t12 = exp(t11);
t15 = 1.0./t13;
t18 = 1.0./t14;
t20 = 1.0./t17;
t9 = t8.^2;
t10 = t8.^3;
t16 = t15.^2;
t19 = t18.^2;
t21 = t20.^2;
t22 = t18.*th2.*x2;
t30 = t15.*t20.*th11.*x4.*2.368518518518519;
t32 = t15.*t20.*th11.*x4.*1.376468877254218e+1;
t34 = t15.*t20.*th11.*x4.*1.376468877254218e+3;
t35 = t15.*t20.*th11.*x4.*2.752937754508435e+3;
t23 = t22.*1.0e+2;
t24 = t22.*2.0e+2;
t27 = t22.*1.720720720720721e-1;
t31 = -t30;
t33 = -t32;
t36 = -t34;
t37 = -t35;
t25 = -t23;
t39 = exp(t36);
t41 = t22+t33;
t42 = t23+t36;
t43 = t24+t37;
t46 = t27+t31;
t26 = exp(t25);
t44 = exp(t42);
t45 = exp(t43);
t49 = t46.^2;
t54 = t8.*t46.*x3.*1.0e+2;
t55 = t8.*t46.*x3.*2.0e+2;
t64 = t15.*t20.*t39.*x4.*1.376468877254218e+1;
t65 = t32.*t39;
t66 = t5.*t16.*t21.*t39.*th11.*1.894666570049486e+4;
t28 = t18.*t26.*x2;
t29 = t22.*t26;
t38 = t3.*t19.*t26.*th2.*1.0e+2;
t47 = t44+1.0;
t48 = t26+t39;
t56 = -t54;
t57 = -t55;
t67 = -t66;
t40 = -t38;
t50 = 1.0./t47;
t52 = 1.0./t48;
t58 = exp(t56);
t59 = exp(t57);
t68 = t29+t65;
t69 = t64+t67;
t51 = t50.^2;
t53 = t52.^2;
t60 = t58+1.0;
t61 = t28+t40;
t62 = 1.0./t60;
t63 = t62.^2;
mt1 = [x1.*(t52.*t61.*th3+t28.*t53.*t68.*th3.*1.0e+2+t18.*t44.*t50.*x2.*3.294e-1+t18.*t41.*t44.*t50.*x2.*3.294e+1-t18.*t41.*t45.*t51.*x2.*3.294e+1-t8.*t18.*t58.*t62.*th10.*x2.*x3.*1.720720720720721e-1+t4.*t9.*t18.*t46.*t58.*t62.*th10.*x2.*1.720720720720721e+1-t4.*t9.*t18.*t46.*t59.*t63.*th10.*x2.*1.720720720720721e+1),-x1.*(t52.*t61+t28.*t53.*t68.*1.0e+2-t18.*t44.*t50.*x2-t18.*t41.*t44.*t50.*x2.*1.0e+2+t18.*t41.*t45.*t51.*x2.*1.0e+2)];
mt2 = [x1.*(t18.*t44.*t50.*x2.*9.674e-1+t18.*t41.*t44.*t50.*x2.*9.674e+1-t18.*t41.*t45.*t51.*x2.*9.674e+1+t8.*t18.*t58.*t62.*x2.*x3.*1.720720720720721e-1-t4.*t9.*t18.*t46.*t58.*t62.*x2.*1.720720720720721e+1+t4.*t9.*t18.*t46.*t59.*t63.*x2.*1.720720720720721e+1)];
mt3 = [-x4.*(t52.*t61.*3.438e-1+t28.*t53.*t68.*3.438e+1+(t18.*t44.*t50.*x2)./1.0e+4+(t18.*t41.*t44.*t50.*x2)./1.0e+2-(t18.*t41.*t45.*t51.*x2)./1.0e+2-t8.*t18.*t58.*t62.*x2.*x3.*3.438e-1+t4.*t9.*t18.*t46.*t58.*t62.*x2.*3.438e+1-t4.*t9.*t18.*t46.*t59.*t63.*x2.*3.438e+1),0.0,x1.*(t52.*t68+t6.*t7.*t12.*3.025e-4),t6.*t7.*t12.*(x2-4.5e+2).*3.025e-4,t6.*t7.*t12.*x3.*3.025e-4,t6.*t7.*t12.*x4.*3.025e-4,t6.*t12.*(-3.025e-4)];
mt4 = [x1.*(t4.*t10.*t49.*t58.*t62.*th10.*-1.0e+2+t4.*t10.*t49.*t59.*t63.*th10.*1.0e+2+t9.*t46.*t58.*t62.*th10.*x3),0.0,-x1.*(t4.*t10.*t49.*t58.*t62.*-1.0e+2+t4.*t10.*t49.*t59.*t63.*1.0e+2+t9.*t46.*t58.*t62.*x3),-x4.*(t4.*t10.*t49.*t58.*t62.*(-9.99e+2./5.0)+t4.*t10.*t49.*t59.*t63.*(9.99e+2./5.0)+t9.*t46.*t58.*t62.*x3.*(9.99e+2./5.0e+2)),0.0,-t8.*t46.*t58.*t62.*x1.*x3,0.0,0.0,0.0,0.0];
mt5 = [x1.*(t52.*t69.*th3-t15.*t20.*t44.*t50.*x4.*4.534088481675393-t15.*t20.*t41.*t44.*t50.*x4.*4.534088481675393e+2+t15.*t20.*t41.*t45.*t51.*x4.*4.534088481675393e+2+t15.*t20.*t39.*t53.*t68.*th3.*x4.*1.376468877254218e+3+t8.*t15.*t20.*t58.*t62.*th10.*x3.*x4.*2.368518518518519-t4.*t9.*t15.*t20.*t46.*t58.*t62.*th10.*x4.*2.368518518518519e+2+t4.*t9.*t15.*t20.*t46.*t59.*t63.*th10.*x4.*2.368518518518519e+2)];
mt6 = [-x1.*(t52.*t69+t15.*t20.*t44.*t50.*x4.*1.376468877254218e+1+t15.*t20.*t41.*t44.*t50.*x4.*1.376468877254218e+3-t15.*t20.*t41.*t45.*t51.*x4.*1.376468877254218e+3+t15.*t20.*t39.*t53.*t68.*x4.*1.376468877254218e+3)];
mt7 = [-x1.*(t15.*t20.*t44.*t50.*x4.*1.33159599185573e+1+t15.*t20.*t41.*t44.*t50.*x4.*1.33159599185573e+3-t15.*t20.*t41.*t45.*t51.*x4.*1.33159599185573e+3+t8.*t15.*t20.*t58.*t62.*x3.*x4.*2.368518518518519-t4.*t9.*t15.*t20.*t46.*t58.*t62.*x4.*2.368518518518519e+2+t4.*t9.*t15.*t20.*t46.*t59.*t63.*x4.*2.368518518518519e+2)];
mt8 = [-x4.*(t52.*t69.*3.438e-1-t15.*t20.*t44.*t50.*x4.*1.376468877254218e-3-t15.*t20.*t41.*t44.*t50.*x4.*1.376468877254218e-1+t15.*t20.*t41.*t45.*t51.*x4.*1.376468877254218e-1+t15.*t20.*t39.*t53.*t68.*x4.*4.7323e+2+t8.*t15.*t20.*t58.*t62.*x3.*x4.*4.7323-t4.*t9.*t15.*t20.*t46.*t58.*t62.*x4.*4.7323e+2+t4.*t9.*t15.*t20.*t46.*t59.*t63.*x4.*4.7323e+2),0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8],5,5);
