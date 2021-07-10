function dfdxSym = dfdxAnane(t,in2,in3)
%DFDXANANE
%    DFDXSYM = DFDXANANE(T,IN2,IN3)

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
t4 = th1.^2;
t5 = x2.^2;
t6 = 1.0./th8;
t7 = 1.0./th9;
t8 = 1.0./th10;
t9 = 1.0./x5;
t17 = x4+1.0e-4;
t18 = t.*(1.97e+2./1.0e+3);
t10 = t9.^2;
t11 = t6.*x3;
t12 = 1.0./t2;
t14 = 1.0./t3;
t19 = exp(t18);
t23 = 1.0./t17;
t13 = t12.^2;
t15 = t14.^2;
t16 = t11+1.0;
t24 = t23.^2;
t40 = t8.*t9.*t19.*5.4175e-4;
t20 = 1.0./t16;
t41 = -t40;
t21 = t20.^2;
t22 = t20.^3;
t25 = t12.*t20.*th1;
t27 = t13.*t20.*th1.*x2;
t26 = t25.*x2;
t29 = -t27;
t30 = t7.*t25;
t31 = t6.*t12.*t21.*th1.*x2;
t33 = t7.*t27;
t28 = t26+th7;
t32 = t7.*t26;
t36 = t7.*t29;
t42 = t25+t29;
t34 = 1.0./t28;
t37 = t32+1.0;
t49 = t30+t36;
t35 = t34.^2;
t38 = 1.0./t37;
t43 = t25.*t34.*th3;
t44 = t26.*t34.*th3;
t45 = t27.*t34.*th3;
t48 = t31.*t34.*th3;
t39 = t38.^2;
t46 = -t43;
t47 = -t44;
t50 = -t48;
t51 = t4.*t5.*t6.*t13.*t22.*t35.*th3;
t57 = t26.*t35.*t42.*th3;
t52 = t26+t47;
t59 = t31+t50+t51;
t61 = t42+t45+t46+t57;
t53 = t23.*t52;
t55 = t24.*t52.*x4;
t54 = t53.*x4;
t56 = -t55;
t58 = t54-6.7e-3;
t60 = t53+t56;
mt1 = [t41+t44.*th11+t58.*th10+t14.*t38.*th2.*th14.*x3,-t26,t44.*th15-t14.*t38.*th2.*x3,-t58.*th12-t14.*t38.*th2.*th13.*x3,0.0,-x1.*(t45.*th11+t46.*th11+t57.*th11-t23.*t61.*th10.*x4+t14.*t39.*t49.*th2.*th14.*x3),t41-t25.*x1+t27.*x1,x1.*(t43.*th15-t57.*th15+t29.*t34.*th3.*th15+t14.*t39.*t49.*th2.*x3),-x1.*(t23.*t61.*th12.*x4-t14.*t39.*t49.*th2.*th13.*x3),0.0,-x1.*(t48.*th11-t51.*th11-t14.*t38.*th2.*th14+t23.*t59.*th10.*x4+t15.*t38.*th2.*th14.*x3-t7.*t11.*t12.*t14.*t21.*t39.*th1.*th2.*th14.*x2),t31.*x1];
mt2 = [t41-x1.*(t48.*th15-t51.*th15+t14.*t38.*th2-t15.*t38.*th2.*x3+t7.*t11.*t12.*t14.*t21.*t39.*th1.*th2.*x2),-x1.*(t14.*t38.*th2.*th13-t23.*t59.*th12.*x4-t15.*t38.*th2.*th13.*x3+t7.*t11.*t12.*t14.*t21.*t39.*th1.*th2.*th13.*x2),0.0,t60.*th10.*x1,0.0,0.0,t8.*t9.*t19.*(-1.0835e-3)-t60.*th12.*x1-1.8e+4,0.0,t8.*t10.*t19.*x1.*5.4175e-4,t8.*t10.*t19.*(x2-4.5e+2).*5.4175e-4,t8.*t10.*t19.*x3.*5.4175e-4,t8.*t10.*t19.*x4.*1.0835e-3,0.0];
dfdxSym = reshape([mt1,mt2],5,5);
