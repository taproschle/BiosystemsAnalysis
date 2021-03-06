function dfdthSym = dfdthAnane(t,in2,in3)
%DFDTHANANE
%    DFDTHSYM = DFDTHANANE(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 21:44:07

th1 = in3(1,:);
th2 = in3(2,:);
th5 = in3(3,:);
th6 = in3(4,:);
th8 = in3(5,:);
th10 = in3(6,:);
th11 = in3(7,:);
th12 = in3(8,:);
th14 = in3(9,:);
th15 = in3(10,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th5+x2;
t3 = th6+x3;
t4 = th1.^2;
t5 = x2.^2;
t6 = x3.^2;
t7 = 1.0./th8;
t9 = 1.0./th10.^2;
t10 = 1.0./x5;
t18 = x4+1.0e-4;
t19 = t.*(1.97e+2./1.0e+3);
t8 = t7.^2;
t11 = t7.*x3;
t12 = 1.0./t2;
t15 = 1.0./t3;
t20 = exp(t19);
t24 = 1.0./t18;
t13 = t12.^2;
t14 = t12.^3;
t16 = t15.^2;
t17 = t11+1.0;
t21 = 1.0./t17;
t22 = t21.^2;
t23 = t21.^3;
t25 = t12.*t21.*x2;
t27 = t13.*t21.*th1.*x2;
t26 = t25.*th1;
t28 = t8.*t12.*t22.*th1.*x2.*x3;
t29 = t26+5.034e-1;
t30 = t26.*4.399665625412469e-1;
t31 = 1.0./t29;
t33 = t30+1.0;
t32 = t31.^2;
t34 = 1.0./t33;
t36 = t25.*t31.*(1.27e+2./2.0e+2);
t38 = t26.*t31.*(1.27e+2./2.0e+2);
t39 = t27.*t31.*(1.27e+2./2.0e+2);
t45 = t28.*t31.*(1.27e+2./2.0e+2);
t35 = t34.^2;
t37 = -t36;
t40 = -t38;
t41 = -t39;
t42 = t38.*x1;
t43 = t5.*t13.*t22.*t32.*th1.*(1.27e+2./2.0e+2);
t44 = t4.*t5.*t14.*t22.*t32.*(1.27e+2./2.0e+2);
t46 = -t45;
t47 = t4.*t5.*t8.*t13.*t23.*t32.*x3.*(1.27e+2./2.0e+2);
t48 = t26+t40;
t50 = t25+t37+t43;
t51 = t27+t41+t44;
t52 = t28+t46+t47;
t49 = t24.*t48.*x4;
mt1 = [x1.*(t36.*th11+t24.*t50.*th10.*x4-t5.*t13.*t22.*t32.*th1.*th11.*(1.27e+2./2.0e+2)-t15.*t25.*t35.*th2.*th14.*x3.*4.399665625412469e-1),-t25.*x1,x1.*(t36.*th15+t15.*t25.*t35.*th2.*x3.*4.399665625412469e-1-t5.*t13.*t22.*t32.*th1.*th15.*(1.27e+2./2.0e+2)),-x1.*(t24.*t50.*th12.*x4-t15.*t25.*t35.*th2.*x3.*3.591447050024198e-1),0.0,t15.*t34.*th14.*x1.*x3,0.0,-t15.*t34.*x1.*x3,t15.*t34.*x1.*x3.*(-8.163e-1),0.0,-x1.*(t39.*th11+t24.*t51.*th10.*x4-t4.*t5.*t14.*t22.*t32.*th11.*(1.27e+2./2.0e+2)-t15.*t27.*t35.*th2.*th14.*x3.*4.399665625412469e-1),t27.*x1];
mt2 = [-x1.*(t39.*th15+t15.*t27.*t35.*th2.*x3.*4.399665625412469e-1-t4.*t5.*t14.*t22.*t32.*th15.*(1.27e+2./2.0e+2)),x1.*(t24.*t51.*th12.*x4-t15.*t27.*t35.*th2.*x3.*3.591447050024198e-1),0.0,-t16.*t34.*th2.*th14.*x1.*x3,0.0,t16.*t34.*th2.*x1.*x3,t16.*t34.*th2.*x1.*x3.*8.163e-1,0.0,x1.*(t45.*th11+t24.*t52.*th10.*x4-t4.*t5.*t8.*t13.*t23.*t32.*th11.*x3.*(1.27e+2./2.0e+2)-t6.*t8.*t12.*t15.*t22.*t35.*th1.*th2.*th14.*x2.*4.399665625412469e-1),-t28.*x1,x1.*(t45.*th15-t4.*t5.*t8.*t13.*t23.*t32.*th15.*x3.*(1.27e+2./2.0e+2)+t6.*t8.*t12.*t15.*t22.*t35.*th1.*th2.*x2.*4.399665625412469e-1)];
mt3 = [-x1.*(t24.*t52.*th12.*x4-t6.*t8.*t12.*t15.*t22.*t35.*th1.*th2.*x2.*3.591447050024198e-1),0.0,x1.*(t49+t9.*t10.*t20.*5.4175e-4-6.7e-3),t9.*t10.*t20.*(x2-4.5e+2).*5.4175e-4,t9.*t10.*t20.*x3.*5.4175e-4,t9.*t10.*t20.*x4.*1.0835e-3,t9.*t20.*(-5.4175e-4),t42,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-x1.*(t49-6.7e-3),0.0,t15.*t34.*th2.*x1.*x3,0.0,0.0,0.0,0.0,0.0,0.0,t42,0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3],5,10);
