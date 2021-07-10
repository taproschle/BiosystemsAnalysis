function dfdthSym = dfdthAnane(t,in2,in3)
%DFDTHANANE
%    DFDTHSYM = DFDTHANANE(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 16:52:58

th1 = in3(1,:);
th2 = in3(2,:);
th5 = in3(3,:);
th6 = in3(4,:);
th7 = in3(5,:);
th8 = in3(6,:);
th10 = in3(7,:);
th11 = in3(8,:);
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
t18 = t.*(1.1e+1./1.0e+2);
t20 = x4+1.0e-4;
t8 = t7.^2;
t11 = t7.*x3;
t12 = 1.0./t2;
t15 = 1.0./t3;
t19 = exp(t18);
t24 = 1.0./t20;
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
t31 = t8.*t12.*t22.*th1.*x2.*x3;
t28 = t26+th7;
t32 = t26.*4.399665625412469e-1;
t29 = 1.0./t28;
t33 = t32+1.0;
t30 = t29.^2;
t34 = 1.0./t33;
t36 = t25.*t29.*(1.27e+2./2.0e+2);
t38 = t27.*t29.*(1.27e+2./2.0e+2);
t39 = t26.*t29.*x1.*(1.27e+2./2.0e+2);
t43 = t29.*t31.*(1.27e+2./2.0e+2);
t35 = t34.^2;
t37 = -t36;
t40 = -t38;
t41 = t5.*t13.*t22.*t30.*th1.*(1.27e+2./2.0e+2);
t42 = t4.*t5.*t14.*t22.*t30.*(1.27e+2./2.0e+2);
t44 = -t43;
t45 = t4.*t5.*t8.*t13.*t23.*t30.*x3.*(1.27e+2./2.0e+2);
t46 = t25+t37+t41;
t47 = t27+t40+t42;
t48 = t31+t44+t45;
mt1 = [x1.*(t36.*th11+t24.*t46.*th10.*x4-t5.*t13.*t22.*t30.*th1.*th11.*(1.27e+2./2.0e+2)-t15.*t25.*t35.*th2.*th14.*x3.*4.399665625412469e-1),-t25.*x1,x1.*(t36.*th15+t15.*t25.*t35.*th2.*x3.*4.399665625412469e-1-t5.*t13.*t22.*t30.*th1.*th15.*(1.27e+2./2.0e+2)),-x1.*(t24.*t46.*x4.*4.684e-1-t15.*t25.*t35.*th2.*x3.*3.591447050024198e-1),0.0,t15.*t34.*th14.*x1.*x3,0.0,-t15.*t34.*x1.*x3,t15.*t34.*x1.*x3.*(-8.163e-1),0.0];
mt2 = [-x1.*(t38.*th11+t24.*t47.*th10.*x4-t4.*t5.*t14.*t22.*t30.*th11.*(1.27e+2./2.0e+2)-t15.*t27.*t35.*th2.*th14.*x3.*4.399665625412469e-1),t27.*x1,-x1.*(t38.*th15+t15.*t27.*t35.*th2.*x3.*4.399665625412469e-1-t4.*t5.*t14.*t22.*t30.*th15.*(1.27e+2./2.0e+2)),x1.*(t24.*t47.*x4.*4.684e-1-t15.*t27.*t35.*th2.*x3.*3.591447050024198e-1),0.0,-t16.*t34.*th2.*th14.*x1.*x3,0.0,t16.*t34.*th2.*x1.*x3,t16.*t34.*th2.*x1.*x3.*8.163e-1,0.0,-x1.*(t26.*t30.*th11.*(1.27e+2./2.0e+2)-t24.*t26.*t30.*th10.*x4.*(1.27e+2./2.0e+2)),0.0];
mt3 = [t26.*t30.*th15.*x1.*(-1.27e+2./2.0e+2),t24.*t26.*t30.*x1.*x4.*(-2.97434e-1),0.0,x1.*(t43.*th11+t24.*t48.*th10.*x4-t4.*t5.*t8.*t13.*t23.*t30.*th11.*x3.*(1.27e+2./2.0e+2)-t6.*t8.*t12.*t15.*t22.*t35.*th1.*th2.*th14.*x2.*4.399665625412469e-1),-t31.*x1,x1.*(t43.*th15-t4.*t5.*t8.*t13.*t23.*t30.*th15.*x3.*(1.27e+2./2.0e+2)+t6.*t8.*t12.*t15.*t22.*t35.*th1.*th2.*x2.*4.399665625412469e-1),-x1.*(t24.*t48.*x4.*4.684e-1-t6.*t8.*t12.*t15.*t22.*t35.*th1.*th2.*x2.*3.591447050024198e-1),0.0];
mt4 = [x1.*(t24.*x4.*(t26-t26.*t29.*(1.27e+2./2.0e+2))+t9.*t10.*t19.*3.025e-4-6.7e-3),t9.*t10.*t19.*(x2-4.5e+2).*3.025e-4,t9.*t10.*t19.*x3.*3.025e-4,t9.*t10.*t19.*x4.*6.05e-4,t9.*t19.*(-3.025e-4),t39,0.0,0.0,0.0,0.0,t15.*t34.*th2.*x1.*x3,0.0,0.0,0.0,0.0,0.0,0.0,t39,0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4],5,10);