function dfdxSym = dfdxAnane(t,in2,in3)
%DFDXANANE
%    DFDXSYM = DFDXANANE(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 23:35:45

th1 = in3(1,:);
th2 = in3(2,:);
th10 = in3(3,:);
th11 = in3(4,:);
th14 = in3(5,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th1.^2;
t3 = x2.^2;
t4 = 1.0./th10;
t5 = 1.0./x5;
t7 = x2+9.0./2.0e+2;
t8 = x4+1.0e-4;
t9 = t.*(1.97e+2./1.0e+3);
t15 = x3+4.267e-1;
t16 = x3.*6.2000124000248e-1;
t6 = t5.^2;
t10 = 1.0./t7;
t12 = exp(t9);
t13 = 1.0./t8;
t17 = t16+1.0;
t18 = 1.0./t15;
t11 = t10.^2;
t14 = t13.^2;
t19 = t18.^2;
t20 = 1.0./t17;
t26 = t4.*t5.*t12.*5.4175e-4;
t21 = t20.^2;
t22 = t20.^3;
t23 = t10.*t20.*th1;
t25 = t11.*t20.*th1.*x2;
t27 = -t26;
t24 = t23.*x2;
t28 = -t25;
t30 = t23.*4.399665625412469e-1;
t31 = t10.*t21.*th1.*x2.*6.2000124000248e-1;
t33 = t25.*4.399665625412469e-1;
t29 = t24+5.034e-1;
t32 = t24.*4.399665625412469e-1;
t34 = -t33;
t40 = t23+t28;
t35 = t32+1.0;
t36 = 1.0./t29;
t48 = t30+t34;
t37 = t36.^2;
t38 = 1.0./t35;
t41 = t23.*t36.*(1.27e+2./2.0e+2);
t42 = t10.*t21.*t36.*th1.*x2.*(5.0e+1./1.27e+2);
t45 = t24.*t36.*(1.27e+2./2.0e+2);
t46 = t25.*t36.*(1.27e+2./2.0e+2);
t39 = t38.^2;
t43 = -t41;
t44 = -t42;
t47 = -t45;
t49 = t2.*t3.*t11.*t22.*t37.*(5.0e+1./1.27e+2);
t51 = t24.*t37.*t40.*(1.27e+2./2.0e+2);
t50 = t24+t47;
t52 = t31+t44+t49;
t53 = t40+t43+t46+t51;
mt1 = [t27+t45.*th11+th10.*(t13.*t50.*x4-6.7e-3)+t18.*t38.*th2.*th14.*x3,-t24,t24.*t36.*1.241425e-1-t18.*t38.*th2.*x3,t13.*t50.*x4.*(-4.684e-1)-t18.*t38.*th2.*x3.*8.163e-1+3.13828e-3,0.0,-x1.*(t46.*th11+t51.*th11-t23.*t36.*th11.*(1.27e+2./2.0e+2)-t13.*t53.*th10.*x4+t18.*t39.*t48.*th2.*th14.*x3),t27-t23.*x1+t25.*x1];
mt2 = [x1.*(t23.*t36.*1.241425e-1-t25.*t36.*1.241425e-1-t24.*t37.*t40.*1.241425e-1+t18.*t39.*t48.*th2.*x3),-x1.*(t13.*t53.*x4.*4.684e-1-t18.*t39.*t48.*th2.*x3.*8.163e-1),0.0,x1.*(t49.*th11+t18.*t38.*th2.*th14-t13.*t52.*th10.*x4-t19.*t38.*th2.*th14.*x3-t10.*t21.*t36.*th1.*th11.*x2.*(5.0e+1./1.27e+2)+t10.*t18.*t21.*t39.*th1.*th2.*th14.*x2.*x3.*2.727798143352017e-1),t31.*x1];
mt3 = [t27-x1.*(t18.*t38.*th2-t19.*t38.*th2.*x3-t2.*t3.*t11.*t22.*t37.*7.696850393700787e-2+t10.*t21.*t36.*th1.*x2.*7.696850393700787e-2+t10.*t18.*t21.*t39.*th1.*th2.*x2.*x3.*2.727798143352017e-1),-x1.*(t18.*t38.*th2.*8.163e-1-t13.*t52.*x4.*4.684e-1-t19.*t38.*th2.*x3.*8.163e-1+t10.*t18.*t21.*t39.*th1.*th2.*x2.*x3.*2.226701624418252e-1),0.0,th10.*x1.*(t13.*t50-t14.*t50.*x4),0.0,0.0];
mt4 = [-x1.*(t13.*t50.*4.684e-1-t14.*t50.*x4.*4.684e-1)-t4.*t5.*t12.*1.0835e-3-1.8e+4,0.0,t4.*t6.*t12.*x1.*5.4175e-4,t4.*t6.*t12.*(x2-4.5e+2).*5.4175e-4,t4.*t6.*t12.*x3.*5.4175e-4,t4.*t6.*t12.*x4.*1.0835e-3,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4],5,5);
