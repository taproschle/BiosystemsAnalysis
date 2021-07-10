function Xdot = Anane(t,in2,in3)
%ANANE
%    XDOT = ANANE(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 21:49:27

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
t2 = 1.0./th10;
t3 = 1.0./x5;
t4 = x2+9.0./2.0e+2;
t5 = x4+1.0e-4;
t6 = t.*(1.97e+2./1.0e+3);
t10 = x3+4.267e-1;
t11 = x3.*6.2000124000248e-1;
t7 = 1.0./t4;
t8 = exp(t6);
t9 = 1.0./t5;
t12 = t11+1.0;
t13 = 1.0./t10;
t14 = 1.0./t12;
t15 = t7.*t14.*th1.*x2;
t16 = t15+5.034e-1;
t17 = t15.*4.399665625412469e-1;
t18 = t17+1.0;
t19 = 1.0./t16;
t20 = 1.0./t18;
t21 = t15.*t19.*(1.27e+2./2.0e+2);
t22 = -t21;
t23 = t15+t22;
mt1 = [x1.*(t21.*th11+th10.*(t9.*t23.*x4-6.7e-3)-t2.*t3.*t8.*5.4175e-4+t13.*t20.*th2.*th14.*x3),-t15.*x1-t2.*t3.*t8.*(x2-4.5e+2).*5.4175e-4,x1.*(t15.*t19.*1.241425e-1-t13.*t20.*th2.*x3)-t2.*t3.*t8.*x3.*5.4175e-4,x4.*-1.8e+4-x1.*(t9.*t23.*x4.*4.684e-1+t13.*t20.*th2.*x3.*8.163e-1-3.13828e-3)-t2.*t3.*t8.*x4.*1.0835e-3+6.3e+2];
mt2 = [t2.*t8.*5.4175e-4];
Xdot = reshape([mt1,mt2],5,1);
