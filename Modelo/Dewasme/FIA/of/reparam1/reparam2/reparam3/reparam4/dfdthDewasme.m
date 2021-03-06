function dfdthSym = dfdthDewasme(t,in2,in3)
%DFDTHDEWASME
%    DFDTHSYM = DFDTHDEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 19:00:29

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th6 = in3(5,:);
th8 = in3(6,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th1+x2;
t3 = th5+x3;
t4 = th2.^2;
t5 = x2.^2;
t6 = x3.^2;
t7 = 1.0./th3.^2;
t8 = 1.0./x5;
t14 = x4+1.0e-4;
t16 = t.*(1.97e+2./1.0e+3);
t17 = x3+5.1e-3;
t9 = 1.0./t2;
t12 = 1.0./t3;
t18 = exp(t16);
t19 = 1.0./t14;
t22 = 1.0./t17;
t10 = t9.^2;
t11 = t9.^3;
t13 = t12.^2;
t15 = t9.*th2.*x2;
t23 = t22.^2;
t33 = t12.*t19.*x4.*1.017161140197789;
t37 = t12.*t19.*x4.*1.75025025025025e-1;
t38 = t12.*t19.*th5.*x4.*(-1.017161140197789);
t42 = t12.*t19.*x4.*8.137289121582315e+1;
t43 = t12.*t19.*th5.*x4.*(-1.75025025025025e-1);
t47 = t12.*t19.*th5.*x4.*1.627457824316463e+2;
t48 = t12.*t19.*th5.*x4.*(-8.137289121582315e+1);
t20 = t15.*8.0e+1;
t21 = t15.*1.6e+2;
t26 = t15.*1.720720720720721e-1;
t35 = t33.*th5;
t36 = t13.*t19.*th5.*x4.*1.017161140197789;
t40 = t37.*th5;
t41 = t13.*t19.*th5.*x4.*1.75025025025025e-1;
t45 = t42.*th5;
t46 = t13.*t19.*th5.*x4.*8.137289121582315e+1;
t50 = -t47;
t51 = exp(t48);
t52 = t15+t38;
t24 = -t20;
t39 = -t36;
t44 = -t41;
t49 = -t46;
t53 = t20+t48;
t54 = t21+t50;
t59 = t26+t43;
t68 = t33.*t51;
t69 = t12.*t19.*t51.*x4.*(-1.017161140197789);
t70 = t35.*t51;
t71 = t36.*t51;
t25 = exp(t24);
t55 = exp(t53);
t56 = exp(t54);
t66 = t33+t39;
t67 = t37+t44;
t72 = t42+t49;
t73 = t22.*t59.*x3.*8.0e+1;
t74 = t22.*t59.*x3.*1.6e+2;
t27 = t9.*t25.*x2;
t28 = t15.*t25;
t29 = t10.*t25.*th2.*x2;
t31 = t5.*t10.*t25.*th2.*8.0e+1;
t34 = t4.*t5.*t11.*t25.*8.0e+1;
t57 = t55+1.0;
t58 = t25+t51;
t75 = -t73;
t76 = -t74;
t83 = t70.*t72;
t30 = -t29;
t32 = -t31;
t60 = 1.0./t57;
t62 = 1.0./t58;
t77 = exp(t75);
t78 = exp(t76);
t82 = t28+t70;
t85 = t69+t71+t83;
t61 = t60.^2;
t63 = t62.^2;
t64 = t27+t32;
t65 = t30+t34;
t79 = t77+1.0;
t84 = t52.*t55.*t60.*x1;
t80 = 1.0./t79;
t81 = t80.^2;
mt1 = [-x1.*(t62.*th3.*(t29-t34)+t29.*t63.*t82.*th3.*8.0e+1+t10.*t55.*t60.*th2.*th8.*x2+t10.*t52.*t55.*t60.*th2.*th8.*x2.*8.0e+1-t10.*t52.*t56.*t61.*th2.*th8.*x2.*8.0e+1-t10.*t22.*t77.*t80.*th2.*x2.*x3.*3.016595495495495e-1+t6.*t10.*t23.*t59.*t77.*t80.*th2.*x2.*2.413276396396396e+1-t6.*t10.*t23.*t59.*t78.*t81.*th2.*x2.*2.413276396396396e+1),x1.*(t62.*(t29-t34)+t29.*t63.*t82.*8.0e+1-t10.*t55.*t60.*th2.*x2-t10.*t52.*t55.*t60.*th2.*x2.*8.0e+1+t10.*t52.*t56.*t61.*th2.*x2.*8.0e+1)];
mt2 = [-x1.*(t10.*t55.*t60.*th2.*th6.*x2+t10.*t52.*t55.*t60.*th2.*th6.*x2.*8.0e+1-t10.*t52.*t56.*t61.*th2.*th6.*x2.*8.0e+1+t10.*t22.*t77.*t80.*th2.*x2.*x3.*1.720720720720721e-1-t6.*t10.*t23.*t59.*t77.*t80.*th2.*x2.*1.376576576576577e+1+t6.*t10.*t23.*t59.*t78.*t81.*th2.*x2.*1.376576576576577e+1)];
mt3 = [x4.*(t62.*(t29-t34).*3.438e-1+t29.*t63.*t82.*2.7504e+1+(t10.*t55.*t60.*th2.*x2)./1.0e+4+(t10.*t52.*t55.*t60.*th2.*x2)./1.25e+2-(t10.*t52.*t56.*t61.*th2.*x2)./1.25e+2-t10.*t22.*t77.*t80.*th2.*x2.*x3.*3.438e-1+t6.*t10.*t23.*t59.*t77.*t80.*th2.*x2.*2.7504e+1-t6.*t10.*t23.*t59.*t78.*t81.*th2.*x2.*2.7504e+1),0.0];
mt4 = [x1.*(t62.*t64.*th3+t27.*t63.*t82.*th3.*8.0e+1+t9.*t55.*t60.*th8.*x2+t9.*t52.*t55.*t60.*th8.*x2.*8.0e+1-t9.*t52.*t56.*t61.*th8.*x2.*8.0e+1-t9.*t22.*t77.*t80.*x2.*x3.*3.016595495495495e-1+t6.*t9.*t23.*t59.*t77.*t80.*x2.*2.413276396396396e+1-t6.*t9.*t23.*t59.*t78.*t81.*x2.*2.413276396396396e+1),-x1.*(t62.*t64+t27.*t63.*t82.*8.0e+1-t9.*t55.*t60.*x2-t9.*t52.*t55.*t60.*x2.*8.0e+1+t9.*t52.*t56.*t61.*x2.*8.0e+1)];
mt5 = [x1.*(t9.*t55.*t60.*th6.*x2+t9.*t52.*t55.*t60.*th6.*x2.*8.0e+1-t9.*t52.*t56.*t61.*th6.*x2.*8.0e+1+t9.*t22.*t77.*t80.*x2.*x3.*1.720720720720721e-1-t6.*t9.*t23.*t59.*t77.*t80.*x2.*1.376576576576577e+1+t6.*t9.*t23.*t59.*t78.*t81.*x2.*1.376576576576577e+1),-x4.*(t62.*t64.*3.438e-1+t27.*t63.*t82.*2.7504e+1+(t9.*t55.*t60.*x2)./1.0e+4+(t9.*t52.*t55.*t60.*x2)./1.25e+2-(t9.*t52.*t56.*t61.*x2)./1.25e+2-t9.*t22.*t77.*t80.*x2.*x3.*3.438e-1+t6.*t9.*t23.*t59.*t77.*t80.*x2.*2.7504e+1-t6.*t9.*t23.*t59.*t78.*t81.*x2.*2.7504e+1),0.0];
mt6 = [x1.*(t62.*t82+t7.*t8.*t18.*5.4175e-4),t7.*t8.*t18.*(x2-4.5e+2).*5.4175e-4,t7.*t8.*t18.*x3.*5.4175e-4,t7.*t8.*t18.*x4.*5.4175e-4,t7.*t18.*(-5.4175e-4),-x1.*(t62.*t85.*th3+t55.*t60.*t66.*th8+t52.*t55.*t60.*t72.*th8-t52.*t56.*t61.*t72.*th8-t51.*t63.*t72.*t82.*th3-t22.*t67.*t77.*t80.*x3.*1.7531+t6.*t23.*t59.*t67.*t77.*t80.*1.40248e+2-t6.*t23.*t59.*t67.*t78.*t81.*1.40248e+2)];
mt7 = [-x1.*(-t62.*t85+t55.*t60.*t66+t52.*t55.*t60.*t72-t52.*t56.*t61.*t72+t51.*t63.*t72.*t82),-x1.*(t55.*t60.*t66.*th6+t52.*t55.*t60.*t72.*th6-t52.*t56.*t61.*t72.*th6+t22.*t67.*t77.*t80.*x3-t6.*t23.*t59.*t67.*t77.*t80.*8.0e+1+t6.*t23.*t59.*t67.*t78.*t81.*8.0e+1),x4.*(t62.*t85.*3.438e-1+(t55.*t60.*t66)./1.0e+4+(t52.*t55.*t60.*t72)./1.0e+4-(t52.*t56.*t61.*t72)./1.0e+4-t51.*t63.*t72.*t82.*3.438e-1-t22.*t67.*t77.*t80.*x3.*(9.99e+2./5.0e+2)+t6.*t23.*t59.*t67.*t77.*t80.*1.5984e+2-t6.*t23.*t59.*t67.*t78.*t81.*1.5984e+2),0.0,0.0,0.0,t84,0.0,0.0,t84,0.0];
mt8 = [0.0,0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8],5,6);
