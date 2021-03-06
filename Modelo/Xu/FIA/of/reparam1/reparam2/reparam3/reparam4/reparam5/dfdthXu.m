function dfdthSym = dfdthXu(t,in2,in3)
%DFDTHXU
%    DFDTHSYM = DFDTHXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    10-Jul-2021 00:14:47

th1 = in3(1,:);
th3 = in3(2,:);
th9 = in3(3,:);
th10 = in3(4,:);
th11 = in3(5,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = x2.^2;
t3 = th9-1.0;
t4 = 1.0./th9.^2;
t5 = 1.0./th11;
t7 = 1.0./x5;
t8 = x3./4.0;
t12 = th9.*(2.3e+1./4.0e+2);
t13 = x4+1.0e-4;
t15 = t.*(1.97e+2./1.0e+3);
t17 = x2+1.034e-1;
t20 = th3.*9.002366863905325e-1;
t21 = th3.*7.20189349112426e+1;
t22 = x3.*1.344483583855441e-1;
t23 = x3+5.019e-1;
t6 = t5.^2;
t9 = 1.0./t3;
t11 = t8+1.0;
t16 = exp(t15);
t18 = 1.0./t13;
t24 = -t21;
t25 = t22+1.0;
t26 = 1.0./t17;
t28 = 1.0./t23;
t10 = t9.^2;
t14 = 1.0./t11;
t19 = t9.*(2.3e+1./5.0);
t27 = t26.^2;
t29 = 1.0./t25;
t31 = t28.*x3.*9.792;
t30 = t29.^2;
t32 = -t31;
t34 = t5.*t14.*t18.*th3.*x4;
t36 = t26.*t29.*x2;
t38 = t14.*t18.*x4.*9.002366863905325e-1;
t39 = t14.*t18.*x4.*7.20189349112426e+1;
t40 = t14.*t18.*t20.*x4;
t41 = t14.*t18.*t21.*x4;
t43 = t14.*t18.*th3.*x4.*(-9.002366863905325e-1);
t33 = exp(t32);
t35 = -t34;
t37 = t36.*th1;
t42 = -t36;
t45 = t36.*8.0e+1;
t52 = t38-9.002366863905325e-1;
t53 = t39-7.20189349112426e+1;
t54 = t20+t43;
t55 = t24+t41;
t44 = -t37;
t46 = t37.*8.0e+1;
t47 = t37.*1.6e+2;
t48 = -t45;
t51 = t12+t35;
t56 = exp(t55);
t59 = t28.*t33.*x3.*1.224e-1;
t49 = -t46;
t57 = t9.*t51.*8.0e+1;
t58 = t10.*t51.*8.0e+1;
t68 = t33+t56;
t78 = t52.*t56;
t79 = t54.*t56;
t50 = exp(t49);
t60 = -t57;
t61 = -t58;
t71 = 1.0./t68;
t86 = t53.*t79;
t89 = t59+t79;
t62 = exp(t60);
t63 = t19+t61;
t66 = t36.*t50;
t67 = t37.*t50;
t69 = t2.*t27.*t30.*t50.*th1.*8.0e+1;
t72 = t71.^2;
t87 = -t86;
t64 = t9.*t62.*(2.3e+1./4.0e+2);
t70 = -t69;
t73 = t5.*t9.*t14.*t18.*t62.*x4;
t74 = t6.*t9.*t14.*t18.*t62.*th3.*x4;
t75 = t50+t62;
t80 = t9.*t51.*t62;
t81 = t10.*t51.*t62;
t82 = t5.*t14.*t18.*t58.*t62.*x4;
t83 = t6.*t14.*t18.*t58.*t62.*th3.*x4;
t94 = t78+t87;
t65 = -t64;
t76 = 1.0./t75;
t84 = t5.*t14.*t18.*t81.*x4.*-8.0e+1;
t85 = t6.*t14.*t18.*t81.*th3.*x4.*-8.0e+1;
t88 = t66+t70;
t90 = t63.*t80;
t91 = t67+t80;
t77 = t76.^2;
t92 = t73+t84;
t93 = t74+t85;
t95 = t76.*t88;
t97 = t76.*t91;
t114 = t65+t81+t90;
t96 = t95.*8.0e+1;
t98 = t97.*8.0e+1;
t99 = t97.*1.6e+2;
t102 = t44+t97;
t108 = t76.*t92;
t111 = t76.*t93;
t115 = t45.*t50.*t77.*t91;
t116 = t66.*t77.*t91.*6.4e+3;
t117 = t73.*t77.*t91.*8.0e+1;
t118 = t74.*t77.*t91.*8.0e+1;
t119 = t73.*t77.*t91.*6.4e+3;
t120 = t74.*t77.*t91.*6.4e+3;
t121 = t62.*t63.*t77.*t91;
t125 = t76.*t114;
t100 = -t98;
t101 = -t99;
t112 = t108.*8.0e+1;
t113 = t111.*8.0e+1;
t122 = -t121;
t123 = t121.*8.0e+1;
t126 = t125.*8.0e+1;
t127 = t42+t95+t115;
t128 = t48+t96+t116;
t129 = t108+t117;
t130 = t111+t118;
t103 = t46+t100;
t104 = t47+t101;
t124 = -t123;
t131 = t112+t119;
t132 = t113+t120;
t133 = t122+t125;
t105 = exp(t103);
t106 = exp(t104);
t134 = t124+t126;
t107 = t105+1.0;
t109 = 1.0./t107;
t110 = t109.^2;
mt1 = [x1.*(th9.*(t95+t115)-t105.*t109.*t127.*th10-t105.*t109.*t128.*th10.*(t37-t97)+t106.*t110.*t128.*th10.*(t37-t97)),t42.*x1,-x1.*(t105.*t109.*t127.*3.652e-1-t105.*t109.*t127.*th10.*3.652e-1+t105.*t109.*t128.*(t37-t97).*3.652e-1-t106.*t110.*t128.*(t37-t97).*3.652e-1-t105.*t109.*t128.*th10.*(t37-t97).*3.652e-1+t106.*t110.*t128.*th10.*(t37-t97).*3.652e-1),0.0,0.0];
mt2 = [-x1.*(t71.*t94.*(3.1e+1./2.0e+2)+t129.*th9+t53.*t56.*t72.*t89.*(3.1e+1./2.0e+2)-t105.*t109.*t129.*th10-t105.*t109.*t131.*th10.*(t37-t97)+t106.*t110.*t131.*th10.*(t37-t97)),0.0,x1.*(t71.*t94+t105.*t109.*t129.*3.652e-1+t53.*t56.*t72.*t89-t105.*t109.*t129.*th10.*3.652e-1+t105.*t109.*t131.*(t37-t97).*3.652e-1-t106.*t110.*t131.*(t37-t97).*3.652e-1-t105.*t109.*t131.*th10.*(t37-t97).*3.652e-1+t106.*t110.*t131.*th10.*(t37-t97).*3.652e-1)];
mt3 = [x1.*(t71.*t94.*6.427915e-1-t14.*t18.*x4+t53.*t56.*t72.*t89.*6.427915e-1),0.0,x1.*(t97+th9.*(t121-t125)+t4.*t7.*t16.*5.4175e-4-t105.*t109.*th10.*(t121-t125)-t105.*t109.*th10.*(t37-t97).*(t123-t126)+t106.*t110.*th10.*(t37-t97).*(t123-t126)-2.3e+1./4.0e+2),t4.*t7.*t16.*(x2-4.5e+2).*5.4175e-4];
mt4 = [-x1.*(t105.*t109.*(t121-t125).*3.652e-1-t105.*t109.*th10.*(t121-t125).*3.652e-1+t105.*t109.*(t37-t97).*(t123-t126).*3.652e-1-t106.*t110.*(t37-t97).*(t123-t126).*3.652e-1-t105.*t109.*th10.*(t37-t97).*(t123-t126).*3.652e-1+t106.*t110.*th10.*(t37-t97).*(t123-t126).*3.652e-1)+t4.*t7.*t16.*x3.*5.4175e-4,t4.*t7.*t16.*x4.*1.0835e-3,t4.*t16.*(-5.4175e-4)];
mt5 = [t105.*t109.*x1.*(t37-t97),0.0,t105.*t109.*x1.*(t37-t97).*(-3.652e-1),0.0,0.0,x1.*(t130.*th9-t105.*t109.*t130.*th10-t105.*t109.*t132.*th10.*(t37-t97)+t106.*t110.*t132.*th10.*(t37-t97)),0.0,-x1.*(t105.*t109.*t130.*3.652e-1-t105.*t109.*t130.*th10.*3.652e-1+t105.*t109.*t132.*(t37-t97).*3.652e-1-t106.*t110.*t132.*(t37-t97).*3.652e-1-t105.*t109.*t132.*th10.*(t37-t97).*3.652e-1+t106.*t110.*t132.*th10.*(t37-t97).*3.652e-1),0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5],5,5);
