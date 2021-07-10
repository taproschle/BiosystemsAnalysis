function dfdthSym = dfdthXu(t,in2,in3)
%DFDTHXU
%    DFDTHSYM = DFDTHXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 21:53:25

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
t12 = t.*(1.1e+1./1.0e+2);
t13 = th9.*(2.3e+1./4.0e+2);
t15 = x4+1.0e-4;
t17 = th3.*9.002366863905325e+1;
t18 = x2+1.034e-1;
t22 = th3.*9.002366863905325e-1;
t23 = x3.*1.344483583855441e-1;
t24 = x3+5.019e-1;
t6 = t5.^2;
t9 = 1.0./t3;
t11 = t8+1.0;
t14 = exp(t12);
t19 = 1.0./t15;
t21 = -t17;
t25 = t23+1.0;
t26 = 1.0./t18;
t28 = 1.0./t24;
t10 = t9.^2;
t16 = 1.0./t11;
t20 = t9.*(2.3e+1./4.0);
t27 = t26.^2;
t29 = 1.0./t25;
t31 = t28.*x3.*(3.06e+2./2.5e+1);
t30 = t29.^2;
t32 = -t31;
t34 = t5.*t16.*t19.*th3.*x4;
t36 = t16.*t19.*x4.*9.002366863905325e+1;
t37 = t26.*t29.*x2;
t38 = t16.*t17.*t19.*x4;
t40 = t16.*t19.*x4.*9.002366863905325e-1;
t41 = t16.*t19.*t22.*x4;
t43 = t16.*t19.*th3.*x4.*(-9.002366863905325e-1);
t33 = exp(t32);
t35 = -t34;
t39 = t37.*th1;
t42 = -t37;
t45 = t37.*1.0e+2;
t52 = t36-9.002366863905325e+1;
t53 = t40-9.002366863905325e-1;
t54 = t21+t38;
t56 = t22+t43;
t44 = -t39;
t46 = t39.*1.0e+2;
t47 = t39.*2.0e+2;
t48 = -t45;
t51 = t13+t35;
t55 = exp(t54);
t57 = t28.*t33.*x3.*1.224e-1;
t49 = -t46;
t58 = t9.*t51.*1.0e+2;
t59 = t10.*t51.*1.0e+2;
t68 = t33+t55;
t76 = t53.*t55;
t79 = t55.*t56;
t50 = exp(t49);
t60 = -t58;
t61 = -t59;
t69 = 1.0./t68;
t86 = t52.*t79;
t88 = t57+t79;
t62 = exp(t60);
t63 = t20+t61;
t66 = t37.*t50;
t67 = t39.*t50;
t70 = t69.^2;
t71 = t2.*t27.*t30.*t50.*th1.*1.0e+2;
t87 = -t86;
t64 = t9.*t62.*(2.3e+1./4.0e+2);
t72 = -t71;
t73 = t5.*t9.*t16.*t19.*t62.*x4;
t74 = t6.*t9.*t16.*t19.*t62.*th3.*x4;
t75 = t50+t62;
t80 = t9.*t51.*t62;
t81 = t10.*t51.*t62;
t82 = t5.*t16.*t19.*t59.*t62.*x4;
t83 = t6.*t16.*t19.*t59.*t62.*th3.*x4;
t94 = t76+t87;
t65 = -t64;
t77 = 1.0./t75;
t84 = t5.*t16.*t19.*t81.*x4.*-1.0e+2;
t85 = t6.*t16.*t19.*t81.*th3.*x4.*-1.0e+2;
t89 = t66+t72;
t90 = t63.*t80;
t91 = t67+t80;
t78 = t77.^2;
t92 = t73+t84;
t93 = t74+t85;
t95 = t77.*t89;
t97 = t77.*t91;
t114 = t65+t81+t90;
t96 = t95.*1.0e+2;
t98 = t97.*1.0e+2;
t99 = t97.*2.0e+2;
t102 = t44+t97;
t108 = t77.*t92;
t111 = t77.*t93;
t115 = t45.*t50.*t78.*t91;
t116 = t66.*t78.*t91.*1.0e+4;
t117 = t73.*t78.*t91.*1.0e+2;
t118 = t74.*t78.*t91.*1.0e+2;
t119 = t73.*t78.*t91.*1.0e+4;
t120 = t74.*t78.*t91.*1.0e+4;
t121 = t62.*t63.*t78.*t91;
t125 = t77.*t114;
t100 = -t98;
t101 = -t99;
t112 = t108.*1.0e+2;
t113 = t111.*1.0e+2;
t122 = -t121;
t123 = t121.*1.0e+2;
t126 = t125.*1.0e+2;
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
mt1 = [x1.*(th9.*(t95+t115)-t105.*t109.*t127.*th10-t105.*t109.*t128.*th10.*(t39-t97)+t106.*t110.*t128.*th10.*(t39-t97)),t42.*x1,-x1.*(t105.*t109.*t127.*3.652e-1-t105.*t109.*t127.*th10.*3.652e-1+t105.*t109.*t128.*(t39-t97).*3.652e-1-t106.*t110.*t128.*(t39-t97).*3.652e-1-t105.*t109.*t128.*th10.*(t39-t97).*3.652e-1+t106.*t110.*t128.*th10.*(t39-t97).*3.652e-1),0.0,0.0];
mt2 = [-x1.*(t69.*t94.*(3.1e+1./2.0e+2)+t129.*th9+t52.*t55.*t70.*t88.*(3.1e+1./2.0e+2)-t105.*t109.*t129.*th10-t105.*t109.*t131.*th10.*(t39-t97)+t106.*t110.*t131.*th10.*(t39-t97)),0.0,x1.*(t69.*t94+t105.*t109.*t129.*3.652e-1+t52.*t55.*t70.*t88-t105.*t109.*t129.*th10.*3.652e-1+t105.*t109.*t131.*(t39-t97).*3.652e-1-t106.*t110.*t131.*(t39-t97).*3.652e-1-t105.*t109.*t131.*th10.*(t39-t97).*3.652e-1+t106.*t110.*t131.*th10.*(t39-t97).*3.652e-1)];
mt3 = [x1.*(t69.*t94.*6.427915e-1-t16.*t19.*x4+t52.*t55.*t70.*t88.*6.427915e-1),0.0,x1.*(t97+th9.*(t121-t125)+t4.*t7.*t14.*3.025e-4-t105.*t109.*th10.*(t121-t125)-t105.*t109.*th10.*(t39-t97).*(t123-t126)+t106.*t110.*th10.*(t39-t97).*(t123-t126)-2.3e+1./4.0e+2),t4.*t7.*t14.*(x2-4.5e+2).*3.025e-4];
mt4 = [-x1.*(t105.*t109.*(t121-t125).*3.652e-1-t105.*t109.*th10.*(t121-t125).*3.652e-1+t105.*t109.*(t39-t97).*(t123-t126).*3.652e-1-t106.*t110.*(t39-t97).*(t123-t126).*3.652e-1-t105.*t109.*th10.*(t39-t97).*(t123-t126).*3.652e-1+t106.*t110.*th10.*(t39-t97).*(t123-t126).*3.652e-1)+t4.*t7.*t14.*x3.*3.025e-4,t4.*t7.*t14.*x4.*6.05e-4,t4.*t14.*(-3.025e-4),t105.*t109.*x1.*(t39-t97),0.0];
mt5 = [t105.*t109.*x1.*(t39-t97).*(-3.652e-1),0.0,0.0,x1.*(t130.*th9-t105.*t109.*t130.*th10-t105.*t109.*t132.*th10.*(t39-t97)+t106.*t110.*t132.*th10.*(t39-t97)),0.0,-x1.*(t105.*t109.*t130.*3.652e-1-t105.*t109.*t130.*th10.*3.652e-1+t105.*t109.*t132.*(t39-t97).*3.652e-1-t106.*t110.*t132.*(t39-t97).*3.652e-1-t105.*t109.*t132.*th10.*(t39-t97).*3.652e-1+t106.*t110.*t132.*th10.*(t39-t97).*3.652e-1),0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5],5,5);
