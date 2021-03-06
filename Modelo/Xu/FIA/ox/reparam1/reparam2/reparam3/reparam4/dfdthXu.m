function dfdthSym = dfdthXu(t,in2,in3)
%DFDTHXU
%    DFDTHSYM = DFDTHXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 19:51:56

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th7 = in3(4,:);
th9 = in3(5,:);
th10 = in3(6,:);
th11 = in3(7,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = x2.^2;
t3 = x3.^2;
t4 = th9-1.0;
t5 = 1.0./th7;
t7 = 1.0./th9.^2;
t8 = 1.0./th11;
t10 = 1.0./x5;
t15 = t.*(1.1e+1./1.0e+2);
t16 = th9.*(2.3e+1./4.0e+2);
t18 = x4+1.0e-4;
t21 = th3.*9.002366863905325e+1;
t22 = x2+1.034e-1;
t26 = th3.*9.002366863905325e-1;
t27 = x3.*1.344483583855441e-1;
t28 = x3+5.019e-1;
t6 = t5.^2;
t9 = t8.^2;
t11 = t5.*x3;
t12 = 1.0./t4;
t17 = exp(t15);
t23 = 1.0./t18;
t25 = -t21;
t29 = t27+1.0;
t30 = 1.0./t22;
t32 = 1.0./t28;
t13 = t12.^2;
t14 = t11+1.0;
t24 = t12.*(2.3e+1./4.0);
t31 = t30.^2;
t33 = t32.^2;
t34 = 1.0./t29;
t36 = t32.*th2.*x3.*1.0e+2;
t19 = 1.0./t14;
t35 = t34.^2;
t37 = -t36;
t41 = t30.*t34.*x2;
t20 = t19.^2;
t38 = exp(t37);
t39 = t8.*t19.*t23.*th3.*x4;
t42 = t41.*th1;
t43 = t19.*t23.*x4.*9.002366863905325e+1;
t44 = t19.*t23.*x4.*9.002366863905325e-1;
t45 = -t41;
t46 = t19.*t21.*t23.*x4;
t48 = t41.*1.0e+2;
t49 = t19.*t23.*t26.*x4;
t53 = t19.*t23.*th3.*x4.*(-9.002366863905325e-1);
t40 = -t39;
t47 = -t42;
t50 = t42.*1.0e+2;
t51 = t42.*2.0e+2;
t52 = -t48;
t56 = t32.*t38.*x3;
t59 = t43-9.002366863905325e+1;
t60 = t3.*t33.*t38.*th2.*1.0e+2;
t62 = t44-9.002366863905325e-1;
t63 = t25+t46;
t65 = t26+t53;
t54 = -t50;
t57 = t56.*th2;
t58 = t16+t40;
t61 = -t60;
t64 = exp(t63);
t55 = exp(t54);
t66 = t12.*t58.*1.0e+2;
t67 = t13.*t58.*1.0e+2;
t76 = t38+t64;
t81 = t56+t61;
t83 = t6.*t20.*t23.*t26.*t64.*x3.*x4;
t89 = t62.*t64;
t90 = t64.*t65;
t68 = -t66;
t69 = -t67;
t74 = t41.*t55;
t75 = t42.*t55;
t77 = 1.0./t76;
t79 = t2.*t31.*t35.*t55.*th1.*1.0e+2;
t97 = t6.*t20.*t21.*t23.*t90.*x3.*x4;
t99 = t6.*t20.*t23.*t90.*th3.*x3.*x4.*(-9.002366863905325e+1);
t101 = t57+t90;
t102 = t59.*t90;
t70 = exp(t68);
t71 = t24+t69;
t78 = t77.^2;
t80 = -t79;
t103 = -t102;
t109 = t83+t99;
t72 = t12.*t70.*(2.3e+1./4.0e+2);
t82 = t8.*t12.*t19.*t23.*t70.*x4;
t84 = t9.*t12.*t19.*t23.*t70.*th3.*x4;
t85 = t55+t70;
t86 = t6.*t8.*t12.*t20.*t23.*t70.*th3.*x3.*x4;
t91 = t12.*t58.*t70;
t92 = t13.*t58.*t70;
t93 = t8.*t19.*t23.*t67.*t70.*x4;
t95 = t9.*t19.*t23.*t67.*t70.*th3.*x4;
t98 = t6.*t8.*t20.*t23.*t67.*t70.*th3.*x3.*x4;
t104 = t74+t80;
t111 = t89+t103;
t73 = -t72;
t87 = 1.0./t85;
t94 = t8.*t19.*t23.*t92.*x4.*-1.0e+2;
t96 = t9.*t19.*t23.*t92.*th3.*x4.*-1.0e+2;
t100 = t6.*t8.*t20.*t23.*t92.*th3.*x3.*x4.*-1.0e+2;
t105 = t71.*t91;
t106 = t75+t91;
t88 = t87.^2;
t107 = t82+t94;
t108 = t84+t96;
t110 = t86+t100;
t112 = t87.*t104;
t114 = t87.*t106;
t133 = t73+t92+t105;
t113 = t112.*1.0e+2;
t115 = t114.*1.0e+2;
t116 = t114.*2.0e+2;
t119 = t47+t114;
t127 = t87.*t107;
t128 = t87.*t108;
t131 = t87.*t110;
t134 = t48.*t55.*t88.*t106;
t135 = t74.*t88.*t106.*1.0e+4;
t136 = t82.*t88.*t106.*1.0e+2;
t137 = t84.*t88.*t106.*1.0e+2;
t138 = t82.*t88.*t106.*1.0e+4;
t139 = t84.*t88.*t106.*1.0e+4;
t140 = t86.*t88.*t106.*1.0e+2;
t141 = t86.*t88.*t106.*1.0e+4;
t142 = t70.*t71.*t88.*t106;
t146 = t87.*t133;
t117 = -t115;
t118 = -t116;
t129 = t127.*1.0e+2;
t130 = t128.*1.0e+2;
t132 = t131.*1.0e+2;
t143 = -t142;
t144 = t142.*1.0e+2;
t147 = t146.*1.0e+2;
t148 = t45+t112+t134;
t149 = t52+t113+t135;
t150 = t127+t136;
t151 = t128+t137;
t154 = t131+t140;
t120 = t50+t117;
t121 = t51+t118;
t145 = -t144;
t152 = t129+t138;
t153 = t130+t139;
t155 = t132+t141;
t156 = t143+t146;
t122 = exp(t120);
t123 = exp(t121);
t157 = t145+t147;
t124 = t122+1.0;
t125 = 1.0./t124;
t126 = t125.^2;
mt1 = [x1.*(th9.*(t112+t134)-t122.*t125.*t148.*th10-t122.*t125.*t149.*th10.*(t42-t114)+t123.*t126.*t149.*th10.*(t42-t114)),t45.*x1,-x1.*(t122.*t125.*t148.*3.652e-1-t122.*t125.*t148.*th10.*3.652e-1+t122.*t125.*t149.*(t42-t114).*3.652e-1-t123.*t126.*t149.*(t42-t114).*3.652e-1-t122.*t125.*t149.*th10.*(t42-t114).*3.652e-1+t123.*t126.*t149.*th10.*(t42-t114).*3.652e-1),0.0,0.0,x1.*(t77.*t81.*(3.1e+1./2.0e+2)+t56.*t78.*t101.*(3.1e+1./2.0)),0.0,-x1.*(t77.*t81+t56.*t78.*t101.*1.0e+2)];
mt2 = [-x1.*(t77.*t81.*6.427915e-1+t56.*t78.*t101.*6.427915e+1),0.0,-x1.*(t77.*t111.*(3.1e+1./2.0e+2)+t150.*th9+t59.*t64.*t78.*t101.*(3.1e+1./2.0e+2)-t122.*t125.*t150.*th10-t122.*t125.*t152.*th10.*(t42-t114)+t123.*t126.*t152.*th10.*(t42-t114)),0.0];
mt3 = [x1.*(t77.*t111+t122.*t125.*t150.*3.652e-1+t59.*t64.*t78.*t101-t122.*t125.*t150.*th10.*3.652e-1+t122.*t125.*t152.*(t42-t114).*3.652e-1-t123.*t126.*t152.*(t42-t114).*3.652e-1-t122.*t125.*t152.*th10.*(t42-t114).*3.652e-1+t123.*t126.*t152.*th10.*(t42-t114).*3.652e-1),x1.*(t77.*t111.*6.427915e-1-t19.*t23.*x4+t59.*t64.*t78.*t101.*6.427915e-1),0.0];
mt4 = [-x1.*(t77.*t109.*(3.1e+1./2.0e+2)+t154.*th9-t122.*t125.*t154.*th10-t122.*t125.*t155.*th10.*(t42-t114)+t123.*t126.*t155.*th10.*(t42-t114)+t6.*t20.*t23.*t64.*t78.*t101.*th3.*x3.*x4.*1.395366863905325e+1),0.0,x1.*(t77.*t109+t122.*t125.*t154.*3.652e-1-t122.*t125.*t154.*th10.*3.652e-1+t122.*t125.*t155.*(t42-t114).*3.652e-1-t123.*t126.*t155.*(t42-t114).*3.652e-1-t122.*t125.*t155.*th10.*(t42-t114).*3.652e-1+t123.*t126.*t155.*th10.*(t42-t114).*3.652e-1+t6.*t20.*t21.*t23.*t64.*t78.*t101.*x3.*x4)];
mt5 = [x1.*(t77.*t109.*6.427915e-1-t6.*t20.*t23.*th3.*x3.*x4+t6.*t20.*t23.*t64.*t78.*t101.*th3.*x3.*x4.*5.7866449e+1),0.0,x1.*(t114+th9.*(t142-t146)+t7.*t10.*t17.*3.025e-4-t122.*t125.*th10.*(t142-t146)-t122.*t125.*th10.*(t42-t114).*(t144-t147)+t123.*t126.*th10.*(t42-t114).*(t144-t147)-2.3e+1./4.0e+2),t7.*t10.*t17.*(x2-4.5e+2).*3.025e-4];
mt6 = [-x1.*(t122.*t125.*(t142-t146).*3.652e-1-t122.*t125.*th10.*(t142-t146).*3.652e-1+t122.*t125.*(t42-t114).*(t144-t147).*3.652e-1-t123.*t126.*(t42-t114).*(t144-t147).*3.652e-1-t122.*t125.*th10.*(t42-t114).*(t144-t147).*3.652e-1+t123.*t126.*th10.*(t42-t114).*(t144-t147).*3.652e-1)+t7.*t10.*t17.*x3.*3.025e-4,t7.*t10.*t17.*x4.*6.05e-4,t7.*t17.*(-3.025e-4),t122.*t125.*x1.*(t42-t114),0.0];
mt7 = [t122.*t125.*x1.*(t42-t114).*(-3.652e-1),0.0,0.0,x1.*(t151.*th9-t122.*t125.*t151.*th10-t122.*t125.*t153.*th10.*(t42-t114)+t123.*t126.*t153.*th10.*(t42-t114)),0.0,-x1.*(t122.*t125.*t151.*3.652e-1-t122.*t125.*t151.*th10.*3.652e-1+t122.*t125.*t153.*(t42-t114).*3.652e-1-t123.*t126.*t153.*(t42-t114).*3.652e-1-t122.*t125.*t153.*th10.*(t42-t114).*3.652e-1+t123.*t126.*t153.*th10.*(t42-t114).*3.652e-1),0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7],5,7);
