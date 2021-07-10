function dfdxSym = dfdxXu(t,in2,in3)
%DFDXXU
%    DFDXSYM = DFDXXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 19:51:53

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
t2 = th1.^2;
t3 = x2.^2;
t4 = th9-1.0;
t5 = 1.0./th7;
t6 = 1.0./th9;
t7 = 1.0./th11;
t8 = 1.0./x5;
t14 = t.*(1.1e+1./1.0e+2);
t15 = th9.*(2.3e+1./4.0e+2);
t17 = x4+1.0e-4;
t20 = th3.*9.002366863905325e+1;
t21 = x2+1.034e-1;
t25 = th3.*9.002366863905325e-1;
t26 = x3.*1.344483583855441e-1;
t27 = x3+5.019e-1;
t9 = t8.^2;
t10 = t5.*x3;
t11 = 1.0./t4;
t16 = exp(t14);
t22 = 1.0./t17;
t24 = -t20;
t28 = t26+1.0;
t29 = 1.0./t21;
t31 = 1.0./t27;
t12 = t11.^2;
t13 = t10+1.0;
t23 = t22.^2;
t30 = t29.^2;
t32 = t31.^2;
t33 = 1.0./t28;
t36 = t31.*th2.*1.0e+2;
t39 = t31.*th2.*x3.*-1.0e+2;
t42 = t6.*t8.*t16.*3.025e-4;
t18 = 1.0./t13;
t34 = t33.^2;
t35 = t33.^3;
t37 = t36.*x3;
t38 = t32.*th2.*x3.*1.0e+2;
t41 = exp(t39);
t43 = -t42;
t49 = t29.*t33.*th1;
t51 = t30.*t33.*th1.*x2;
t19 = t18.^2;
t40 = -t38;
t44 = t7.*t18.*t22.*th3;
t46 = t7.*t18.*t23.*th3.*x4;
t50 = t49.*x2;
t52 = t18.*t20.*t22;
t53 = t18.*t22.*t25;
t55 = t18.*t20.*t23.*x4;
t57 = -t51;
t58 = t49.*1.0e+2;
t59 = t18.*t23.*th3.*x4.*(-9.002366863905325e+1);
t61 = t18.*t23.*t25.*x4;
t63 = t51.*1.0e+2;
t65 = t18.*t22.*th3.*x4.*(-9.002366863905325e-1);
t66 = t18.*t23.*th3.*x4.*(-9.002366863905325e-1);
t70 = t31.*t41.*th2;
t72 = t32.*t41.*th2.*x3;
t76 = t29.*t34.*th1.*x2.*1.344483583855441e-1;
t77 = t29.*t34.*th1.*x2.*1.344483583855441e+1;
t45 = t44.*x4;
t48 = -t46;
t54 = t52.*x4;
t56 = -t50;
t60 = t53.*x4;
t62 = t50.*1.0e+2;
t64 = t50.*2.0e+2;
t68 = -t63;
t71 = t70.*x3;
t74 = -t70;
t75 = t36+t40;
t80 = t25+t65;
t86 = t52+t59;
t90 = t53+t66;
t47 = -t45;
t67 = -t62;
t78 = t24+t54;
t84 = t44+t48;
t85 = t41.*t75;
t93 = t58+t68;
t97 = t71.*t75;
t69 = exp(t67);
t73 = t15+t47;
t79 = exp(t78);
t81 = t11.*t73.*1.0e+2;
t87 = t49.*t69;
t88 = t50.*t69;
t89 = t51.*t69;
t91 = t41+t79;
t96 = t69.*t76;
t98 = t29.*t34.*t69.*th1.*x2.*(-1.344483583855441e-1);
t99 = t69.*t77;
t100 = t5.*t19.*t20.*t22.*t79.*x4;
t101 = t5.*t19.*t22.*t25.*t79.*x4;
t102 = t2.*t3.*t30.*t35.*t69.*1.344483583855441e+1;
t103 = t5.*t19.*t22.*t79.*th3.*x4.*(-9.002366863905325e-1);
t110 = t79.*t80;
t112 = t79.*t90;
t82 = -t81;
t92 = -t87;
t94 = 1.0./t91;
t114 = t80.*t100;
t117 = t71+t110;
t118 = t88.*t93;
t119 = t85+t100;
t122 = t86.*t110;
t83 = exp(t82);
t95 = t94.^2;
t123 = -t122;
t132 = t89+t92+t118;
t149 = t72+t74+t97+t103+t114;
t104 = t69+t83;
t105 = t5.*t7.*t11.*t19.*t22.*t83.*th3.*x4;
t111 = t11.*t73.*t83;
t113 = t11.*t83.*t84;
t115 = t5.*t7.*t12.*t19.*t22.*t73.*t83.*th3.*x4.*1.0e+2;
t120 = t12.*t73.*t83.*t84.*1.0e+2;
t133 = t112+t123;
t106 = 1.0./t104;
t108 = t105.*1.0e+2;
t116 = -t115;
t121 = -t120;
t124 = t88+t111;
t107 = t106.^2;
t109 = -t108;
t126 = t106.*t124;
t131 = t113+t121;
t147 = t106.*t132;
t153 = t98+t102+t105+t116;
t125 = t99+t109;
t127 = t126.*1.0e+2;
t128 = t126.*2.0e+2;
t134 = t56+t126;
t142 = t69.*t93.*t107.*t124;
t145 = t106.*t131;
t150 = t147.*1.0e+2;
t151 = t107.*t113.*t124.*1.0e+2;
t152 = t107.*t113.*t124.*1.0e+4;
t158 = t106.*t153;
t129 = -t127;
t130 = -t128;
t143 = -t142;
t144 = t142.*1.0e+2;
t148 = t145.*1.0e+2;
t154 = t107.*t124.*t125;
t159 = t158.*1.0e+2;
t160 = t145+t151;
t135 = t62+t129;
t136 = t64+t130;
t146 = -t144;
t155 = -t154;
t156 = t154.*1.0e+2;
t161 = t148+t152;
t162 = t49+t57+t143+t147;
t137 = exp(t135);
t138 = exp(t136);
t157 = -t156;
t163 = t93+t146+t150;
t164 = t76+t155+t158;
t139 = t137+1.0;
t165 = t77+t157+t159;
t140 = 1.0./t139;
t141 = t140.^2;
mt1 = [t43+t94.*t117.*(3.1e+1./2.0e+2)+th9.*(t126-2.3e+1./4.0e+2)+t137.*t140.*th10.*(t50-t126),t56,-t94.*t117+t137.*t140.*(t50-t126).*3.652e-1-t137.*t140.*th10.*(t50-t126).*3.652e-1,t94.*t117.*(-6.427915e-1)-t18.*t22.*th3.*x4,0.0,x1.*(th9.*(t142-t147)+t137.*t140.*t162.*th10+t137.*t140.*t163.*th10.*(t50-t126)-t138.*t141.*t163.*th10.*(t50-t126)),t43-t49.*x1+t51.*x1];
mt2 = [x1.*(t137.*t140.*t162.*3.652e-1-t137.*t140.*t162.*th10.*3.652e-1+t137.*t140.*t163.*(t50-t126).*3.652e-1-t138.*t141.*t163.*(t50-t126).*3.652e-1-t137.*t140.*t163.*th10.*(t50-t126).*3.652e-1+t138.*t141.*t163.*th10.*(t50-t126).*3.652e-1),0.0,0.0,-x1.*(t94.*t149.*(3.1e+1./2.0e+2)+th9.*(t154-t158)-t95.*t117.*t119.*(3.1e+1./2.0e+2)+t137.*t140.*t164.*th10+t137.*t140.*t165.*th10.*(t50-t126)-t138.*t141.*t165.*th10.*(t50-t126)),t76.*x1];
mt3 = [t43+x1.*(t94.*t149-t95.*t117.*t119-t137.*t140.*t164.*3.652e-1+t137.*t140.*t164.*th10.*3.652e-1-t137.*t140.*t165.*(t50-t126).*3.652e-1+t138.*t141.*t165.*(t50-t126).*3.652e-1+t137.*t140.*t165.*th10.*(t50-t126).*3.652e-1-t138.*t141.*t165.*th10.*(t50-t126).*3.652e-1),x1.*(t94.*t149.*6.427915e-1-t95.*t117.*t119.*6.427915e-1+t5.*t19.*t22.*th3.*x4),0.0];
mt4 = [-x1.*(t94.*t133.*(3.1e+1./2.0e+2)+t160.*th9+t79.*t86.*t95.*t117.*(3.1e+1./2.0e+2)-t137.*t140.*t160.*th10-t137.*t140.*t161.*th10.*(t50-t126)+t138.*t141.*t161.*th10.*(t50-t126)),0.0,x1.*(t94.*t133+t137.*t140.*t160.*3.652e-1+t79.*t86.*t95.*t117-t137.*t140.*t160.*th10.*3.652e-1+t137.*t140.*t161.*(t50-t126).*3.652e-1-t138.*t141.*t161.*(t50-t126).*3.652e-1-t137.*t140.*t161.*th10.*(t50-t126).*3.652e-1+t138.*t141.*t161.*th10.*(t50-t126).*3.652e-1)];
mt5 = [x1.*(t94.*t133.*6.427915e-1-t18.*t22.*th3+t79.*t86.*t95.*t117.*6.427915e-1+t18.*t23.*th3.*x4)-t6.*t8.*t16.*6.05e-4-1.8e+4,0.0,t6.*t9.*t16.*x1.*3.025e-4,t6.*t9.*t16.*(x2-4.5e+2).*3.025e-4,t6.*t9.*t16.*x3.*3.025e-4,t6.*t9.*t16.*x4.*6.05e-4,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5],5,5);
