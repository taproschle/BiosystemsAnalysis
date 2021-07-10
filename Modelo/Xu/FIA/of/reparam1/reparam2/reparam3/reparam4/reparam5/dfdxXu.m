function dfdxSym = dfdxXu(t,in2,in3)
%DFDXXU
%    DFDXSYM = DFDXXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 21:52:07

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
t2 = th1.^2;
t3 = x2.^2;
t4 = th9-1.0;
t5 = 1.0./th9;
t6 = 1.0./th11;
t7 = 1.0./x5;
t9 = x3./4.0;
t13 = th9.*(2.3e+1./4.0e+2);
t14 = x4+1.0e-4;
t17 = t.*(1.97e+2./1.0e+3);
t19 = x2+1.034e-1;
t22 = th3.*9.002366863905325e-1;
t23 = th3.*7.20189349112426e+1;
t24 = x3.*1.344483583855441e-1;
t25 = x3+5.019e-1;
t8 = t7.^2;
t10 = 1.0./t4;
t12 = t9+1.0;
t18 = exp(t17);
t20 = 1.0./t14;
t26 = -t23;
t27 = t24+1.0;
t28 = 1.0./t19;
t30 = 1.0./t25;
t11 = t10.^2;
t15 = 1.0./t12;
t21 = t20.^2;
t29 = t28.^2;
t31 = t30.^2;
t32 = 1.0./t27;
t35 = t30.*9.792;
t39 = t30.*x3.*(-9.792);
t46 = t5.*t7.*t18.*5.4175e-4;
t16 = t15.^2;
t33 = t32.^2;
t34 = t32.^3;
t36 = t35.*x3;
t37 = t31.*x3.*9.792;
t38 = -t35;
t40 = exp(t39);
t41 = t6.*t15.*t20.*th3;
t43 = t6.*t15.*t21.*th3.*x4;
t47 = -t46;
t48 = t28.*t32.*th1;
t50 = t29.*t32.*th1.*x2;
t51 = t15.*t20.*t22;
t52 = t15.*t20.*t23;
t54 = t15.*t21.*t22.*x4;
t56 = t15.*t21.*t23.*x4;
t57 = t15.*t20.*th3.*x4.*(-9.002366863905325e-1);
t58 = t15.*t21.*th3.*x4.*(-9.002366863905325e-1);
t61 = t15.*t21.*th3.*x4.*(-7.20189349112426e+1);
t42 = t41.*x4;
t45 = -t43;
t49 = t48.*x2;
t53 = t51.*x4;
t55 = t52.*x4;
t60 = -t50;
t62 = t48.*8.0e+1;
t64 = t50.*8.0e+1;
t70 = t28.*t33.*th1.*x2.*1.344483583855441e-1;
t71 = t28.*t33.*th1.*x2.*1.075586867084353e+1;
t72 = t22+t57;
t75 = t37+t38;
t76 = t30.*t40.*1.224e-1;
t79 = t31.*t40.*x3.*1.224e-1;
t84 = t51+t58;
t85 = t52+t61;
t92 = -t40.*(t35-t37);
t110 = t30.*t40.*x3.*(t35-t37).*(-1.224e-1);
t44 = -t42;
t59 = -t49;
t63 = t49.*8.0e+1;
t65 = t49.*1.6e+2;
t67 = -t64;
t73 = t26+t55;
t78 = t76.*x3;
t82 = -t79;
t83 = t41+t45;
t66 = -t63;
t69 = t13+t44;
t74 = exp(t73);
t90 = t62+t67;
t68 = exp(t66);
t77 = t10.*t69.*8.0e+1;
t91 = t40+t74;
t98 = t16.*t20.*t74.*th3.*x4.*1.800473372781065e+1;
t99 = t16.*t20.*t74.*th3.*x4.*2.250591715976331e-1;
t108 = t72.*t74;
t112 = t74.*t84;
t80 = -t77;
t86 = t48.*t68;
t87 = t49.*t68;
t88 = t50.*t68;
t93 = 1.0./t91;
t95 = t68.*t70;
t96 = t28.*t33.*t68.*th1.*x2.*(-1.344483583855441e-1);
t97 = t68.*t71;
t100 = -t98;
t101 = t2.*t3.*t29.*t34.*t68.*1.075586867084353e+1;
t115 = t72.*t98;
t116 = t16.*t20.*t108.*th3.*x4.*(-1.800473372781065e+1);
t118 = t78+t108;
t122 = t85.*t108;
t81 = exp(t80);
t89 = -t86;
t94 = t93.^2;
t117 = t87.*t90;
t121 = t92+t100;
t123 = -t122;
t153 = t76+t82+t99+t110+t116;
t102 = t6.*t10.*t16.*t20.*t81.*th3.*x4.*2.0e+1;
t103 = (t6.*t10.*t16.*t20.*t81.*th3.*x4)./4.0;
t105 = t68+t81;
t109 = t10.*t69.*t81;
t111 = t10.*t81.*t83;
t113 = t6.*t11.*t16.*t20.*t69.*t81.*th3.*x4.*2.0e+1;
t119 = t11.*t69.*t81.*t83.*8.0e+1;
t132 = t112+t123;
t133 = t88+t89+t117;
t104 = -t102;
t106 = 1.0./t105;
t114 = -t113;
t120 = -t119;
t124 = t87+t109;
t107 = t106.^2;
t125 = t97+t104;
t126 = t106.*t124;
t131 = t111+t120;
t148 = t106.*t133;
t152 = t96+t101+t103+t114;
t127 = t126.*8.0e+1;
t128 = t126.*1.6e+2;
t134 = t59+t126;
t142 = t106.*t131;
t143 = t68.*t90.*t107.*t124;
t149 = t107.*t111.*t124.*8.0e+1;
t150 = t148.*8.0e+1;
t151 = t107.*t111.*t124.*6.4e+3;
t154 = t107.*t124.*t125;
t158 = t106.*t152;
t129 = -t127;
t130 = -t128;
t144 = t142.*8.0e+1;
t145 = -t143;
t146 = t143.*8.0e+1;
t155 = -t154;
t156 = t154.*8.0e+1;
t159 = t158.*8.0e+1;
t160 = t142+t149;
t135 = t63+t129;
t136 = t65+t130;
t147 = -t146;
t157 = -t156;
t161 = t144+t151;
t162 = t48+t60+t145+t148;
t164 = t70+t155+t158;
t137 = exp(t135);
t138 = exp(t136);
t163 = t90+t147+t150;
t165 = t71+t157+t159;
t139 = t137+1.0;
t140 = 1.0./t139;
t141 = t140.^2;
mt1 = [t47+t93.*t118.*(3.1e+1./2.0e+2)+th9.*(t126-2.3e+1./4.0e+2)+t137.*t140.*th10.*(t49-t126),t59,-t93.*t118+t137.*t140.*(t49-t126).*3.652e-1-t137.*t140.*th10.*(t49-t126).*3.652e-1,t93.*t118.*(-6.427915e-1)-t15.*t20.*th3.*x4,0.0,x1.*(th9.*(t143-t148)+t137.*t140.*t162.*th10+t137.*t140.*t163.*th10.*(t49-t126)-t138.*t141.*t163.*th10.*(t49-t126)),t47-t48.*x1+t50.*x1];
mt2 = [x1.*(t137.*t140.*t162.*3.652e-1-t137.*t140.*t162.*th10.*3.652e-1+t137.*t140.*t163.*(t49-t126).*3.652e-1-t138.*t141.*t163.*(t49-t126).*3.652e-1-t137.*t140.*t163.*th10.*(t49-t126).*3.652e-1+t138.*t141.*t163.*th10.*(t49-t126).*3.652e-1),0.0,0.0];
mt3 = [-x1.*(t93.*(-t76+t79-t99+t115+t78.*(t35-t37)).*(3.1e+1./2.0e+2)+th9.*(t154-t158)-t94.*t118.*(t98+t40.*(t35-t37)).*(3.1e+1./2.0e+2)+t137.*t140.*t164.*th10+t137.*t140.*t165.*th10.*(t49-t126)-t138.*t141.*t165.*th10.*(t49-t126)),t70.*x1];
mt4 = [t47+x1.*(t93.*(-t76+t79-t99+t115+t78.*(t35-t37))-t94.*t118.*(t98+t40.*(t35-t37))-t137.*t140.*t164.*3.652e-1+t137.*t140.*t164.*th10.*3.652e-1-t137.*t140.*t165.*(t49-t126).*3.652e-1+t138.*t141.*t165.*(t49-t126).*3.652e-1+t137.*t140.*t165.*th10.*(t49-t126).*3.652e-1-t138.*t141.*t165.*th10.*(t49-t126).*3.652e-1)];
mt5 = [x1.*(t93.*(-t76+t79-t99+t115+t78.*(t35-t37)).*6.427915e-1-t94.*t118.*(t98+t40.*(t35-t37)).*6.427915e-1+(t16.*t20.*th3.*x4)./4.0),0.0,-x1.*(t93.*t132.*(3.1e+1./2.0e+2)+t160.*th9+t74.*t85.*t94.*t118.*(3.1e+1./2.0e+2)-t137.*t140.*t160.*th10-t137.*t140.*t161.*th10.*(t49-t126)+t138.*t141.*t161.*th10.*(t49-t126)),0.0];
mt6 = [x1.*(t93.*t132+t137.*t140.*t160.*3.652e-1+t74.*t85.*t94.*t118-t137.*t140.*t160.*th10.*3.652e-1+t137.*t140.*t161.*(t49-t126).*3.652e-1-t138.*t141.*t161.*(t49-t126).*3.652e-1-t137.*t140.*t161.*th10.*(t49-t126).*3.652e-1+t138.*t141.*t161.*th10.*(t49-t126).*3.652e-1),x1.*(t93.*t132.*6.427915e-1-t15.*t20.*th3+t74.*t85.*t94.*t118.*6.427915e-1+t15.*t21.*th3.*x4)-t5.*t7.*t18.*1.0835e-3-1.8e+4,0.0];
mt7 = [t5.*t8.*t18.*x1.*5.4175e-4,t5.*t8.*t18.*(x2-4.5e+2).*5.4175e-4,t5.*t8.*t18.*x3.*5.4175e-4,t5.*t8.*t18.*x4.*1.0835e-3,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7],5,5);
