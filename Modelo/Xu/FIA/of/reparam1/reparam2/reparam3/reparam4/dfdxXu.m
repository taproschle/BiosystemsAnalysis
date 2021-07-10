function dfdxSym = dfdxXu(t,in2,in3)
%DFDXXU
%    DFDXSYM = DFDXXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 20:32:14

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th9 = in3(4,:);
th10 = in3(5,:);
th11 = in3(6,:);
th14 = in3(7,:);
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
t35 = t30.*th2.*8.0e+1;
t38 = t30.*th2.*x3.*-8.0e+1;
t46 = t5.*t7.*t18.*5.4175e-4;
t16 = t15.^2;
t33 = t32.^2;
t34 = t32.^3;
t36 = t35.*x3;
t37 = t31.*th2.*x3.*8.0e+1;
t40 = exp(t38);
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
t39 = -t37;
t42 = t41.*x4;
t45 = -t43;
t49 = t48.*x2;
t53 = t51.*x4;
t55 = t52.*x4;
t60 = -t50;
t62 = t48.*8.0e+1;
t64 = t50.*8.0e+1;
t70 = t30.*t40.*th2;
t72 = t31.*t40.*th2.*x3;
t75 = t28.*t33.*th1.*x2.*1.344483583855441e-1;
t76 = t28.*t33.*th1.*x2.*1.075586867084353e+1;
t77 = t22+t57;
t85 = t51+t58;
t86 = t52+t61;
t44 = -t42;
t59 = -t49;
t63 = t49.*8.0e+1;
t65 = t49.*1.6e+2;
t67 = -t64;
t71 = t70.*x3;
t73 = -t70;
t74 = t35+t39;
t78 = t26+t55;
t83 = t41+t45;
t66 = -t63;
t69 = t13+t44;
t79 = exp(t78);
t84 = t40.*t74;
t92 = t62+t67;
t96 = t71.*t74;
t68 = exp(t66);
t80 = t10.*t69.*8.0e+1;
t91 = t40+t79;
t99 = t16.*t20.*t79.*th3.*x4.*1.800473372781065e+1;
t100 = t16.*t20.*t79.*th3.*x4.*2.250591715976331e-1;
t109 = t77.*t79;
t112 = t79.*t85;
t81 = -t80;
t87 = t48.*t68;
t88 = t49.*t68;
t89 = t50.*t68;
t93 = 1.0./t91;
t95 = t68.*t75;
t97 = t28.*t33.*t68.*th1.*x2.*(-1.344483583855441e-1);
t98 = t68.*t76;
t101 = -t100;
t102 = t2.*t3.*t29.*t34.*t68.*1.075586867084353e+1;
t115 = t77.*t99;
t116 = t71+t109;
t118 = t84+t99;
t121 = t86.*t109;
t82 = exp(t81);
t90 = -t87;
t94 = t93.^2;
t117 = t88.*t92;
t122 = -t121;
t145 = t72+t73+t96+t101+t115;
t103 = t6.*t10.*t16.*t20.*t82.*th3.*x4.*2.0e+1;
t104 = (t6.*t10.*t16.*t20.*t82.*th3.*x4)./4.0;
t106 = t68+t82;
t110 = t10.*t69.*t82;
t111 = t10.*t82.*t83;
t113 = t6.*t11.*t16.*t20.*t69.*t82.*th3.*x4.*2.0e+1;
t119 = t11.*t69.*t82.*t83.*8.0e+1;
t131 = t112+t122;
t132 = t89+t90+t117;
t105 = -t103;
t107 = 1.0./t106;
t114 = -t113;
t120 = -t119;
t123 = t88+t110;
t108 = t107.^2;
t124 = t98+t105;
t125 = t107.*t123;
t130 = t111+t120;
t148 = t107.*t132;
t152 = t97+t102+t104+t114;
t126 = t125.*8.0e+1;
t127 = t125.*1.6e+2;
t133 = t59+t125;
t141 = t107.*t130;
t142 = t68.*t92.*t108.*t123;
t149 = t108.*t111.*t123.*8.0e+1;
t150 = t148.*8.0e+1;
t151 = t108.*t111.*t123.*6.4e+3;
t153 = t108.*t123.*t124;
t157 = t107.*t152;
t128 = -t126;
t129 = -t127;
t143 = t141.*8.0e+1;
t144 = -t142;
t146 = t142.*8.0e+1;
t154 = -t153;
t155 = t153.*8.0e+1;
t158 = t157.*8.0e+1;
t159 = t141+t149;
t134 = t63+t128;
t135 = t65+t129;
t147 = -t146;
t156 = -t155;
t160 = t143+t151;
t161 = t48+t60+t144+t148;
t165 = t75+t154+t157;
t136 = exp(t134);
t137 = exp(t135);
t162 = t92+t147+t150;
t166 = t76+t156+t158;
t138 = t136+1.0;
t139 = 1.0./t138;
t140 = t139.^2;
t163 = -t136.*t139.*th10.*(t49-t125);
t164 = t136.*t139.*th10.*(t49-t125);
t167 = t136.*t139.*t159.*th10;
t169 = t136.*t139.*t161.*th10;
t171 = t136.*t139.*t165.*th10;
t168 = -t167;
t170 = -t169;
t172 = -t171;
t173 = t160.*t163;
t174 = -t137.*t140.*t160.*th10.*(t49-t125);
t175 = t137.*t140.*t160.*th10.*(t49-t125);
t176 = t162.*t163;
t177 = -t137.*t140.*t162.*th10.*(t49-t125);
t178 = t137.*t140.*t162.*th10.*(t49-t125);
t179 = t163.*t166;
t180 = -t137.*t140.*t166.*th10.*(t49-t125);
t181 = t137.*t140.*t166.*th10.*(t49-t125);
mt1 = [t47+t164+t93.*t116.*(3.1e+1./2.0e+2)+th9.*(t125-2.3e+1./4.0e+2),t59,-t93.*t116+th14.*(t163+t136.*t139.*(t49-t125)),t93.*t116.*(-6.427915e-1)-t15.*t20.*th3.*x4,0.0,x1.*(t169+t177+t162.*t164+th9.*(t142-t148)),t47-t48.*x1+t50.*x1,th14.*x1.*(t170+t176+t178+t136.*t139.*t161+t136.*t139.*t162.*(t49-t125)-t137.*t140.*t162.*(t49-t125)),0.0,0.0,-x1.*(t171+t180+t93.*t145.*(3.1e+1./2.0e+2)+t164.*t166+th9.*(t153-t157)-t94.*t116.*t118.*(3.1e+1./2.0e+2)),t75.*x1];
mt2 = [t47-x1.*(th14.*(t172+t179+t181+t136.*t139.*t165+t136.*t139.*t166.*(t49-t125)-t137.*t140.*t166.*(t49-t125))-t93.*t145+t94.*t116.*t118),x1.*(t93.*t145.*6.427915e-1-t94.*t116.*t118.*6.427915e-1+(t16.*t20.*th3.*x4)./4.0),0.0,-x1.*(t168+t173+t175+t93.*t131.*(3.1e+1./2.0e+2)+t159.*th9+t79.*t86.*t94.*t116.*(3.1e+1./2.0e+2)),0.0,x1.*(th14.*(t168+t173+t175+t136.*t139.*t159+t136.*t139.*t160.*(t49-t125)-t137.*t140.*t160.*(t49-t125))+t93.*t131+t79.*t86.*t94.*t116)];
mt3 = [x1.*(t93.*t131.*6.427915e-1-t15.*t20.*th3+t79.*t86.*t94.*t116.*6.427915e-1+t15.*t21.*th3.*x4)-t5.*t7.*t18.*1.0835e-3-1.8e+4,0.0,t5.*t8.*t18.*x1.*5.4175e-4,t5.*t8.*t18.*(x2-4.5e+2).*5.4175e-4,t5.*t8.*t18.*x3.*5.4175e-4,t5.*t8.*t18.*x4.*1.0835e-3,0.0];
dfdxSym = reshape([mt1,mt2,mt3],5,5);