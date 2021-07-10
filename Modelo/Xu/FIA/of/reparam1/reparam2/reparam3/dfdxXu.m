function dfdxSym = dfdxXu(t,in2,in3)
%DFDXXU
%    DFDXSYM = DFDXXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 20:22:36

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th8 = in3(5,:);
th9 = in3(6,:);
th10 = in3(7,:);
th11 = in3(8,:);
th14 = in3(9,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th5+x2;
t3 = th1.^2;
t4 = x2.^2;
t5 = th9-1.0;
t6 = 1.0./th8;
t7 = 1.0./th9;
t8 = 1.0./th11;
t9 = 1.0./x5;
t11 = x3./4.0;
t19 = th9.*(2.3e+1./4.0e+2);
t20 = x4+1.0e-4;
t23 = t.*(1.97e+2./1.0e+3);
t30 = th3.*9.002366863905325e-1;
t31 = th3.*7.20189349112426e+1;
t32 = x3+5.019e-1;
t10 = t9.^2;
t12 = t6.*x3;
t13 = 1.0./t2;
t15 = 1.0./t5;
t17 = t11+1.0;
t24 = exp(t23);
t28 = 1.0./t20;
t33 = -t31;
t34 = 1.0./t32;
t14 = t13.^2;
t16 = t15.^2;
t18 = t12+1.0;
t21 = 1.0./t17;
t29 = t28.^2;
t35 = t34.^2;
t37 = t34.*th2.*8.0e+1;
t44 = t34.*th2.*x3.*-8.0e+1;
t61 = t7.*t9.*t24.*5.4175e-4;
t22 = t21.^2;
t25 = 1.0./t18;
t40 = t37.*x3;
t41 = t35.*th2.*x3.*8.0e+1;
t50 = exp(t44);
t56 = t8.*t21.*t28.*th3;
t58 = t8.*t21.*t29.*th3.*x4;
t62 = -t61;
t63 = t21.*t28.*t30;
t64 = t21.*t28.*t31;
t66 = t21.*t29.*t30.*x4;
t68 = t21.*t29.*t31.*x4;
t69 = t21.*t28.*th3.*x4.*(-9.002366863905325e-1);
t70 = t21.*t29.*th3.*x4.*(-9.002366863905325e-1);
t71 = t21.*t29.*th3.*x4.*(-7.20189349112426e+1);
t26 = t25.^2;
t27 = t25.^3;
t36 = t13.*t25.*th1;
t39 = t14.*t25.*th1.*x2;
t45 = -t41;
t57 = t56.*x4;
t60 = -t58;
t65 = t63.*x4;
t67 = t64.*x4;
t73 = t34.*t50.*th2;
t75 = t35.*t50.*th2.*x3;
t83 = t30+t69;
t95 = t63+t70;
t96 = t64+t71;
t38 = t36.*x2;
t43 = -t39;
t46 = t36.*8.0e+1;
t48 = t39.*8.0e+1;
t51 = t6.*t13.*t26.*th1.*x2;
t59 = -t57;
t74 = t73.*x3;
t76 = -t73;
t77 = t37+t45;
t84 = t33+t67;
t93 = t56+t60;
t42 = -t38;
t47 = t38.*8.0e+1;
t49 = t38.*1.6e+2;
t53 = -t48;
t55 = t51.*8.0e+1;
t72 = t19+t59;
t86 = exp(t84);
t94 = t50.*t77;
t103 = t74.*t77;
t52 = -t47;
t82 = t46+t53;
t87 = t15.*t72.*8.0e+1;
t97 = t50+t86;
t104 = t22.*t28.*t86.*th3.*x4.*1.800473372781065e+1;
t105 = t22.*t28.*t86.*th3.*x4.*2.250591715976331e-1;
t111 = t83.*t86;
t114 = t86.*t95;
t54 = exp(t52);
t89 = -t87;
t98 = 1.0./t97;
t106 = -t105;
t117 = t83.*t104;
t118 = t74+t111;
t121 = t94+t104;
t124 = t96.*t111;
t78 = t36.*t54;
t79 = t38.*t54;
t80 = t39.*t54;
t85 = t51.*t54;
t90 = t54.*t55;
t91 = exp(t89);
t92 = t3.*t4.*t6.*t14.*t27.*t54.*8.0e+1;
t99 = t98.^2;
t125 = -t124;
t157 = t75+t76+t103+t106+t117;
t81 = -t78;
t88 = -t85;
t100 = t54+t91;
t107 = t8.*t15.*t22.*t28.*t91.*th3.*x4.*2.0e+1;
t108 = (t8.*t15.*t22.*t28.*t91.*th3.*x4)./4.0;
t110 = t79.*t82;
t112 = t15.*t72.*t91;
t113 = t15.*t91.*t93;
t115 = t8.*t16.*t22.*t28.*t72.*t91.*th3.*x4.*2.0e+1;
t122 = t16.*t72.*t91.*t93.*8.0e+1;
t143 = t114+t125;
t101 = 1.0./t100;
t109 = -t107;
t116 = -t115;
t120 = t79+t112;
t123 = -t122;
t126 = t80+t81+t110;
t102 = t101.^2;
t119 = t90+t109;
t127 = t101.*t120;
t140 = t101.*t126;
t142 = t113+t123;
t148 = t88+t92+t108+t116;
t158 = -t101.*(t85-t92-t108+t115);
t159 = t101.*(t85-t92-t108+t115).*-8.0e+1;
t128 = t127.*8.0e+1;
t129 = t127.*1.6e+2;
t132 = t42+t127;
t141 = t140.*8.0e+1;
t144 = t54.*t82.*t102.*t120;
t149 = t102.*t113.*t120.*8.0e+1;
t150 = t102.*t113.*t120.*6.4e+3;
t151 = t101.*t142;
t153 = t102.*t119.*t120;
t130 = -t128;
t131 = -t129;
t145 = -t144;
t146 = t144.*8.0e+1;
t152 = t151.*8.0e+1;
t154 = -t153;
t155 = t153.*8.0e+1;
t162 = t149+t151;
t133 = t47+t130;
t134 = t49+t131;
t147 = -t146;
t156 = -t155;
t160 = t36+t43+t140+t145;
t163 = t150+t152;
t166 = t51+t154+t158;
t135 = exp(t133);
t136 = exp(t134);
t161 = t82+t141+t147;
t167 = t55+t156+t159;
t137 = t135+1.0;
t138 = 1.0./t137;
t139 = t138.^2;
t164 = -t135.*t138.*th10.*(t38-t127);
t165 = t135.*t138.*th10.*(t38-t127);
t168 = t135.*t138.*t160.*th10;
t170 = t135.*t138.*t162.*th10;
t172 = -t135.*t138.*th10.*(-t51+t153+t101.*(t85-t92-t108+t115));
t169 = -t168;
t171 = -t170;
t173 = t161.*t164;
t174 = -t136.*t139.*t161.*th10.*(t38-t127);
t175 = t136.*t139.*t161.*th10.*(t38-t127);
t176 = t163.*t164;
t177 = -t136.*t139.*t163.*th10.*(t38-t127);
t178 = t136.*t139.*t163.*th10.*(t38-t127);
t179 = t165.*(-t55+t155+t101.*(t85-t92-t108+t115).*8.0e+1);
t180 = t136.*t139.*th10.*(t38-t127).*(-t55+t155+t101.*(t85-t92-t108+t115).*8.0e+1);
mt1 = [t62+t165+t98.*t118.*(3.1e+1./2.0e+2)+th9.*(t127-2.3e+1./4.0e+2),t42,-t98.*t118+th14.*(t164+t135.*t138.*(t38-t127)),t98.*t118.*(-6.427915e-1)-t21.*t28.*th3.*x4,0.0,-x1.*(t169+t173+t175+th9.*(t140+t145)),t62-t36.*x1+t39.*x1,th14.*x1.*(t169+t173+t175+t135.*t138.*t160+t135.*t138.*t161.*(t38-t127)-t136.*t139.*t161.*(t38-t127)),0.0,0.0,-x1.*(t172+t180+t164.*(-t55+t155+t101.*(t85-t92-t108+t115).*8.0e+1)+t98.*t157.*(3.1e+1./2.0e+2)+th9.*(t153+t101.*(t85-t92-t108+t115))-t99.*t118.*t121.*(3.1e+1./2.0e+2)),t51.*x1];
mt2 = [t62+x1.*(t98.*t157+th14.*(t172+t180+t164.*(-t55+t155+t101.*(t85-t92-t108+t115).*8.0e+1)+t135.*t138.*(-t51+t153+t101.*(t85-t92-t108+t115))+t135.*t138.*(t38-t127).*(-t55+t155+t101.*(t85-t92-t108+t115).*8.0e+1)-t136.*t139.*(t38-t127).*(-t55+t155+t101.*(t85-t92-t108+t115).*8.0e+1))-t99.*t118.*t121),x1.*(t98.*t157.*6.427915e-1-t99.*t118.*t121.*6.427915e-1+(t22.*t28.*th3.*x4)./4.0),0.0];
mt3 = [-x1.*(t171+t176+t178+t98.*t143.*(3.1e+1./2.0e+2)+t162.*th9+t86.*t96.*t99.*t118.*(3.1e+1./2.0e+2)),0.0,x1.*(th14.*(t171+t176+t178+t135.*t138.*t162+t135.*t138.*t163.*(t38-t127)-t136.*t139.*t163.*(t38-t127))+t98.*t143+t86.*t96.*t99.*t118),x1.*(t98.*t143.*6.427915e-1-t21.*t28.*th3+t86.*t96.*t99.*t118.*6.427915e-1+t21.*t29.*th3.*x4)-t7.*t9.*t24.*1.0835e-3-1.8e+4,0.0,t7.*t10.*t24.*x1.*5.4175e-4,t7.*t10.*t24.*(x2-4.5e+2).*5.4175e-4];
mt4 = [t7.*t10.*t24.*x3.*5.4175e-4,t7.*t10.*t24.*x4.*1.0835e-3,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4],5,5);
