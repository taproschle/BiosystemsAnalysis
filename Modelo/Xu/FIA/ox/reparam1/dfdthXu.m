function dfdthSym = dfdthXu(t,in2,in3)
%DFDTHXU
%    DFDTHSYM = DFDTHXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 16:21:22

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th7 = in3(7,:);
th8 = in3(8,:);
th9 = in3(9,:);
th10 = in3(10,:);
th11 = in3(11,:);
th12 = in3(12,:);
th14 = in3(13,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th4.*th9;
t3 = th5+x2;
t4 = th6+x3;
t5 = th1.^2;
t6 = th2.^2;
t7 = th12.^2;
t8 = x2.^2;
t9 = x3.^2;
t10 = th9-1.0;
t11 = 1.0./th7;
t13 = 1.0./th8;
t15 = 1.0./th9.^2;
t16 = 1.0./th11;
t18 = 1.0./x5;
t31 = t.*(1.1e+1./1.0e+2);
t33 = x4+1.0e-4;
t39 = th3.*1.183431952662722e+2;
t12 = t11.^2;
t14 = t13.^2;
t17 = t16.^2;
t19 = t11.*x3;
t20 = t13.*x3;
t21 = 1.0./t3;
t24 = 1.0./t4;
t27 = 1.0./t10;
t32 = exp(t31);
t41 = 1.0./t33;
t22 = t21.^2;
t23 = t21.^3;
t25 = t24.^2;
t26 = t24.^3;
t28 = t27.^2;
t29 = t19+1.0;
t30 = t20+1.0;
t40 = t27.*th4.*1.0e+2;
t42 = t24.*th2.*x3.*1.0e+2;
t34 = 1.0./t29;
t36 = 1.0./t30;
t43 = -t42;
t35 = t34.^2;
t37 = t36.^2;
t38 = t36.^3;
t44 = exp(t43);
t45 = t21.*t36.*x2;
t47 = t22.*t36.*th1.*x2;
t62 = t34.*t41.*x4;
t46 = t45.*th1;
t48 = -t45;
t50 = -t47;
t51 = t24.*t44.*x3;
t52 = t45.*1.0e+2;
t54 = t25.*t44.*th2.*x3;
t56 = t47.*1.0e+2;
t63 = t14.*t21.*t37.*th1.*x2.*x3;
t65 = t62.*th3;
t67 = t9.*t25.*t44.*th2.*1.0e+2;
t68 = t62-1.0;
t71 = t6.*t9.*t26.*t44.*1.0e+2;
t75 = t39.*t62;
t49 = -t46;
t53 = t51.*th2;
t55 = t46.*1.0e+2;
t57 = t46.*2.0e+2;
t58 = -t52;
t60 = -t56;
t61 = -t54;
t66 = t63.*1.0e+2;
t69 = -t65;
t70 = -t67;
t73 = t16.*t65;
t77 = t65.*(-1.183431952662722e+2);
t78 = th12.*(t65-th3).*(-1.183431952662722e+2);
t79 = th12.*(t65-th3).*1.183431952662722e+2;
t59 = -t55;
t72 = t69+th3;
t74 = t16.*t69;
t80 = exp(t79);
t86 = t39+t77;
t95 = t51+t70;
t97 = t61+t71;
t64 = exp(t59);
t76 = t2+t74;
t101 = t44+t80;
t109 = t68.*t80.*th12.*(2.0e+2./1.69e+2);
t110 = t80.*(t65-th3).*(-2.0e+2./1.69e+2);
t112 = t12.*t35.*t41.*t80.*th3.*th12.*x3.*x4.*(2.0e+2./1.69e+2);
t124 = t7.*t68.*t80.*(t65-th3).*(-1.400511186583103e+2);
t125 = t7.*t68.*t80.*(t65-th3).*1.400511186583103e+2;
t126 = t7.*t12.*t35.*t41.*t80.*th3.*x3.*x4.*(t65-th3).*(-1.400511186583103e+2);
t127 = t7.*t12.*t35.*t41.*t80.*th3.*x3.*x4.*(t65-th3).*1.400511186583103e+2;
t135 = t80.*t86.*th12.*(t65-th3).*(2.0e+2./1.69e+2);
t81 = t45.*t64;
t82 = t46.*t64;
t83 = t47.*t64;
t84 = t27.*t76.*1.0e+2;
t85 = t28.*t76.*1.0e+2;
t90 = t63.*t64;
t91 = t8.*t22.*t37.*t64.*th1.*1.0e+2;
t94 = t5.*t8.*t23.*t37.*t64.*1.0e+2;
t98 = t5.*t8.*t14.*t22.*t38.*t64.*x3.*1.0e+2;
t103 = 1.0./t101;
t111 = t110.*th12;
t145 = t109+t125;
t146 = t112+t127;
t148 = t110+t135;
t87 = -t84;
t88 = -t85;
t92 = -t91;
t93 = -t90;
t96 = -t94;
t104 = t103.^2;
t122 = t53+t111;
t133 = t86.*t111;
t89 = exp(t87);
t105 = t40+t88;
t118 = t81+t92;
t121 = t83+t96;
t123 = t93+t98;
t99 = t27.*t89.*th4;
t100 = t27.*t89.*th9;
t106 = t64+t89;
t113 = t16.*t27.*t62.*t89;
t114 = t17.*t27.*t65.*t89;
t115 = t12.*t16.*t27.*t35.*t41.*t89.*th3.*x3.*x4;
t116 = t27.*t76.*t89;
t117 = t28.*t76.*t89;
t119 = t85.*t89.*th9;
t128 = t16.*t62.*t85.*t89;
t130 = t17.*t65.*t85.*t89;
t132 = t12.*t16.*t35.*t41.*t85.*t89.*th3.*x3.*x4;
t102 = -t99;
t107 = 1.0./t106;
t120 = t117.*th9.*-1.0e+2;
t129 = t16.*t62.*t117.*-1.0e+2;
t131 = t17.*t65.*t117.*-1.0e+2;
t134 = t12.*t16.*t35.*t41.*t117.*th3.*x3.*x4.*-1.0e+2;
t136 = t82+t116;
t137 = t105.*t116;
t108 = t107.^2;
t138 = t100+t120;
t139 = t107.*t118;
t140 = t107.*t121;
t143 = -t107.*(t90-t98);
t144 = t107.*(t90-t98).*-1.0e+2;
t147 = t113+t129;
t149 = t114+t131;
t150 = t107.*t136;
t155 = t115+t134;
t186 = t102+t117+t137;
t141 = t139.*1.0e+2;
t142 = t140.*1.0e+2;
t151 = t150.*1.0e+2;
t152 = t150.*2.0e+2;
t156 = t107.*t138;
t158 = t49+t150;
t166 = t52.*t64.*t108.*t136;
t167 = t56.*t64.*t108.*t136;
t168 = t81.*t108.*t136.*1.0e+4;
t169 = t83.*t108.*t136.*1.0e+4;
t170 = t64.*t66.*t108.*t136;
t171 = t90.*t108.*t136.*-1.0e+2;
t172 = t90.*t108.*t136.*1.0e+4;
t174 = t100.*t108.*t136.*1.0e+2;
t175 = t100.*t108.*t136.*1.0e+4;
t176 = t107.*t147;
t177 = t107.*t149;
t180 = t107.*t155;
t182 = t108.*t113.*t136.*1.0e+2;
t183 = t108.*t114.*t136.*1.0e+2;
t184 = t108.*t113.*t136.*1.0e+4;
t185 = t108.*t114.*t136.*1.0e+4;
t187 = t108.*t115.*t136.*1.0e+2;
t188 = t108.*t115.*t136.*1.0e+4;
t189 = t89.*t105.*t108.*t136;
t193 = t107.*t186;
t153 = -t151;
t154 = -t152;
t157 = t156.*1.0e+2;
t173 = -t172;
t178 = t176.*1.0e+2;
t179 = t177.*1.0e+2;
t181 = t180.*1.0e+2;
t190 = -t189;
t191 = t189.*1.0e+2;
t194 = t193.*1.0e+2;
t195 = t48+t139+t166;
t196 = t50+t140+t167;
t197 = t58+t141+t168;
t198 = t60+t142+t169;
t199 = t63+t143+t171;
t201 = t156+t174;
t203 = t176+t182;
t204 = t177+t183;
t207 = t180+t187;
t159 = t55+t153;
t160 = t57+t154;
t192 = -t191;
t200 = t66+t144+t173;
t202 = t157+t175;
t205 = t178+t184;
t206 = t179+t185;
t208 = t181+t188;
t209 = t190+t193;
t161 = exp(t159);
t162 = exp(t160);
t210 = t192+t194;
t163 = t161+1.0;
t164 = 1.0./t163;
t165 = t164.^2;
t211 = t161.*t164.*t195.*th10;
t213 = t161.*t164.*t196.*th10;
t215 = -t161.*t164.*th10.*(-t63+t170+t107.*(t90-t98));
t216 = t161.*t164.*th10.*(-t63+t170+t107.*(t90-t98));
t217 = t161.*t164.*t201.*th10;
t219 = t161.*t164.*t203.*th10;
t221 = t161.*t164.*t204.*th10;
t223 = t161.*t164.*t207.*th10;
t225 = -t161.*t164.*t197.*th10.*(t46-t150);
t228 = -t161.*t164.*t198.*th10.*(t46-t150);
t231 = -t161.*t164.*th10.*(t189-t193);
t232 = t161.*t164.*th10.*(t46-t150).*(-t66+t172+t107.*(t90-t98).*1.0e+2);
t235 = -t161.*t164.*t202.*th10.*(t46-t150);
t238 = -t161.*t164.*t205.*th10.*(t46-t150);
t241 = -t161.*t164.*t206.*th10.*(t46-t150);
t244 = -t161.*t164.*t208.*th10.*(t46-t150);
t247 = t161.*t164.*th10.*(t46-t150).*(t191-t194);
t212 = -t211;
t214 = -t213;
t218 = -t217;
t220 = -t219;
t222 = -t221;
t224 = -t223;
t226 = -t162.*t165.*t197.*th10.*(t46-t150);
t227 = t162.*t165.*t197.*th10.*(t46-t150);
t229 = -t162.*t165.*t198.*th10.*(t46-t150);
t230 = t162.*t165.*t198.*th10.*(t46-t150);
t233 = t162.*t165.*th10.*(t46-t150).*(-t66+t172+t107.*(t90-t98).*1.0e+2);
t236 = -t162.*t165.*t202.*th10.*(t46-t150);
t237 = t162.*t165.*t202.*th10.*(t46-t150);
t239 = -t162.*t165.*t205.*th10.*(t46-t150);
t240 = t162.*t165.*t205.*th10.*(t46-t150);
t242 = -t162.*t165.*t206.*th10.*(t46-t150);
t243 = t162.*t165.*t206.*th10.*(t46-t150);
t245 = -t162.*t165.*t208.*th10.*(t46-t150);
t246 = t162.*t165.*t208.*th10.*(t46-t150);
t248 = t162.*t165.*th10.*(t46-t150).*(t191-t194);
t234 = -t233;
mt1 = [-x1.*(t211+t226-th9.*(t139+t166)+t161.*t164.*t197.*th10.*(t46-t150)),t48.*x1,th14.*x1.*(t211+t226-t161.*t164.*t195-t161.*t164.*t197.*(t46-t150)+t162.*t165.*t197.*(t46-t150)+t161.*t164.*t197.*th10.*(t46-t150)),0.0,0.0,x1.*(t95.*t103.*(3.1e+1./2.0e+2)+t51.*t104.*t122.*(3.1e+1./2.0)),0.0,-x1.*(t95.*t103+t51.*t104.*t122.*1.0e+2),-x1.*(t95.*t103.*th12.*(1.69e+2./2.0e+2)+t51.*t104.*t122.*th12.*(1.69e+2./2.0)),0.0,-x1.*(t220+t238+t240+t103.*t145.*(3.1e+1./2.0e+2)+t203.*th9+t68.*t80.*t104.*t122.*th12.*1.834319526627219e+1),0.0];
mt2 = [x1.*(t103.*t145-th14.*(t219+t239-t161.*t164.*t203-t161.*t164.*t205.*(t46-t150)+t162.*t165.*t205.*(t46-t150)+t161.*t164.*t205.*th10.*(t46-t150))+t68.*t80.*t104.*t122.*th12.*1.183431952662722e+2),x1.*(-t62+t103.*t145.*th12.*(1.69e+2./2.0e+2)+t7.*t68.*t80.*t104.*t122.*1.0e+2),0.0,-x1.*(t217+t236-th9.*(t201-1.0)+t161.*t164.*t202.*th10.*(t46-t150)),0.0,th14.*x1.*(t217+t236-t161.*t164.*t201-t161.*t164.*t202.*(t46-t150)+t162.*t165.*t202.*(t46-t150)+t161.*t164.*t202.*th10.*(t46-t150)),0.0,0.0,x1.*(t213+t229-th9.*(t140+t167)+t161.*t164.*t198.*th10.*(t46-t150)),t47.*x1];
mt3 = [-th14.*x1.*(t213+t229-t161.*t164.*t196-t161.*t164.*t198.*(t46-t150)+t162.*t165.*t198.*(t46-t150)+t161.*t164.*t198.*th10.*(t46-t150)),0.0,0.0,-x1.*(t103.*(t54-t71).*(3.1e+1./2.0e+2)+t54.*t104.*t122.*(3.1e+1./2.0)),0.0,x1.*(t103.*(t54-t71)+t54.*t104.*t122.*1.0e+2),x1.*(t103.*th12.*(t54-t71).*(1.69e+2./2.0e+2)+t54.*t104.*t122.*th12.*(1.69e+2./2.0)),0.0,-x1.*(t224+t244+t246+t103.*t146.*(3.1e+1./2.0e+2)+t207.*th9+t12.*t35.*t41.*t80.*t104.*t122.*th3.*th12.*x3.*x4.*1.834319526627219e+1),0.0];
mt4 = [x1.*(t103.*t146-th14.*(t223+t245-t161.*t164.*t207-t161.*t164.*t208.*(t46-t150)+t162.*t165.*t208.*(t46-t150)+t161.*t164.*t208.*th10.*(t46-t150))+t12.*t35.*t39.*t41.*t80.*t104.*t122.*th12.*x3.*x4),x1.*(t103.*t146.*th12.*(1.69e+2./2.0e+2)-t12.*t35.*t41.*th3.*x3.*x4+t7.*t12.*t35.*t41.*t80.*t104.*t122.*th3.*x3.*x4.*1.0e+2),0.0,x1.*(t215-t232+t233+th9.*(t170+t107.*(t90-t98))),-t63.*x1];
mt5 = [-th14.*x1.*(t215-t232+t233+t161.*t164.*(-t63+t170+t107.*(t90-t98))+t161.*t164.*(t46-t150).*(-t66+t172+t107.*(t90-t98).*1.0e+2)-t162.*t165.*(t46-t150).*(-t66+t172+t107.*(t90-t98).*1.0e+2)),0.0,0.0,x1.*(t150+t231-t247+t248-th4+th9.*(t189-t193)+t15.*t18.*t32.*3.025e-4),t15.*t18.*t32.*(x2-4.5e+2).*3.025e-4];
mt6 = [-th14.*x1.*(t231-t247+t248+t161.*t164.*(t189-t193)+t161.*t164.*(t46-t150).*(t191-t194)-t162.*t165.*(t46-t150).*(t191-t194))+t15.*t18.*t32.*x3.*3.025e-4,t15.*t18.*t32.*x4.*6.05e-4,t15.*t32.*(-3.025e-4),t161.*t164.*x1.*(t46-t150),0.0,-t161.*t164.*th14.*x1.*(t46-t150),0.0,0.0,-x1.*(t221+t242-t204.*th9+t161.*t164.*t206.*th10.*(t46-t150)),0.0,th14.*x1.*(t221+t242-t161.*t164.*t204-t161.*t164.*t206.*(t46-t150)+t162.*t165.*t206.*(t46-t150)+t161.*t164.*t206.*th10.*(t46-t150)),0.0,0.0];
mt7 = [-x1.*(t103.*(t133+t80.*(t65-th3).*(2.0e+2./1.69e+2)).*(3.1e+1./2.0e+2)-t80.*t86.*t104.*t122.*(3.1e+1./2.0e+2)),0.0,x1.*(t103.*(t133+t80.*(t65-th3).*(2.0e+2./1.69e+2))-t80.*t86.*t104.*t122),-x1.*(t103.*t122.*(1.69e+2./2.0e+2)-t103.*th12.*(t133+t80.*(t65-th3).*(2.0e+2./1.69e+2)).*(1.69e+2./2.0e+2)+t80.*t86.*t104.*t122.*th12.*(1.69e+2./2.0e+2)),0.0,0.0,0.0,x1.*(t161.*t164.*(t46-t150)-t161.*t164.*th10.*(t46-t150)),0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7],5,13);
