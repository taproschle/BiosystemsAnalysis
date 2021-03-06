function dfdthSym = dfdthXu(t,in2,in3)
%DFDTHXU
%    DFDTHSYM = DFDTHXU(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 16:07:23

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
th13 = in3(13,:);
th14 = in3(14,:);
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
t11 = th13-1.0;
t12 = 1.0./th7;
t14 = 1.0./th8;
t16 = 1.0./th9.^2;
t17 = 1.0./th11;
t19 = 1.0./x5;
t35 = t.*(1.1e+1./1.0e+2);
t37 = x4+1.0e-4;
t13 = t12.^2;
t15 = t14.^2;
t18 = t17.^2;
t20 = t12.*x3;
t21 = t14.*x3;
t22 = 1.0./t3;
t25 = 1.0./t4;
t28 = 1.0./t10;
t30 = 1.0./t11;
t36 = exp(t35);
t44 = 1.0./t37;
t23 = t22.^2;
t24 = t22.^3;
t26 = t25.^2;
t27 = t25.^3;
t29 = t28.^2;
t31 = t30.^2;
t32 = t30.^3;
t33 = t20+1.0;
t34 = t21+1.0;
t43 = t28.*th4.*1.0e+2;
t45 = t25.*th2.*x3.*1.0e+2;
t38 = 1.0./t33;
t40 = 1.0./t34;
t46 = -t45;
t39 = t38.^2;
t41 = t40.^2;
t42 = t40.^3;
t47 = exp(t46);
t48 = t22.*t40.*x2;
t50 = t23.*t40.*th1.*x2;
t65 = t38.*t44.*x4;
t49 = t48.*th1;
t51 = -t48;
t53 = -t50;
t54 = t25.*t47.*x3;
t55 = t48.*1.0e+2;
t57 = t26.*t47.*th2.*x3;
t59 = t50.*1.0e+2;
t66 = t15.*t22.*t41.*th1.*x2.*x3;
t68 = t65.*th3;
t70 = t9.*t26.*t47.*th2.*1.0e+2;
t71 = t65-1.0;
t74 = t6.*t9.*t27.*t47.*1.0e+2;
t52 = -t49;
t56 = t54.*th2;
t58 = t49.*1.0e+2;
t60 = t49.*2.0e+2;
t61 = -t55;
t63 = -t59;
t64 = -t57;
t69 = t66.*1.0e+2;
t72 = -t68;
t73 = -t70;
t76 = t17.*t68;
t77 = (t68-th3).^2;
t80 = t30.*th12.*(t68-th3).*-1.0e+2;
t62 = -t58;
t75 = t72+th3;
t78 = t17.*t72;
t81 = exp(t80);
t95 = t54+t73;
t97 = t64+t74;
t67 = exp(t62);
t79 = t2+t78;
t101 = t47+t81;
t109 = t30.*t71.*t81.*th12;
t110 = -t30.*t81.*(t68-th3);
t111 = t13.*t30.*t39.*t44.*t81.*th3.*th12.*x3.*x4;
t113 = -t31.*t81.*th12.*(t68-th3);
t114 = t30.*t81.*th12.*(t68-th3);
t117 = t31.*t77.*t81.*th12.*1.0e+2;
t118 = t7.*t32.*t77.*t81.*1.0e+2;
t130 = t7.*t31.*t71.*t81.*(t68-th3).*-1.0e+2;
t131 = t7.*t13.*t31.*t39.*t44.*t81.*th3.*x3.*x4.*(t68-th3).*-1.0e+2;
t82 = t48.*t67;
t83 = t49.*t67;
t84 = t50.*t67;
t85 = t28.*t79.*1.0e+2;
t86 = t29.*t79.*1.0e+2;
t90 = t66.*t67;
t91 = t8.*t23.*t41.*t67.*th1.*1.0e+2;
t94 = t5.*t8.*t24.*t41.*t67.*1.0e+2;
t98 = t5.*t8.*t15.*t23.*t42.*t67.*x3.*1.0e+2;
t103 = 1.0./t101;
t112 = t110.*th12;
t126 = t56+t114;
t143 = t110+t117;
t145 = t113+t118;
t156 = t109+t130;
t158 = t111+t131;
t87 = -t85;
t88 = -t86;
t92 = -t91;
t93 = -t90;
t96 = -t94;
t104 = t103.^2;
t128 = t95.*t103;
t132 = -t103.*(t57-t74);
t144 = t103.*t126;
t167 = t103.*(t117-t30.*t81.*(t68-th3));
t169 = t103.*t145;
t189 = t103.*t156;
t191 = t103.*t158;
t89 = exp(t87);
t105 = t43+t88;
t122 = t82+t92;
t125 = t84+t96;
t127 = t93+t98;
t129 = t128.*th13;
t133 = t132.*th13;
t152 = t54.*t104.*t126.*1.0e+2;
t153 = t57.*t104.*t126.*1.0e+2;
t168 = t167.*th13;
t170 = t169.*th13;
t181 = t104.*t109.*t126.*1.0e+2;
t182 = t30.*t81.*t104.*t126.*(t68-th3).*-1.0e+2;
t184 = t31.*t81.*t104.*t126.*th12.*(t68-th3).*-1.0e+2;
t186 = t104.*t111.*t126.*1.0e+2;
t190 = t189.*th13;
t192 = t191.*th13;
t99 = t28.*t89.*th4;
t100 = t28.*t89.*th9;
t106 = t67+t89;
t115 = t17.*t28.*t65.*t89;
t116 = t18.*t28.*t68.*t89;
t119 = t13.*t17.*t28.*t39.*t44.*t89.*th3.*x3.*x4;
t120 = t28.*t79.*t89;
t121 = t29.*t79.*t89;
t123 = t86.*t89.*th9;
t134 = t17.*t65.*t86.*t89;
t136 = t18.*t68.*t86.*t89;
t138 = t13.*t17.*t39.*t44.*t86.*t89.*th3.*x3.*x4;
t154 = t152.*th13;
t155 = -t153;
t157 = t153.*th13;
t183 = t181.*th13;
t185 = t182.*th13;
t187 = t184.*th13;
t188 = t186.*th13;
t102 = -t99;
t107 = 1.0./t106;
t124 = t121.*th9.*-1.0e+2;
t135 = t17.*t65.*t121.*-1.0e+2;
t137 = t18.*t68.*t121.*-1.0e+2;
t139 = t13.*t17.*t39.*t44.*t121.*th3.*x3.*x4.*-1.0e+2;
t140 = t83+t120;
t141 = t105.*t120;
t108 = t107.^2;
t142 = t100+t124;
t146 = t107.*t122;
t147 = t107.*t125;
t150 = -t107.*(t90-t98);
t151 = t107.*(t90-t98).*-1.0e+2;
t159 = t115+t135;
t160 = t116+t137;
t161 = t107.*t140;
t166 = t119+t139;
t213 = t102+t121+t141;
t148 = t146.*1.0e+2;
t149 = t147.*1.0e+2;
t162 = t161.*1.0e+2;
t163 = t161.*2.0e+2;
t171 = t107.*t142;
t173 = t52+t161;
t193 = t55.*t67.*t108.*t140;
t194 = t59.*t67.*t108.*t140;
t195 = t82.*t108.*t140.*1.0e+4;
t196 = t84.*t108.*t140.*1.0e+4;
t197 = t67.*t69.*t108.*t140;
t198 = t90.*t108.*t140.*-1.0e+2;
t199 = t90.*t108.*t140.*1.0e+4;
t201 = t100.*t108.*t140.*1.0e+2;
t202 = t100.*t108.*t140.*1.0e+4;
t203 = t107.*t159;
t204 = t107.*t160;
t207 = t107.*t166;
t209 = t108.*t115.*t140.*1.0e+2;
t210 = t108.*t116.*t140.*1.0e+2;
t211 = t108.*t115.*t140.*1.0e+4;
t212 = t108.*t116.*t140.*1.0e+4;
t214 = t108.*t119.*t140.*1.0e+2;
t215 = t108.*t119.*t140.*1.0e+4;
t216 = t89.*t105.*t108.*t140;
t220 = t107.*t213;
t164 = -t162;
t165 = -t163;
t172 = t171.*1.0e+2;
t200 = -t199;
t205 = t203.*1.0e+2;
t206 = t204.*1.0e+2;
t208 = t207.*1.0e+2;
t217 = -t216;
t218 = t216.*1.0e+2;
t221 = t220.*1.0e+2;
t222 = t51+t146+t193;
t223 = t53+t147+t194;
t224 = t61+t148+t195;
t225 = t63+t149+t196;
t226 = t66+t150+t198;
t228 = t171+t201;
t230 = t203+t209;
t231 = t204+t210;
t234 = t207+t214;
t174 = t58+t164;
t175 = t60+t165;
t219 = -t218;
t227 = t69+t151+t200;
t229 = t172+t202;
t232 = t205+t211;
t233 = t206+t212;
t235 = t208+t215;
t236 = t217+t220;
t176 = exp(t174);
t177 = exp(t175);
t237 = t219+t221;
t178 = t176+1.0;
t179 = 1.0./t178;
t180 = t179.^2;
t238 = t176.*t179.*t222.*th10;
t240 = t176.*t179.*t223.*th10;
t242 = -t176.*t179.*th10.*(-t66+t197+t107.*(t90-t98));
t243 = t176.*t179.*th10.*(-t66+t197+t107.*(t90-t98));
t244 = t176.*t179.*t228.*th10;
t246 = t176.*t179.*t230.*th10;
t247 = t176.*t179.*t231.*th10;
t249 = t176.*t179.*t234.*th10;
t250 = -t176.*t179.*t224.*th10.*(t49-t161);
t253 = -t176.*t179.*t225.*th10.*(t49-t161);
t256 = -t176.*t179.*th10.*(t216-t220);
t257 = t176.*t179.*th10.*(t49-t161).*(-t69+t199+t107.*(t90-t98).*1.0e+2);
t260 = -t176.*t179.*t229.*th10.*(t49-t161);
t263 = -t176.*t179.*t232.*th10.*(t49-t161);
t265 = -t176.*t179.*t233.*th10.*(t49-t161);
t268 = -t176.*t179.*t235.*th10.*(t49-t161);
t270 = t176.*t179.*th10.*(t49-t161).*(t218-t221);
t239 = -t238;
t241 = -t240;
t245 = -t244;
t248 = -t247;
t251 = -t177.*t180.*t224.*th10.*(t49-t161);
t252 = t177.*t180.*t224.*th10.*(t49-t161);
t254 = -t177.*t180.*t225.*th10.*(t49-t161);
t255 = t177.*t180.*t225.*th10.*(t49-t161);
t258 = t177.*t180.*th10.*(t49-t161).*(-t69+t199+t107.*(t90-t98).*1.0e+2);
t261 = -t177.*t180.*t229.*th10.*(t49-t161);
t262 = t177.*t180.*t229.*th10.*(t49-t161);
t264 = -t177.*t180.*t232.*th10.*(t49-t161);
t266 = -t177.*t180.*t233.*th10.*(t49-t161);
t267 = t177.*t180.*t233.*th10.*(t49-t161);
t269 = -t177.*t180.*t235.*th10.*(t49-t161);
t271 = t177.*t180.*th10.*(t49-t161).*(t218-t221);
t259 = -t258;
mt1 = [-x1.*(t238+t251-th9.*(t146+t193)+t176.*t179.*t224.*th10.*(t49-t161)),t51.*x1,th14.*x1.*(t238+t251-t176.*t179.*t222-t176.*t179.*t224.*(t49-t161)+t177.*t180.*t224.*(t49-t161)+t176.*t179.*t224.*th10.*(t49-t161)),0.0,0.0,x1.*(t129+t154),0.0,-x1.*(t128+t152),-th12.*x1.*(t128-t129+t152-t54.*t104.*t126.*th13.*1.0e+2),0.0,x1.*(t183+t190+t246+t264-t230.*th9+t176.*t179.*t232.*th10.*(t49-t161)),0.0,-x1.*(t181+t189+th14.*(t246+t264-t176.*t179.*t230-t176.*t179.*t232.*(t49-t161)+t177.*t180.*t232.*(t49-t161)+t176.*t179.*t232.*th10.*(t49-t161)))];
mt2 = [-x1.*(t65+th12.*(t181+t189-t190-t104.*t109.*t126.*th13.*1.0e+2)),0.0,-x1.*(t244+t261-th9.*(t228-1.0)+t176.*t179.*t229.*th10.*(t49-t161)),0.0,th14.*x1.*(t244+t261-t176.*t179.*t228-t176.*t179.*t229.*(t49-t161)+t177.*t180.*t229.*(t49-t161)+t176.*t179.*t229.*th10.*(t49-t161)),0.0,0.0,x1.*(t240+t254-th9.*(t147+t194)+t176.*t179.*t225.*th10.*(t49-t161)),t50.*x1,-th14.*x1.*(t240+t254-t176.*t179.*t223-t176.*t179.*t225.*(t49-t161)+t177.*t180.*t225.*(t49-t161)+t176.*t179.*t225.*th10.*(t49-t161)),0.0,0.0,-x1.*(t157+t103.*th13.*(t57-t74)),0.0];
mt3 = [x1.*(t153+t103.*(t57-t74)),th12.*x1.*(t133+t153+t103.*(t57-t74)-t57.*t104.*t126.*th13.*1.0e+2),0.0,x1.*(t188+t192+t249+t269-t234.*th9+t176.*t179.*t235.*th10.*(t49-t161)),0.0,-x1.*(t186+t191+th14.*(t249+t269-t176.*t179.*t234-t176.*t179.*t235.*(t49-t161)+t177.*t180.*t235.*(t49-t161)+t176.*t179.*t235.*th10.*(t49-t161))),-x1.*(th12.*(t186+t191-t192-t104.*t111.*t126.*th13.*1.0e+2)+t13.*t39.*t44.*th3.*x3.*x4),0.0,x1.*(t242-t257+t258+th9.*(t197+t107.*(t90-t98))),-t66.*x1];
mt4 = [-th14.*x1.*(t242-t257+t258+t176.*t179.*(-t66+t197+t107.*(t90-t98))+t176.*t179.*(t49-t161).*(-t69+t199+t107.*(t90-t98).*1.0e+2)-t177.*t180.*(t49-t161).*(-t69+t199+t107.*(t90-t98).*1.0e+2)),0.0,0.0,x1.*(t161+t256-t270+t271-th4+th9.*(t216-t220)+t16.*t19.*t36.*3.025e-4),t16.*t19.*t36.*(x2-4.5e+2).*3.025e-4];
mt5 = [-th14.*x1.*(t256-t270+t271+t176.*t179.*(t216-t220)+t176.*t179.*(t49-t161).*(t218-t221)-t177.*t180.*(t49-t161).*(t218-t221))+t16.*t19.*t36.*x3.*3.025e-4,t16.*t19.*t36.*x4.*6.05e-4,t16.*t36.*(-3.025e-4),t176.*t179.*x1.*(t49-t161),0.0,-t176.*t179.*th14.*x1.*(t49-t161),0.0,0.0,-x1.*(t247+t266-t231.*th9+t176.*t179.*t233.*th10.*(t49-t161)),0.0,th14.*x1.*(t247+t266-t176.*t179.*t231-t176.*t179.*t233.*(t49-t161)+t177.*t180.*t233.*(t49-t161)+t176.*t179.*t233.*th10.*(t49-t161)),0.0,0.0];
mt6 = [-x1.*(t168+t185),0.0,x1.*(t167+t182),x1.*(-t144+th12.*(t167-t168+t182+t30.*t81.*t104.*t126.*th13.*(t68-th3).*1.0e+2)+t144.*th13),0.0,x1.*(t144+t170+t187),0.0,-x1.*(t169+t184),th12.*x1.*(t144-t169+t170+t187+t31.*t81.*t104.*t126.*th12.*(t68-th3).*1.0e+2),0.0,0.0,0.0,x1.*(t176.*t179.*(t49-t161)-t176.*t179.*th10.*(t49-t161)),0.0,0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6],5,14);
