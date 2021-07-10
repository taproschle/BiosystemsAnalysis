function dfdxSym = dfdxDewasme(t,in2,in3)
%DFDXDEWASME
%    DFDXSYM = DFDXDEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    06-Jul-2021 16:38:40

th2 = in3(1,:);
th3 = in3(2,:);
th5 = in3(3,:);
th7 = in3(4,:);
th10 = in3(5,:);
th11 = in3(6,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th5+x3;
t3 = th7+x3;
t4 = th5.^2;
t5 = th11.^2;
t6 = x3.^2;
t7 = x4.^2;
t8 = 1.0./th3;
t9 = 1.0./x5;
t16 = th5.*1.719e+3;
t17 = x3.*1.719e+3;
t18 = x2.*5.0e+3;
t19 = x4.*1.0e+4;
t20 = t.*(1.1e+1./1.0e+2);
t23 = th2.*x2.*5.0e+5;
t24 = x4+1.0e-4;
t27 = th5.*x4.*1.719e+7;
t28 = x3.*x4.*1.719e+7;
t29 = x2+1.414e-1;
t39 = th5.*th11.*x4.*5.0e+9;
t10 = t9.^2;
t11 = 1.0./t2;
t14 = 1.0./t3;
t21 = t19+1.0;
t22 = exp(t20);
t25 = -t23;
t26 = t18+7.07e+2;
t30 = 1.0./t24;
t35 = 1.0./t29;
t40 = -t39;
t12 = t11.^2;
t13 = t11.^3;
t15 = t14.^2;
t31 = t30.^2;
t32 = 1.0./t21;
t36 = 1.0./t26;
t41 = t35.*th2.*x2;
t45 = t25+t26;
t54 = t8.*t9.*t22.*3.025e-4;
t56 = t11.*t30.*th5.*th11.*x4.*(5.0e+2./9.99e+2);
t58 = t11.*t30.*th5.*th11.*x4.*2.908667830133799;
t60 = t11.*t30.*th5.*th11.*x4.*2.908667830133799e+2;
t63 = t11.*t30.*th5.*th11.*x4.*5.817335660267597e+2;
t78 = t2.*1.719e+3+t27+t28+t40;
t33 = t32.^2;
t34 = t32.^3;
t37 = t36.^2;
t38 = t36.^3;
t42 = t41.*1.0e+2;
t43 = t41.*2.0e+2;
t46 = t18.*t36.*th2;
t48 = t23.*t36;
t49 = t36.*th2.*x2.*-5.0e+5;
t50 = t36.*th2.*x2.*1.0e+6;
t52 = t41.*1.720720720720721e-1;
t53 = t36.*th2.*x2.*8.603603603603604e+2;
t55 = -t54;
t57 = -t56;
t59 = -t58;
t61 = t11.*t32.*th5.*th11.*x4.*5.005005005005005e+3;
t62 = -t60;
t65 = t11.*t32.*th5.*th11.*x4.*2.908667830133799e+4;
t67 = -t63;
t69 = t11.*t32.*th5.*th11.*x4.*2.908667830133799e+6;
t71 = t11.*t32.*th5.*th11.*x4.*5.817335660267597e+6;
t75 = t12.*t14.*t30.*th5.*th11.*x3.*x4.*5.005005005005005e+1;
t44 = -t42;
t51 = exp(t49);
t64 = -t61;
t66 = exp(t62);
t68 = -t65;
t70 = -t69;
t72 = -t71;
t77 = t41+t59;
t79 = t42+t62;
t80 = t43+t67;
t82 = t52+t57;
t47 = exp(t44);
t74 = exp(t70);
t76 = t46.*t51;
t81 = exp(t79);
t83 = exp(t80);
t84 = t46+t68;
t91 = t48+t70;
t92 = t53+t64;
t94 = t50+t72;
t98 = t58.*t66;
t99 = t12.*t30.*t66.*th5.*th11.*x4.*2.908667830133799;
t100 = t14.*t82.*1.0e+2;
t102 = t15.*t82.*x3.*1.0e+2;
t103 = t14.*t82.*x3.*2.0e+2;
t108 = t14.*t82.*x3.*-1.0e+2;
t124 = t4.*t5.*t7.*t13.*t31.*t66.*8.460348546055261e+2;
t73 = t41.*t47;
t85 = t81+1.0;
t86 = t47+t66;
t93 = exp(t91);
t96 = t51+t74;
t97 = exp(t94);
t101 = t100.*x3;
t109 = -t102;
t110 = -t103;
t111 = exp(t108);
t114 = t65.*t74;
t115 = t14.*t92.*x3.*1.0e+2;
t116 = t14.*t92.*x3.*2.0e+2;
t127 = -t124;
t87 = 1.0./t85;
t89 = 1.0./t86;
t95 = t93+1.0;
t106 = 1.0./t96;
t112 = exp(t110);
t113 = t111+1.0;
t119 = -t115;
t120 = -t116;
t128 = t73+t98;
t129 = t76+t114;
t130 = t99+t127;
t131 = t75+t100+t109;
t88 = t87.^2;
t90 = t89.^2;
t104 = 1.0./t95;
t107 = t106.^2;
t117 = 1.0./t113;
t121 = exp(t119);
t122 = exp(t120);
t105 = t104.^2;
t118 = t117.^2;
t123 = t121+1.0;
t125 = 1.0./t123;
t126 = t125.^2;
mt1 = [t55+t77.*t81.*t87.*3.294e-1+t89.*t128.*th3-t14.*t82.*t111.*t117.*th10.*x3,-t89.*t128+t77.*t81.*t87,t77.*t81.*t87.*9.674e-1+t14.*t82.*t111.*t117.*x3,0.0,0.0];
mt2 = [x1.*(t37.*t93.*t104.*th2.*1.164429e+6+t37.*t84.*t93.*t104.*th2.*1.164429e+8-t37.*t84.*t97.*t105.*th2.*1.164429e+8+t38.*t45.*t51.*t106.*th2.*th3.*3.535e+6+t37.*t51.*t107.*t129.*th2.*th3.*3.535e+8-t14.*t37.*t121.*t125.*th2.*th10.*x3.*6.082747747747748e+5+t6.*t15.*t37.*t92.*t121.*t125.*th2.*th10.*6.082747747747748e+7-t6.*t15.*t37.*t92.*t122.*t126.*th2.*th10.*6.082747747747748e+7)];
mt3 = [t55-x1.*(t37.*t93.*t104.*th2.*-3.535e+6+t38.*t45.*t51.*t106.*th2.*3.535e+6-t37.*t84.*t93.*t104.*th2.*3.535e+8+t37.*t84.*t97.*t105.*th2.*3.535e+8+t37.*t51.*t107.*t129.*th2.*3.535e+8),x1.*(t37.*t93.*t104.*th2.*3.419759e+6+t37.*t84.*t93.*t104.*th2.*3.419759e+8-t37.*t84.*t97.*t105.*th2.*3.419759e+8+t14.*t37.*t121.*t125.*th2.*x3.*6.082747747747748e+5-t6.*t15.*t37.*t92.*t121.*t125.*th2.*6.082747747747748e+7+t6.*t15.*t37.*t92.*t122.*t126.*th2.*6.082747747747748e+7)];
mt4 = [-x4.*(t37.*t93.*t104.*th2.*(7.07e+2./2.0)+t38.*t45.*t51.*t106.*th2.*1.215333e+6+t37.*t84.*t93.*t104.*th2.*3.535e+4-t37.*t84.*t97.*t105.*th2.*3.535e+4+t37.*t51.*t107.*t129.*th2.*1.215333e+8-t14.*t37.*t121.*t125.*th2.*x3.*1.215333e+6+t6.*t15.*t37.*t92.*t121.*t125.*th2.*1.215333e+8-t6.*t15.*t37.*t92.*t122.*t126.*th2.*1.215333e+8),0.0];
mt5 = [-x1.*(t89.*t130.*th3+t14.*t82.*t111.*t117.*th10-t15.*t82.*t111.*t117.*th10.*x3-t14.*t82.*t111.*t117.*t131.*th10.*x3+t14.*t82.*t112.*t118.*t131.*th10.*x3-t12.*t30.*t81.*t87.*th5.*th11.*x4.*(1.83e+2./1.91e+2)-t12.*t30.*t77.*t81.*t87.*th5.*th11.*x4.*9.581151832460733e+1+t12.*t30.*t77.*t83.*t88.*th5.*th11.*x4.*9.581151832460733e+1+t12.*t30.*t66.*t90.*t128.*th3.*th5.*th11.*x4.*2.908667830133799e+2+t12.*t14.*t30.*t111.*t117.*th5.*th10.*th11.*x3.*x4.*(5.0e+2./9.99e+2))];
mt6 = [x1.*(t12.*t32.*t93.*t104.*th5.*th11.*x4.*2.908667830133799e+4+t13.*t33.*t74.*t78.*t106.*th5.*th11.*x4.*1.692069709211052e+1+t12.*t32.*t84.*t93.*t104.*th5.*th11.*x4.*2.908667830133799e+6-t12.*t32.*t84.*t97.*t105.*th5.*th11.*x4.*2.908667830133799e+6+t12.*t32.*t74.*t107.*t129.*th5.*th11.*x4.*2.908667830133799e+6)];
mt7 = [t55+x1.*(t14.*t82.*t111.*t117-t15.*t82.*t111.*t117.*x3-t14.*t82.*t111.*t117.*t131.*x3+t14.*t82.*t112.*t118.*t131.*x3+t12.*t30.*t81.*t87.*th5.*th11.*x4.*2.813845258871437+t12.*t30.*t77.*t81.*t87.*th5.*th11.*x4.*2.813845258871437e+2-t12.*t30.*t77.*t83.*t88.*th5.*th11.*x4.*2.813845258871437e+2+t12.*t14.*t30.*t111.*t117.*th5.*th11.*x3.*x4.*(5.0e+2./9.99e+2))];
mt8 = [x4.*(t89.*t130.*3.438e-1+t14.*t82.*t111.*t117.*(9.99e+2./5.0e+2)-t15.*t82.*t111.*t117.*x3.*(9.99e+2./5.0e+2)-t14.*t82.*t111.*t117.*t131.*x3.*(9.99e+2./5.0e+2)+t14.*t82.*t112.*t118.*t131.*x3.*(9.99e+2./5.0e+2)-(t12.*t30.*t81.*t87.*th5.*th11.*x4)./3.438e+3-t12.*t30.*t77.*t81.*t87.*th5.*th11.*x4.*2.908667830133799e-2+t12.*t30.*t77.*t83.*t88.*th5.*th11.*x4.*2.908667830133799e-2+t12.*t30.*t66.*t90.*t128.*th5.*th11.*x4.*1.0e+2+t12.*t14.*t30.*t111.*t117.*th5.*th11.*x3.*x4),0.0];
mt9 = [x1.*(t11.*t33.*t93.*t104.*th5.*th11.*(-9.581151832460733e+3)-t11.*t33.*t84.*t93.*t104.*th5.*th11.*9.581151832460733e+5+t11.*t33.*t84.*t97.*t105.*th5.*th11.*9.581151832460733e+5+t12.*t34.*t74.*t78.*t106.*th3.*th5.*th11.*1.692069709211052e+1+t11.*t33.*t74.*t107.*t129.*th3.*th5.*th11.*2.908667830133799e+6+t11.*t14.*t33.*t121.*t125.*th5.*th10.*th11.*x3.*5.005005005005005e+3-t6.*t11.*t15.*t33.*t92.*t121.*t125.*th5.*th10.*th11.*5.005005005005005e+5+t6.*t11.*t15.*t33.*t92.*t122.*t126.*th5.*th10.*th11.*5.005005005005005e+5)];
mt10 = [-x1.*(t11.*t33.*t93.*t104.*th5.*th11.*2.908667830133799e+4+t12.*t34.*t74.*t78.*t106.*th5.*th11.*1.692069709211052e+1+t11.*t33.*t84.*t93.*t104.*th5.*th11.*2.908667830133799e+6-t11.*t33.*t84.*t97.*t105.*th5.*th11.*2.908667830133799e+6+t11.*t33.*t74.*t107.*t129.*th5.*th11.*2.908667830133799e+6)];
mt11 = [-x1.*(t11.*t33.*t93.*t104.*th5.*th11.*2.813845258871437e+4+t11.*t33.*t84.*t93.*t104.*th5.*th11.*2.813845258871437e+6-t11.*t33.*t84.*t97.*t105.*th5.*th11.*2.813845258871437e+6+t11.*t14.*t33.*t121.*t125.*th5.*th11.*x3.*5.005005005005005e+3-t6.*t11.*t15.*t33.*t92.*t121.*t125.*th5.*th11.*5.005005005005005e+5+t6.*t11.*t15.*t33.*t92.*t122.*t126.*th5.*th11.*5.005005005005005e+5)];
mt12 = [t55-t106.*t129.*3.438e-1-x4.*(t11.*t33.*t93.*t104.*th5.*th11.*(-2.908667830133799)+t12.*t34.*t74.*t78.*t106.*th5.*th11.*5.817335660267597-t11.*t33.*t84.*t93.*t104.*th5.*th11.*2.908667830133799e+2+t11.*t33.*t84.*t97.*t105.*th5.*th11.*2.908667830133799e+2+t11.*t33.*t74.*t107.*t129.*th5.*th11.*1.0e+6+t11.*t14.*t33.*t121.*t125.*th5.*th11.*x3.*1.0e+4-t6.*t11.*t15.*t33.*t92.*t121.*t125.*th5.*th11.*1.0e+6+t6.*t11.*t15.*t33.*t92.*t122.*t126.*th5.*th11.*1.0e+6)-(t84.*t93.*t104)./1.0e+4+t14.*t92.*t121.*t125.*x3.*(9.99e+2./5.0e+2)-1.8e+4,0.0];
mt13 = [t8.*t10.*t22.*x1.*3.025e-4,t8.*t10.*t22.*(x2-4.5e+2).*3.025e-4,t8.*t10.*t22.*x3.*3.025e-4,t8.*t10.*t22.*x4.*3.025e-4,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12,mt13],5,5);