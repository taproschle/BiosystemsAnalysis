function dfdthSym = dfdthDewasme(t,in2,in3)
%DFDTHDEWASME
%    DFDTHSYM = DFDTHDEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    06-Jul-2021 16:38:42

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
t5 = x2.^2;
t6 = x3.^2;
t7 = x3.^3;
t8 = x4.^2;
t9 = 1.0./th3.^2;
t10 = 1.0./x5;
t17 = th5.*1.719e+3;
t18 = x3.*1.719e+3;
t19 = x2.*5.0e+3;
t20 = x4.*1.0e+4;
t21 = t.*(1.1e+1./1.0e+2);
t24 = x4+1.0e-4;
t26 = th5.*x4.*1.719e+7;
t27 = x3.*x4.*1.719e+7;
t28 = x2+1.414e-1;
t36 = th5.*th11.*x4.*5.0e+9;
t11 = 1.0./t2;
t14 = 1.0./t3;
t22 = t20+1.0;
t23 = exp(t21);
t25 = t19+7.07e+2;
t29 = 1.0./t24;
t33 = 1.0./t28;
t37 = -t36;
t12 = t11.^2;
t13 = t11.^3;
t15 = t14.^2;
t16 = t14.^3;
t30 = t29.^2;
t31 = 1.0./t22;
t34 = t33.^2;
t35 = 1.0./t25;
t38 = t33.*th2.*x2;
t50 = t11.*t29.*th5.*th11.*x4.*(5.0e+2./9.99e+2);
t52 = t11.*t29.*th5.*th11.*x4.*2.908667830133799;
t54 = t11.*t29.*th5.*th11.*x4.*2.908667830133799e+2;
t57 = t11.*t29.*th5.*th11.*x4.*5.817335660267597e+2;
t74 = t2.*1.719e+3+t26+t27+t37;
t32 = t31.^2;
t39 = t38.*1.0e+2;
t40 = t38.*2.0e+2;
t42 = t19.*t35.*th2;
t44 = t35.*th2.*x2.*5.0e+5;
t46 = t35.*th2.*x2.*1.0e+6;
t48 = t38.*1.720720720720721e-1;
t49 = t35.*th2.*x2.*8.603603603603604e+2;
t51 = -t50;
t53 = -t52;
t55 = t11.*t31.*th5.*th11.*x4.*5.005005005005005e+3;
t56 = -t54;
t59 = t11.*t31.*th5.*th11.*x4.*2.908667830133799e+4;
t61 = -t57;
t63 = t11.*t31.*th5.*th11.*x4.*2.908667830133799e+6;
t65 = t11.*t31.*th5.*th11.*x4.*5.817335660267597e+6;
t41 = -t39;
t45 = -t44;
t58 = -t55;
t60 = exp(t56);
t62 = -t59;
t64 = -t63;
t67 = -t65;
t73 = t38+t53;
t75 = t39+t56;
t76 = t40+t61;
t78 = t48+t51;
t43 = exp(t41);
t47 = exp(t45);
t69 = exp(t64);
t77 = exp(t75);
t79 = exp(t76);
t80 = t42+t62;
t83 = t78.^2;
t88 = t44+t64;
t89 = t49+t58;
t91 = t46+t67;
t93 = t11.*t29.*t60.*th5.*x4.*2.908667830133799;
t96 = t52.*t60;
t97 = t14.*t78.*x3.*1.0e+2;
t98 = t14.*t78.*x3.*2.0e+2;
t118 = t4.*t8.*t12.*t30.*t60.*th11.*8.460348546055261e+2;
t66 = t33.*t43.*x2;
t68 = t38.*t43;
t70 = t5.*t34.*t43.*th2.*1.0e+2;
t72 = t42.*t47;
t81 = t77+1.0;
t82 = t43+t60;
t90 = exp(t88);
t94 = t47+t69;
t95 = exp(t91);
t103 = -t97;
t104 = -t98;
t108 = t59.*t69;
t109 = t14.*t89.*x3.*1.0e+2;
t110 = t14.*t89.*x3.*2.0e+2;
t119 = -t118;
t71 = -t70;
t84 = 1.0./t81;
t86 = 1.0./t82;
t92 = t90+1.0;
t101 = 1.0./t94;
t105 = exp(t103);
t106 = exp(t104);
t113 = -t109;
t114 = -t110;
t123 = t68+t96;
t124 = t72+t108;
t125 = t93+t119;
t85 = t84.^2;
t87 = t86.^2;
t99 = 1.0./t92;
t102 = t101.^2;
t107 = t105+1.0;
t115 = exp(t113);
t116 = exp(t114);
t122 = t66+t71;
t100 = t99.^2;
t111 = 1.0./t107;
t117 = t115+1.0;
t112 = t111.^2;
t120 = 1.0./t117;
t121 = t120.^2;
mt1 = [x1.*(t86.*t122.*th3+t66.*t87.*t123.*th3.*1.0e+2+t33.*t77.*t84.*x2.*3.294e-1+t33.*t73.*t77.*t84.*x2.*3.294e+1-t33.*t73.*t79.*t85.*x2.*3.294e+1-t14.*t33.*t105.*t111.*th10.*x2.*x3.*1.720720720720721e-1+t6.*t15.*t33.*t78.*t105.*t111.*th10.*x2.*1.720720720720721e+1-t6.*t15.*t33.*t78.*t106.*t112.*th10.*x2.*1.720720720720721e+1),-x1.*(t86.*t122+t66.*t87.*t123.*1.0e+2-t33.*t77.*t84.*x2-t33.*t73.*t77.*t84.*x2.*1.0e+2+t33.*t73.*t79.*t85.*x2.*1.0e+2)];
mt2 = [x1.*(t33.*t77.*t84.*x2.*9.674e-1+t33.*t73.*t77.*t84.*x2.*9.674e+1-t33.*t73.*t79.*t85.*x2.*9.674e+1+t14.*t33.*t105.*t111.*x2.*x3.*1.720720720720721e-1-t6.*t15.*t33.*t78.*t105.*t111.*x2.*1.720720720720721e+1+t6.*t15.*t33.*t78.*t106.*t112.*x2.*1.720720720720721e+1)];
mt3 = [-x4.*(t86.*t122.*3.438e-1+t66.*t87.*t123.*3.438e+1+(t33.*t77.*t84.*x2)./1.0e+4+(t33.*t73.*t77.*t84.*x2)./1.0e+2-(t33.*t73.*t79.*t85.*x2)./1.0e+2-t14.*t33.*t105.*t111.*x2.*x3.*3.438e-1+t6.*t15.*t33.*t78.*t105.*t111.*x2.*3.438e+1-t6.*t15.*t33.*t78.*t106.*t112.*x2.*3.438e+1),0.0,x1.*(t86.*t123+t9.*t10.*t23.*3.025e-4),t9.*t10.*t23.*(x2-4.5e+2).*3.025e-4,t9.*t10.*t23.*x3.*3.025e-4,t9.*t10.*t23.*x4.*3.025e-4,t9.*t23.*(-3.025e-4)];
mt4 = [x1.*(t12.*t31.*t90.*t99.*th11.*x3.*x4.*(-9.581151832460733e+3)-t12.*t31.*t80.*t90.*t99.*th11.*x3.*x4.*9.581151832460733e+5+t12.*t31.*t80.*t95.*t100.*th11.*x3.*x4.*9.581151832460733e+5+t6.*t12.*t14.*t31.*t115.*t120.*th10.*th11.*x4.*5.005005005005005e+3+t13.*t32.*t69.*t74.*t101.*th3.*th11.*x3.*x4.*1.692069709211052e+1+t12.*t31.*t69.*t102.*t124.*th3.*th11.*x3.*x4.*2.908667830133799e+6-t7.*t12.*t15.*t31.*t89.*t115.*t120.*th10.*th11.*x4.*5.005005005005005e+5+t7.*t12.*t15.*t31.*t89.*t116.*t121.*th10.*th11.*x4.*5.005005005005005e+5)];
mt5 = [-x1.*(t12.*t31.*t90.*t99.*th11.*x3.*x4.*2.908667830133799e+4+t13.*t32.*t69.*t74.*t101.*th11.*x3.*x4.*1.692069709211052e+1+t12.*t31.*t80.*t90.*t99.*th11.*x3.*x4.*2.908667830133799e+6-t12.*t31.*t80.*t95.*t100.*th11.*x3.*x4.*2.908667830133799e+6+t12.*t31.*t69.*t102.*t124.*th11.*x3.*x4.*2.908667830133799e+6)];
mt6 = [-x1.*(t12.*t31.*t90.*t99.*th11.*x3.*x4.*2.813845258871437e+4+t6.*t12.*t14.*t31.*t115.*t120.*th11.*x4.*5.005005005005005e+3+t12.*t31.*t80.*t90.*t99.*th11.*x3.*x4.*2.813845258871437e+6-t12.*t31.*t80.*t95.*t100.*th11.*x3.*x4.*2.813845258871437e+6-t7.*t12.*t15.*t31.*t89.*t115.*t120.*th11.*x4.*5.005005005005005e+5+t7.*t12.*t15.*t31.*t89.*t116.*t121.*th11.*x4.*5.005005005005005e+5)];
mt7 = [-x4.*(t12.*t31.*t90.*t99.*th11.*x3.*x4.*(-2.908667830133799)+t13.*t32.*t69.*t74.*t101.*th11.*x3.*x4.*5.817335660267597-t12.*t31.*t80.*t90.*t99.*th11.*x3.*x4.*2.908667830133799e+2+t12.*t31.*t80.*t95.*t100.*th11.*x3.*x4.*2.908667830133799e+2+t12.*t31.*t69.*t102.*t124.*th11.*x3.*x4.*1.0e+6+t6.*t12.*t14.*t20.*t31.*t115.*t120.*th11-t7.*t12.*t15.*t31.*t89.*t115.*t120.*th11.*x4.*1.0e+6+t7.*t12.*t15.*t31.*t89.*t116.*t121.*th11.*x4.*1.0e+6),0.0,x1.*(t6.*t16.*t83.*t105.*t111.*th10.*-1.0e+2+t6.*t16.*t83.*t106.*t112.*th10.*1.0e+2+t15.*t78.*t105.*t111.*th10.*x3),0.0];
mt8 = [-x1.*(t6.*t16.*t83.*t105.*t111.*-1.0e+2+t6.*t16.*t83.*t106.*t112.*1.0e+2+t15.*t78.*t105.*t111.*x3),-x4.*(t6.*t16.*t83.*t105.*t111.*(-9.99e+2./5.0)+t6.*t16.*t83.*t106.*t112.*(9.99e+2./5.0)+t15.*t78.*t105.*t111.*x3.*(9.99e+2./5.0e+2)),0.0,-t14.*t78.*t105.*t111.*x1.*x3,0.0,0.0,0.0,0.0];
mt9 = [x1.*(t86.*t125.*th3-t11.*t29.*t77.*t84.*th5.*x4.*(1.83e+2./1.91e+2)-t11.*t29.*t73.*t77.*t84.*th5.*x4.*9.581151832460733e+1+t11.*t29.*t73.*t79.*t85.*th5.*x4.*9.581151832460733e+1+t11.*t29.*t60.*t87.*t123.*th3.*th5.*x4.*2.908667830133799e+2+t11.*t14.*t29.*t105.*t111.*th5.*th10.*x3.*x4.*(5.0e+2./9.99e+2)-t6.*t11.*t15.*t29.*t78.*t105.*t111.*th5.*th10.*x4.*5.005005005005005e+1+t6.*t11.*t15.*t29.*t78.*t106.*t112.*th5.*th10.*x4.*5.005005005005005e+1)];
mt10 = [-x1.*(t11.*t31.*t90.*t99.*th5.*x4.*2.908667830133799e+4+t12.*t32.*t69.*t74.*t101.*th5.*x4.*1.692069709211052e+1+t11.*t31.*t80.*t90.*t99.*th5.*x4.*2.908667830133799e+6-t11.*t31.*t80.*t95.*t100.*th5.*x4.*2.908667830133799e+6+t11.*t31.*t69.*t102.*t124.*th5.*x4.*2.908667830133799e+6)];
mt11 = [-x1.*(t11.*t29.*t77.*t84.*th5.*x4.*2.813845258871437+t11.*t29.*t73.*t77.*t84.*th5.*x4.*2.813845258871437e+2-t11.*t29.*t73.*t79.*t85.*th5.*x4.*2.813845258871437e+2+t11.*t14.*t29.*t105.*t111.*th5.*x3.*x4.*(5.0e+2./9.99e+2)-t6.*t11.*t15.*t29.*t78.*t105.*t111.*th5.*x4.*5.005005005005005e+1+t6.*t11.*t15.*t29.*t78.*t106.*t112.*th5.*x4.*5.005005005005005e+1)];
mt12 = [-x4.*(t86.*t125.*3.438e-1-(t11.*t29.*t77.*t84.*th5.*x4)./3.438e+3-t11.*t29.*t73.*t77.*t84.*th5.*x4.*2.908667830133799e-2+t11.*t29.*t73.*t79.*t85.*th5.*x4.*2.908667830133799e-2+t11.*t29.*t60.*t87.*t123.*th5.*x4.*1.0e+2+t11.*t14.*t29.*t105.*t111.*th5.*x3.*x4-t6.*t11.*t15.*t29.*t78.*t105.*t111.*th5.*x4.*1.0e+2+t6.*t11.*t15.*t29.*t78.*t106.*t112.*th5.*x4.*1.0e+2),0.0];
dfdthSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12],5,6);
