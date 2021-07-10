function dfdxSym = dfdxDewasme(t,in2,in3)
%DFDXDEWASME
%    DFDXSYM = DFDXDEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    06-Jul-2021 16:32:18

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th5 = in3(4,:);
th7 = in3(5,:);
th10 = in3(6,:);
th11 = in3(7,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th1+x2;
t3 = th5+x3;
t4 = th7+x3;
t5 = x3.^2;
t6 = 1.0./th3;
t7 = 1.0./x5;
t9 = th2.*x2.*1.0e+2;
t18 = th5.*1.719e+3;
t19 = x3.*1.719e+3;
t21 = x4.*1.0e+4;
t22 = t.*(1.1e+1./1.0e+2);
t25 = x4+1.0e-4;
t28 = th5.*x4.*1.719e+7;
t29 = x3.*x4.*1.719e+7;
t38 = th5.*th11.*x4.*5.0e+9;
t8 = t7.^2;
t10 = 1.0./t2;
t13 = 1.0./t3;
t16 = 1.0./t4;
t20 = -t9;
t23 = t21+1.0;
t24 = exp(t22);
t30 = 1.0./t25;
t39 = -t38;
t11 = t10.^2;
t12 = t10.^3;
t14 = t13.^2;
t15 = t13.^3;
t17 = t16.^2;
t26 = t10.*th2.*x2;
t27 = t2+t20;
t31 = 1.0./t23;
t34 = t9.*t10;
t42 = t6.*t7.*t24.*3.025e-4;
t44 = t13.*t30.*th5.*th11.*x4.*(5.0e+2./9.99e+2);
t46 = t13.*t30.*th5.*th11.*x4.*2.908667830133799;
t48 = t13.*t30.*th5.*th11.*x4.*2.908667830133799e+2;
t51 = t13.*t30.*th5.*th11.*x4.*5.817335660267597e+2;
t71 = t3.*1.719e+3+t28+t29+t39;
t32 = t31.^2;
t33 = t31.^3;
t35 = t26.*2.0e+2;
t36 = t26.*-1.0e+2;
t40 = t26.*1.720720720720721e-1;
t43 = -t42;
t45 = -t44;
t47 = -t46;
t49 = t13.*t31.*th5.*th11.*x4.*5.005005005005005e+3;
t50 = -t48;
t53 = t13.*t31.*th5.*th11.*x4.*2.908667830133799e+4;
t55 = -t51;
t57 = t13.*t31.*th5.*th11.*x4.*2.908667830133799e+6;
t59 = t13.*t31.*th5.*th11.*x4.*5.817335660267597e+6;
t62 = t14.*t16.*t30.*th5.*th11.*x3.*x4.*5.005005005005005e+1;
t64 = t14.*t16.*t31.*th5.*th11.*x3.*x4.*5.005005005005005e+5;
t37 = exp(t36);
t52 = -t49;
t54 = exp(t50);
t56 = -t53;
t58 = -t57;
t60 = -t59;
t63 = t26+t47;
t66 = t34+t50;
t68 = t35+t55;
t70 = t40+t45;
t41 = t26.*t37;
t61 = exp(t58);
t65 = t26+t56;
t67 = exp(t66);
t69 = exp(t68);
t73 = t37+t54;
t74 = t34+t58;
t75 = t35+t60;
t81 = t40+t52;
t88 = t16.*t70.*1.0e+2;
t90 = t17.*t70.*x3.*1.0e+2;
t91 = t16.*t70.*x3.*2.0e+2;
t92 = t16.*t70.*x3.*-1.0e+2;
t107 = t46.*t54;
t72 = t67+1.0;
t76 = exp(t74);
t77 = exp(t75);
t80 = 1.0./t73;
t83 = t37+t61;
t89 = t88.*x3;
t93 = -t90;
t94 = -t91;
t95 = exp(t92);
t97 = t16.*t81.*1.0e+2;
t100 = t17.*t81.*x3.*1.0e+2;
t101 = t16.*t81.*x3.*2.0e+2;
t102 = t16.*t81.*x3.*-1.0e+2;
t113 = t53.*t61;
t114 = t41+t107;
t78 = 1.0./t72;
t82 = t76+1.0;
t86 = 1.0./t83;
t96 = exp(t94);
t98 = t95+1.0;
t99 = t97.*x3;
t103 = -t100;
t104 = -t101;
t105 = exp(t102);
t115 = t41+t113;
t116 = t62+t88+t93;
t79 = t78.^2;
t84 = 1.0./t82;
t87 = t86.^2;
t106 = exp(t104);
t108 = 1.0./t98;
t110 = t105+1.0;
t117 = t64+t97+t103;
t85 = t84.^2;
t109 = t108.^2;
t111 = 1.0./t110;
t112 = t111.^2;
mt1 = [t43+t63.*t67.*t78.*3.294e-1+t80.*t114.*th3-t16.*t70.*t95.*t108.*th10.*x3,-t80.*t114+t63.*t67.*t78,t63.*t67.*t78.*9.674e-1+t16.*t70.*t95.*t108.*x3,0.0,0.0,x1.*(t11.*t76.*t84.*th1.*th2.*3.294e-1+t11.*t65.*t76.*t84.*th1.*th2.*3.294e+1-t11.*t65.*t77.*t85.*th1.*th2.*3.294e+1+t12.*t27.*t37.*t86.*th1.*th2.*th3+t11.*t37.*t87.*t115.*th1.*th2.*th3.*1.0e+2-t11.*t16.*t105.*t111.*th1.*th2.*th10.*x3.*1.720720720720721e-1+t5.*t11.*t17.*t81.*t105.*t111.*th1.*th2.*th10.*1.720720720720721e+1-t5.*t11.*t17.*t81.*t106.*t112.*th1.*th2.*th10.*1.720720720720721e+1)];
mt2 = [t43-x1.*(-t11.*t76.*t84.*th1.*th2+t12.*t27.*t37.*t86.*th1.*th2-t11.*t65.*t76.*t84.*th1.*th2.*1.0e+2+t11.*t65.*t77.*t85.*th1.*th2.*1.0e+2+t11.*t37.*t87.*t115.*th1.*th2.*1.0e+2),x1.*(t11.*t76.*t84.*th1.*th2.*9.674e-1+t11.*t65.*t76.*t84.*th1.*th2.*9.674e+1-t11.*t65.*t77.*t85.*th1.*th2.*9.674e+1+t11.*t16.*t105.*t111.*th1.*th2.*x3.*1.720720720720721e-1-t5.*t11.*t17.*t81.*t105.*t111.*th1.*th2.*1.720720720720721e+1+t5.*t11.*t17.*t81.*t106.*t112.*th1.*th2.*1.720720720720721e+1)];
mt3 = [-x4.*((t11.*t76.*t84.*th1.*th2)./1.0e+4+t12.*t27.*t37.*t86.*th1.*th2.*3.438e-1+(t11.*t65.*t76.*t84.*th1.*th2)./1.0e+2-(t11.*t65.*t77.*t85.*th1.*th2)./1.0e+2+t11.*t37.*t87.*t115.*th1.*th2.*3.438e+1-t11.*t16.*t105.*t111.*th1.*th2.*x3.*3.438e-1+t5.*t11.*t17.*t81.*t105.*t111.*th1.*th2.*3.438e+1-t5.*t11.*t17.*t81.*t106.*t112.*th1.*th2.*3.438e+1),0.0];
mt4 = [-x1.*(t80.*th3.*(t14.*t30.*t54.*th5.*th11.*x4.*2.908667830133799-t15.*t30.^2.*t54.*th5.^2.*th11.^2.*x4.^2.*8.460348546055261e+2)+t16.*t70.*t95.*t108.*th10-t17.*t70.*t95.*t108.*th10.*x3-t16.*t70.*t95.*t108.*t116.*th10.*x3+t16.*t70.*t96.*t109.*t116.*th10.*x3-t14.*t30.*t67.*t78.*th5.*th11.*x4.*(1.83e+2./1.91e+2)-t14.*t30.*t63.*t67.*t78.*th5.*th11.*x4.*9.581151832460733e+1+t14.*t30.*t63.*t69.*t79.*th5.*th11.*x4.*9.581151832460733e+1+t14.*t30.*t54.*t80.^2.*t114.*th3.*th5.*th11.*x4.*2.908667830133799e+2+t14.*t16.*t30.*t95.*t108.*th5.*th10.*th11.*x3.*x4.*(5.0e+2./9.99e+2))];
mt5 = [x1.*(t14.*t31.*t76.*t84.*th5.*th11.*x4.*2.908667830133799e+4+t15.*t32.*t61.*t71.*t86.*th5.*th11.*x4.*1.692069709211052e+1+t14.*t31.*t65.*t76.*t84.*th5.*th11.*x4.*2.908667830133799e+6-t14.*t31.*t65.*t77.*t85.*th5.*th11.*x4.*2.908667830133799e+6+t14.*t31.*t61.*t87.*t115.*th5.*th11.*x4.*2.908667830133799e+6)];
mt6 = [t43+x1.*(t16.*t70.*t95.*t108-t17.*t70.*t95.*t108.*x3-t16.*t70.*t95.*t108.*t116.*x3+t16.*t70.*t96.*t109.*t116.*x3+t14.*t30.*t67.*t78.*th5.*th11.*x4.*2.813845258871437+t14.*t30.*t63.*t67.*t78.*th5.*th11.*x4.*2.813845258871437e+2-t14.*t30.*t63.*t69.*t79.*th5.*th11.*x4.*2.813845258871437e+2+t14.*t16.*t30.*t95.*t108.*th5.*th11.*x3.*x4.*(5.0e+2./9.99e+2))];
mt7 = [x4.*(t16.*t81.*t105.*t111.*(9.99e+2./5.0e+2)-t17.*t81.*t105.*t111.*x3.*(9.99e+2./5.0e+2)-t16.*t81.*t105.*t111.*t117.*x3.*(9.99e+2./5.0e+2)+t16.*t81.*t106.*t112.*t117.*x3.*(9.99e+2./5.0e+2)-t14.*t31.*t76.*t84.*th5.*th11.*x4.*2.908667830133799+t15.*t32.*t61.*t71.*t86.*th5.*th11.*x4.*5.817335660267597-t14.*t31.*t65.*t76.*t84.*th5.*th11.*x4.*2.908667830133799e+2+t14.*t31.*t65.*t77.*t85.*th5.*th11.*x4.*2.908667830133799e+2+t14.*t31.*t61.*t87.*t115.*th5.*th11.*x4.*1.0e+6+t14.*t16.*t21.*t31.*t105.*t111.*th5.*th11.*x3),0.0];
mt8 = [x1.*(t13.*t32.*t76.*t84.*th5.*th11.*(-9.581151832460733e+3)-t13.*t32.*t65.*t76.*t84.*th5.*th11.*9.581151832460733e+5+t13.*t32.*t65.*t77.*t85.*th5.*th11.*9.581151832460733e+5+t14.*t33.*t61.*t71.*t86.*th3.*th5.*th11.*1.692069709211052e+1+t13.*t32.*t61.*t87.*t115.*th3.*th5.*th11.*2.908667830133799e+6+t13.*t16.*t32.*t105.*t111.*th5.*th10.*th11.*x3.*5.005005005005005e+3-t5.*t13.*t17.*t32.*t81.*t105.*t111.*th5.*th10.*th11.*5.005005005005005e+5+t5.*t13.*t17.*t32.*t81.*t106.*t112.*th5.*th10.*th11.*5.005005005005005e+5)];
mt9 = [-x1.*(t13.*t32.*t76.*t84.*th5.*th11.*2.908667830133799e+4+t14.*t33.*t61.*t71.*t86.*th5.*th11.*1.692069709211052e+1+t13.*t32.*t65.*t76.*t84.*th5.*th11.*2.908667830133799e+6-t13.*t32.*t65.*t77.*t85.*th5.*th11.*2.908667830133799e+6+t13.*t32.*t61.*t87.*t115.*th5.*th11.*2.908667830133799e+6)];
mt10 = [-x1.*(t13.*t32.*t76.*t84.*th5.*th11.*2.813845258871437e+4+t13.*t32.*t65.*t76.*t84.*th5.*th11.*2.813845258871437e+6-t13.*t32.*t65.*t77.*t85.*th5.*th11.*2.813845258871437e+6+t13.*t16.*t32.*t105.*t111.*th5.*th11.*x3.*5.005005005005005e+3-t5.*t13.*t17.*t32.*t81.*t105.*t111.*th5.*th11.*5.005005005005005e+5+t5.*t13.*t17.*t32.*t81.*t106.*t112.*th5.*th11.*5.005005005005005e+5)];
mt11 = [t43-t86.*t115.*3.438e-1-x4.*(t13.*t32.*t76.*t84.*th5.*th11.*(-2.908667830133799)+t14.*t33.*t61.*t71.*t86.*th5.*th11.*5.817335660267597-t13.*t32.*t65.*t76.*t84.*th5.*th11.*2.908667830133799e+2+t13.*t32.*t65.*t77.*t85.*th5.*th11.*2.908667830133799e+2+t13.*t32.*t61.*t87.*t115.*th5.*th11.*1.0e+6+t13.*t16.*t32.*t105.*t111.*th5.*th11.*x3.*1.0e+4-t5.*t13.*t17.*t32.*t81.*t105.*t111.*th5.*th11.*1.0e+6+t5.*t13.*t17.*t32.*t81.*t106.*t112.*th5.*th11.*1.0e+6)-(t65.*t76.*t84)./1.0e+4+t16.*t81.*t105.*t111.*x3.*(9.99e+2./5.0e+2)-1.8e+4,0.0];
mt12 = [t6.*t8.*t24.*x1.*3.025e-4,t6.*t8.*t24.*(x2-4.5e+2).*3.025e-4,t6.*t8.*t24.*x3.*3.025e-4,t6.*t8.*t24.*x4.*3.025e-4,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12],5,5);