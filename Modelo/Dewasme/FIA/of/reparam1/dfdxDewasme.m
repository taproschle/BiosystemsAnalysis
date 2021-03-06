function dfdxSym = dfdxDewasme(t,in2,in3)
%DFDXDEWASME
%    DFDXSYM = DFDXDEWASME(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    09-Jul-2021 18:50:51

th1 = in3(1,:);
th2 = in3(2,:);
th3 = in3(3,:);
th4 = in3(4,:);
th5 = in3(5,:);
th6 = in3(6,:);
th8 = in3(7,:);
th9 = in3(8,:);
th10 = in3(9,:);
th11 = in3(10,:);
th12 = in3(11,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
t2 = th1+x2;
t3 = th5+x3;
t4 = th4.^2;
t5 = th5.^2;
t6 = th11.^2;
t7 = x3.^2;
t8 = x4.^2;
t9 = 1.0./th3;
t10 = 1.0./th4;
t12 = 1.0./th9;
t14 = 1.0./x5;
t22 = x4+1.0e-4;
t25 = t.*(1.97e+2./1.0e+3);
t26 = x3+5.1e-3;
t11 = 1.0./t4;
t13 = t12.^2;
t15 = t14.^2;
t16 = 1.0./t2;
t18 = 1.0./t3;
t27 = exp(t25);
t30 = 1.0./t22;
t35 = 1.0./t26;
t17 = t16.^2;
t19 = t18.^2;
t20 = t18.^3;
t21 = t16.*th2;
t31 = t30.^2;
t36 = t35.^2;
t46 = t10.*t18.*t30.*th5.*th11;
t59 = t9.*t14.*t27.*5.4175e-4;
t23 = t21.*x2;
t24 = t17.*th2.*x2;
t28 = t21.*8.0e+1;
t47 = t46.*x4;
t48 = t10.*t18.*t31.*th5.*th11.*x4;
t51 = t46.*8.0e+1;
t60 = -t59;
t62 = t12.*t19.*t30.*t35.*th5.*th11.*x3.*x4.*8.0e+1;
t29 = -t24;
t32 = t23.*8.0e+1;
t33 = t24.*8.0e+1;
t34 = t23.*1.6e+2;
t49 = -t47;
t50 = -t48;
t52 = t47.*8.0e+1;
t53 = t48.*8.0e+1;
t54 = t47.*1.6e+2;
t37 = -t32;
t38 = -t33;
t40 = t21+t29;
t55 = -t52;
t56 = -t53;
t57 = -t54;
t61 = t23+t49;
t74 = t46+t50;
t39 = exp(t37);
t45 = t28+t38;
t58 = exp(t55);
t63 = t32+t55;
t64 = t34+t57;
t80 = t51+t56;
t81 = t12.*t35.*t61.*th4.*8.0e+1;
t83 = t12.*t36.*t61.*th4.*x3.*8.0e+1;
t84 = t12.*t35.*t61.*th4.*x3.*1.6e+2;
t85 = t12.*t35.*t61.*th4.*x3.*-8.0e+1;
t41 = t21.*t39;
t42 = t23.*t39;
t43 = t24.*t39;
t65 = exp(t63);
t66 = exp(t64);
t68 = t39+t58;
t75 = t46.*t58;
t76 = t47.*t58;
t77 = t48.*t58;
t78 = t10.*t19.*t30.*t58.*th5.*th11.*x4;
t82 = t81.*x3;
t86 = -t83;
t87 = -t84;
t88 = exp(t85);
t91 = t5.*t6.*t8.*t11.*t20.*t31.*t58.*8.0e+1;
t44 = -t41;
t67 = t65+1.0;
t69 = t42.*t45;
t72 = 1.0./t68;
t79 = -t75;
t89 = exp(t87);
t90 = t88+1.0;
t92 = -t91;
t95 = t42+t76;
t97 = t76.*t80;
t99 = t62+t81+t86;
t70 = 1.0./t67;
t73 = t72.^2;
t93 = 1.0./t90;
t96 = t43+t44+t69;
t98 = t78+t92;
t100 = t77+t79+t97;
t71 = t70.^2;
t94 = t93.^2;
mt1 = [t60+t72.*t95.*th3+t61.*t65.*t70.*th8-t12.*t35.*t61.*t88.*t93.*th4.*th10.*x3,-t72.*t95+t61.*t65.*t70,t61.*t65.*t70.*th6+t12.*t35.*t61.*t88.*t93.*th4.*x3,0.0,0.0,x1.*(t72.*th3.*(t41-t69+t29.*t39)+t40.*t65.*t70.*th8+t45.*t61.*t65.*t70.*th8-t45.*t61.*t66.*t71.*th8+t39.*t45.*t73.*t95.*th3-t12.*t35.*t40.*t88.*t93.*th4.*th10.*x3+t4.*t7.*t13.*t36.*t40.*t61.*t88.*t93.*th10.*8.0e+1-t4.*t7.*t13.*t36.*t40.*t61.*t89.*t94.*th10.*8.0e+1),t60-x1.*(t72.*(t41-t69+t29.*t39)-t40.*t65.*t70-t45.*t61.*t65.*t70+t45.*t61.*t66.*t71+t39.*t45.*t73.*t95)];
mt2 = [x1.*(t40.*t65.*t70.*th6+t45.*t61.*t65.*t70.*th6-t45.*t61.*t66.*t71.*th6+t12.*t35.*t40.*t88.*t93.*th4.*x3-t4.*t7.*t13.*t36.*t40.*t61.*t88.*t93.*8.0e+1+t4.*t7.*t13.*t36.*t40.*t61.*t89.*t94.*8.0e+1),-x4.*(t72.*th4.*(t41-t69+t29.*t39)+t40.*t65.*t70.*th12+t45.*t61.*t65.*t70.*th12-t45.*t61.*t66.*t71.*th12+t39.*t45.*t73.*t95.*th4-t35.*t40.*t88.*t93.*th4.*x3+t4.*t7.*t12.*t36.*t40.*t61.*t88.*t93.*8.0e+1-t4.*t7.*t12.*t36.*t40.*t61.*t89.*t94.*8.0e+1),0.0];
mt3 = [-x1.*(t72.*t98.*th3+t73.*t78.*t95.*th3.*8.0e+1+t12.*t35.*t61.*t88.*t93.*th4.*th10-t12.*t36.*t61.*t88.*t93.*th4.*th10.*x3-t12.*t35.*t61.*t88.*t93.*t99.*th4.*th10.*x3+t12.*t35.*t61.*t89.*t94.*t99.*th4.*th10.*x3-t10.*t19.*t30.*t65.*t70.*th5.*th8.*th11.*x4-t10.*t19.*t30.*t61.*t65.*t70.*th5.*th8.*th11.*x4.*8.0e+1+t10.*t19.*t30.*t61.*t66.*t71.*th5.*th8.*th11.*x4.*8.0e+1+t12.*t19.*t30.*t35.*t88.*t93.*th5.*th10.*th11.*x3.*x4),x1.*(t72.*t98+t73.*t78.*t95.*8.0e+1+t10.*t19.*t30.*t65.*t70.*th5.*th11.*x4+t10.*t19.*t30.*t61.*t65.*t70.*th5.*th11.*x4.*8.0e+1-t10.*t19.*t30.*t61.*t66.*t71.*th5.*th11.*x4.*8.0e+1)];
mt4 = [t60+x1.*(t12.*t35.*t61.*t88.*t93.*th4-t12.*t36.*t61.*t88.*t93.*th4.*x3-t12.*t35.*t61.*t88.*t93.*t99.*th4.*x3+t12.*t35.*t61.*t89.*t94.*t99.*th4.*x3+t10.*t19.*t30.*t65.*t70.*th5.*th6.*th11.*x4+t10.*t19.*t30.*t61.*t65.*t70.*th5.*th6.*th11.*x4.*8.0e+1-t10.*t19.*t30.*t61.*t66.*t71.*th5.*th6.*th11.*x4.*8.0e+1+t12.*t19.*t30.*t35.*t88.*t93.*th5.*th11.*x3.*x4),x4.*(t72.*t98.*th4+t35.*t61.*t88.*t93.*th4-t36.*t61.*t88.*t93.*th4.*x3-t35.*t61.*t88.*t93.*t99.*th4.*x3+t35.*t61.*t89.*t94.*t99.*th4.*x3+t19.*t30.*t58.*t73.*t95.*th5.*th11.*x4.*8.0e+1-t10.*t19.*t30.*t65.*t70.*th5.*th11.*th12.*x4+t19.*t30.*t35.*t88.*t93.*th5.*th11.*x3.*x4-t10.*t19.*t30.*t61.*t65.*t70.*th5.*th11.*th12.*x4.*8.0e+1+t10.*t19.*t30.*t61.*t66.*t71.*th5.*th11.*th12.*x4.*8.0e+1),0.0];
mt5 = [-x1.*(t72.*t100.*th3+t65.*t70.*t74.*th8+t61.*t65.*t70.*t80.*th8-t61.*t66.*t71.*t80.*th8-t58.*t73.*t80.*t95.*th3-t12.*t35.*t74.*t88.*t93.*th4.*th10.*x3+t4.*t7.*t13.*t36.*t61.*t74.*t88.*t93.*th10.*8.0e+1-t4.*t7.*t13.*t36.*t61.*t74.*t89.*t94.*th10.*8.0e+1),-x1.*(-t72.*t100+t65.*t70.*t74+t61.*t65.*t70.*t80-t61.*t66.*t71.*t80+t58.*t73.*t80.*t95),-x1.*(t65.*t70.*t74.*th6+t61.*t65.*t70.*t80.*th6-t61.*t66.*t71.*t80.*th6+t12.*t35.*t74.*t88.*t93.*th4.*x3-t4.*t7.*t13.*t36.*t61.*t74.*t88.*t93.*8.0e+1+t4.*t7.*t13.*t36.*t61.*t74.*t89.*t94.*8.0e+1)];
mt6 = [t60+x4.*(t72.*t100.*th4+t65.*t70.*t74.*th12+t61.*t65.*t70.*t80.*th12-t61.*t66.*t71.*t80.*th12-t58.*t73.*t80.*t95.*th4-t35.*t74.*t88.*t93.*th4.*x3+t4.*t7.*t12.*t36.*t61.*t74.*t88.*t93.*8.0e+1-t4.*t7.*t12.*t36.*t61.*t74.*t89.*t94.*8.0e+1)-t72.*t95.*th4-t61.*t65.*t70.*th12+t35.*t61.*t88.*t93.*th4.*x3-1.8e+4,0.0,t9.*t15.*t27.*x1.*5.4175e-4,t9.*t15.*t27.*(x2-4.5e+2).*5.4175e-4,t9.*t15.*t27.*x3.*5.4175e-4,t9.*t15.*t27.*x4.*5.4175e-4,0.0];
dfdxSym = reshape([mt1,mt2,mt3,mt4,mt5,mt6],5,5);
