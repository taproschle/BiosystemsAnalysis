function dydt= fn_Ecoli(t, y, data)
%{ 
Emmanuel Anane 
Institute fur Biotechnologie, TU-Berlin.
09/11/2016

this function is the set of ODEs for a fed-batch culture of e coli. y(1) = dV/dt,
y(2) = dx/dt, y(3) = dS/dt, y(4) = dP/dt, y(5) = dF/dt. The script that 
calls it implements a fed-batch process of E. coli fermentation, with
product formation after iptg induction.
%}
%% For calculating Jacobians

if isempty(t)       
    
    syms('V', 'X', 'S', 'A', 'DOTm', 'F','Cs', 'Cx', 'H',...
        'Si','DOTstar','Kla','tau')
    
    Cs          = 0.391;                   
    Cx          = 0.488;                     
    H           = 14000;                
else
    %% State variables
    X       = y(1);           %biomass
    S       = y(2);           %substrate
    A       = y(3);           %acetate
    DOTm    = y(4);           %dissolved oxygen, no probe response
    P       = y(5);           %dissolved oxygen, with probe response time           %feed
    %%  Fermentation constants
    Cs  = 0.391;               %carbon content of glucose (ref:enfors)
    Cx  = 0.488;               %carbon content of biomas(e coli)  
    H   = 14000;               %Henry's law constant to convert from mol to % DO
    %% INPUTS
    Si          = data.u(1);        %feed concentration (g/L)
    DOTstar     = data.u(2);     
    Kla         = data.u(3);
    tau         = data.u(4);
    F           = data.u(6);
    V           = data.u(7);


end
%%  Parameters

if isempty(t)
    syms('Kap', 'Ksa', 'Ks', 'Kia', 'Kis', 'pAmax', 'qAmax',...
        'qm', 'qSmax', 'Yas', 'Yoa', 'Yxa', 'Yem', 'Yos', 'Yxsof')
else
    
    Kap     = data.Par(1);   %monod-type saturation constant, intracellular acetate production
    Ksa     = data.Par(2);   %affinity constant, acetate consumption (0.05 g/L)
    Ks      = data.Par(3);   %affinity constant, substrate consumption (0.05 gglu./L)
    Kia     = data.Par(4);   %inhibition constant, inhib. of glucose uptake by acetate
    Kis     = data.Par(5);   %inhibition constant, inhib. of ace.uptake by glucose(4g/L)
    pAmax   = data.Par(6);   %max spec acetate production rate (0.15 g_ace/(gx.h)))
    qAmax   = data.Par(7);   %max spec acetate consumption rate (0.15 g_ace/(gx.h)))
    qm      = data.Par(8);   %spec maintenance coefficient (0.04 g_glu/(gx.h))
    qSmax   = data.Par(9);   %max spec glucose uptake rate (1.3 g_glu./gx.h)
    Yas     = data.Par(10);  %yield of acetate on substrate (g.ace/g.glu)
    Yoa     = data.Par(11);  %yield of oxygen on acetate g/g)
    Yxa     = data.Par(12);  %yield of biomass on acetate
    Yem     = data.Par(13);  %yield exclusive maintenance (Yem = Yxs-qm)
    Yxsof   = data.Par(14);  %yield of biomass on other overflow products, such as formate (g/g)
    Yos     = data.Par(15);  %yield of oxygen on glucose (g/g)
    Ypx     = data.Par(16);  %yield of product on biomass (g/g)
end
%% Explicit algebraic equations
qS      = (qSmax./(1+(A/Kia)))*(S/(S+Ks));
qSof    = pAmax*(qS/(qS+Kap));
pA      = qSof*Yas;                             %acetate production from the overflow flux, and not from total qS
qSox    = (qS-qSof);                            %When there is no oxygen, only the qsof is operational            
qSan    = (qSox-qm)*Yem*Cx/Cs;      
qsA     = (qAmax./(1+(qS/Kis)))*(A/(A+Ksa));    %acetate consumption, what happens if A is replaced with pA?
qA      = pA-qsA;                               %net acetate accumulation
mu      = (qSox-qm)*Yem + qsA*Yxa + qSof*Yxsof;
qO      = Yos*(qSox-qSan)+qsA*Yoa;
Kp      = (1/tau)*3600;
DOT     = (Kla*DOTstar-qO*X*H)./Kla;
qp      = mu.*Ypx;

%% ODE system
dXdt    = X*(mu - (F/V));             %dX/dt
dSdt    = F/V*(Si - S) - (qS*X);      %dS/dt
dAdt    = qA*X-(F/V)*A;               %dA/dt
dDOTmdt = Kp*(DOT-DOTm);              %dissolved oxygen with probe response 

    if data.flag == 1              %batch phase
      dPdt = 0;
      elseif data.flag == 2     %exponential fed-batch phase, non-production.
      dPdt = qp-mu*P;
    end
            
if ~isempty(t)
        if length(y) == 5 %normal without Jacobians
            dydt = [dXdt;dSdt;dAdt;dDOTmdt;dPdt];
        else
            m = length(JacX); % number of state variables
            n = length(JacP);%./length(JacX); % number of Parameters
            Sr = reshape(y(7:end),[m,n]); 
            Sm = JacX*Sr;       %relative sensitivity matrix
            SLi = reshape(Sm,m*n,1);
            JacPe = reshape(JacP,[m*n,1]);
            dydt = [dXdt;dSdt;dAdt;dDOTmdt;dPdt;JacPe+SLi];

        end
    else
      dydt = [dXdt;dSdt;dAdt;dDOTmdt;dPdt];
end
end


