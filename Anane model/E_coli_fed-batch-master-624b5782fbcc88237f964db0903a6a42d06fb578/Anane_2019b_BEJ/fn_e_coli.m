function dydt= fn_e_coli(t, y, data)

%Emmanuel Anane, 
%Institute fur Biotechnologie, TU-Berlin.
%28/07/2016

%this function is the set of ODEs for a fed-batch culture of e coli. y(1) = dV/dt,
%y(2) = dx/dt, y(3) = dS/dt, y(4) = dP/dt, y(5) = dF/dt. The script that 
%calls it implements a fed-batch process of e. coli fermentation, with
%induction at a later stage in the fermentation.
%%
if isempty(t)       
    
    syms('V', 'X', 'S', 'A', 'DOT','DOTm', 'F',...
        'Si','DOTstar','Kla','tau','mu_set')
    
else

    V       = y(1);           %volume
    X       = y(2);           %biomass
    S       = y(3);           %substrate
    A       = y(4);           %acetate
    DOT     = y(5);         %dissolved oxygen, no probe response
    DOTm    = y(6);        %dissolved oxygen, with probe response time
    F       = y(7);           % feed
    %% INPUTS
    Si          = data.u(1);                      %feed concentration (g/L)
    mu_set      = data.u(2);
    DOTstar     = data.u(3);     
    Kla         = data.u(4);
    tau         = data.u(5);
end

%%  Fermentation constants
Cs = 0.391;%                    %carbon content of glucose (ref:enfors)
Cx = 0.488;                     %carbon content of biomas(e coli)  
H = 14000;                      %Henry's law constant to convert from mol to % DO
Ko = 0.0001;

%%
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
    Yos     = data.Par(14);  %yield of oxygen on glucose (g/g)
    Yxsof   = data.Par(15);  %biomass yield from the overflow route
end
%% Explicit algebraic equations

qS      = (qSmax./(1+(A/Kia)))*(S/(S+Ks));
qSof    = pAmax*(qS/(qS+Kap));                     %move Yas to acetate production
pA      = qSof*Yas;                 %acetate production from the overflow flux, and not from total qS
qSox    = (qS-qSof)*(DOT/(DOT+Ko)); % ;            %When there is no oxygen, only the qsof is operational            
qSan    = (qSox-qm)*Yem*Cx/Cs;      
qsA     = (qAmax./(1+(qS/Kis)))*(A/(A+Ksa));
qA      = pA-qsA;                            %net acetate accumulation
mu      = (qSox-qm)*Yem + qsA*Yxa + qSof*Yxsof;
qO      = Yos*(qSox-qSan)+qsA*Yoa;
Kp      = (1/tau)*3600;
%% ODE system
dVdt    = F;                          %dV/dt
dXdt    = X*(mu - (F/V));             %dX/dt
dSdt    = F/V*(Si - S) - (qS*X);      %dS/dt
dAdt    = qA*X-(F/V)*A;               %dA/dt
dDOTdt  = Kla*(DOTstar-DOT)-qO*X*H;   %dDOT/dt
dDOTmdt = Kp*(DOT-DOTm);              %dissolved oxygen with probe response 
    if data.flag == true                %batch phase
        dFdt = 0;                  %feed rate not changing.
    else
        dFdt = mu_set*F;            %dF/dt = feed rate changing exponentially
    end

  if ~isempty(t)
        if length(y) == 7 %normal without Jacobians
            dydt = [dVdt;dXdt;dSdt;dAdt;dDOTdt;dDOTmdt;dFdt];
        else
               if data.flag == true              %batch phase
                  ModJacX_batch
                  ModJacP_batch
               else    
                  ModJacX_fb
                  ModJacP_fb
               end
            m = length(JacX); % number of state variables
            n = length(JacP);%./length(JacX); % number of Parameters
            Sr = reshape(y(8:end),[m,n]); 
            Sm = JacX*Sr;       %relative sensitivity matrix
            SLi = reshape(Sm,m*n,1);
            JacPe = reshape(JacP,[m*n,1]);
            dydt = [dVdt;dXdt;dSdt;dAdt;dDOTdt;dDOTmdt;dFdt;JacPe+SLi];
        end
    else
      dydt = [dVdt;dXdt;dSdt;dAdt;dDOTdt;dDOTmdt;dFdt];
   end
end


