function dydt= fn_Ecoli(t, y, data)


%Emmanuel Anane
%Chair of Bioprocess Enginering,
%Institute fur Biotechnologie, TU-Berlin.
%2/07/2018

%This function is the set of ODEs for a fed-batch culture of e coli. The script that 
%calls it implements a fed-batch process of recombinant E. coli fermentation
%with recombinant protein production. The model is also used to calculate
%the sensitivities of outputs wrt parameters using the forward sensitivity
%method.
%% for generating jacobians
if isempty(t)       
    
    syms('V', 'X', 'S', 'A', 'DOT','F', 'P',...
        'Si','mu_set')
    
else

    %% State variables

    X       = y(1);           %biomass
    S       = y(2);           %substrate
    A       = y(3);           %acetate
    DOTm    = y(4);           % dissolved oxygem with probe response
    P       = y(5);
    GR      = y(6);
    
    %% Experimental INPUTS
    Si          = data.u(1);                     
    V           = data.u(2);           %volume
    F           = data.u(3);
    pH          = data.u(5);

end
%%  Fermentation constants
Cs      = 0.391;%             % carbon content of glucose (ref:enfors)
Cx      = 0.488;              % carbon content of biomas(e coli)  
H       = 14000;              % Henry's law constant to convert dissolved oxygen from mol to % DO
Ko      = 0.0001;
DOTstar = 98.9;
tau     = 26/3600;            %from Kla determination, tau = 26seconds.
%% Glucose release by Enbase

% The glucose release equation is valid for Enbase defined medium "Enpresso B
% defined", where 1 tablet is dissolved in 25 ml distilled water to make the
% growth medium. Operating pH should be in the range 6-7.5. 

% GR =  28.6542698856704 * pH + 3.26024668806902 * U - 0.351712 * pH * U - 2.1029093224566 * pH^2 - 0.04190619499737 * U^2 - 99.95486283157305;

% GR = 0;
%% Parameters

if isempty(t)
    syms('Kap', 'Ksa', 'Ks', 'Kis', 'Kla','pAmax', 'qAmax', 'qm','qSmax',...
    'Yas', 'Yxa', 'Yem', 'Yxsof','Ypx', 'Yoa', 'Yos', 'Kip')
else

    Kap     = data.Par(1);   %monod-type saturation constant, intracellular acetate production
    Ksa     = data.Par(2);   %affinity constant, acetate consumption (0.05 g/L)
    Ks      = data.Par(3);   %affinity constant, substrate consumption (0.05 gglu./L)
    Kip     = data.Par(4);
    Kis     = data.Par(5);   %inhibition constant, inhib. of ace.uptake by glucose(4g/L)
    pAmax   = data.Par(6);   %max spec acetate production rate (0.15 g_ace/(gx.h)))
    qAmax   = data.Par(7);   %max spec acetate consumption rate (0.15 g_ace/(gx.h)))
    qm      = data.Par(8);   %spec maintenance coefficient (0.04 g_glu/(gx.h))
    qSmax   = data.Par(9);   %max spec glucose uptake rate (1.3 g_glu./gx.h)
    Yas     = data.Par(10);  %yield of intracellular acetate on qsof (g.ace/g.glu) intracellular. How much of the total overflow flux is converted to acetate. The remaining goes to PPP and mixed acid pathways.
    Yxa     = data.Par(11);  %yield of biomass on acetate
    Yem     = data.Par(12);  %yield exclusive maintenance (Yem = Yxs-qm)
    Yxsof   = data.Par(13);  %biomass yield from the overflow route
    Ypx     = data.Par(14);
    Yoa     = data.Par(15);
    Yos     = data.Par(16);
    U1      = data.Par(17);
    U2      = data.Par(18);
    
    switch data.change
%         case 1
%             Kla     = data.Par(13);
%         case 2
%             Kla     = data.Par(14);
%         case 3
%             Kla     = data.Par(15);
%         case 4
%             Kla     = data.Par(16);
        case 5
            Kla     = data.Par(19);
        case 6
            Kla     = data.Par(20);

    end
end
%% Explicit algebraic equations (See text in paper for explanations)

qS      = (qSmax.*S)./(S+Ks).*exp(-(Kip*P));

qSox    = (qS-(pAmax*(qS/(qS+Kap)))).*(DOTm/(DOTm+Ko));

qSof    = qS-qSox;            

pA      = qSof*Yas;      

qSan    = (qSox-qm)*Yem*Cx/Cs;  

qsA     = (qAmax./(1+(qS/Kis)))*(A/(A+Ksa));

qA      = pA-qsA;

mu      = (qSox-qm)*Yem + qsA*Yxa + (qSof-pA)*Yxsof;

qO      = Yos*(qSox-qSan)+qsA*Yoa;

qp      = mu.*Ypx;   

Kp      = 1./tau;

DOT     = (Kla*DOTstar-qO*X*H)./Kla;

%% ODEs
if data.flag == 0 % batch phase
       dPdt    = 0;
       dGRdt =  0;

elseif data.flag == 1                    %exponential fed-batch phase
        dPdt    = 0;
        dGRdt   =  U1*GR;

elseif data.flag == 2;           %constant feed fed-batch phase, protein induction
        dPdt    = qp-mu*P;
        GR      = 28.6542698856704 * pH + 3.26024668806902 * U2 - 0.351712 * pH * U2 - 2.1029093224566 * pH^2 - 0.04190619499737 * U2^2 - 99.95486283157305;
        dGRdt   = 0;
end

dXdt        = X*(mu - (F/V));             %
dSdt        = F/V*(Si - S) + GR - (qS*X);      %dS/dt
dAdt        = qA*X-(F/V)*A;               %dA/dt
dDOTmdt     = Kp.*(DOT-DOTm);
    

if ~isempty(t)
    if length(y) == 6 %normal solution of ODE without Senstitivities
        dydt = [dXdt; dSdt; dAdt;dDOTmdt;dPdt;dGRdt];    
    else
        if data.flag == 1
            ModJacX_exp
            ModJacP_exp
        else
            ModJacX_const
            ModJacP_const
        end
        m = length(JacX); % number of state variables
        n = length(JacP);%./length(JacX); % number of Parameters
        Sr = reshape(y(8:end),[m,n]); 
        Sm = JacX*Sr;       %relative sensitivity matrix
        SLi = reshape(Sm,m*n,1);
        JacPe = reshape(JacP,[m*n,1]);
        dydt = [dXdt;dSdt;dAdt;dDOTmdt;dPdt;JacPe+SLi];
    end
else
    
            dydt = [dXdt; dSdt; dAdt;dDOTmdt;dPdt];    %for generating jacobians if t is empty
end

