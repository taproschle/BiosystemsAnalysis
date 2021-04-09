function fn_celstate(T,Y,data)
%Emmanuel Anane
%06/01/2019

%Calculation of different specific uptake and growth rates in the course of
%the ferementation. Time-based variation of metabolic states..


%% State variables

% X       = Y(1);           %biomass
S       = Y(:,2);           %substrate
A       = Y(:,3);           %acetate
DOTm    = Y(:,4);           % dissolved oxygem with probe response
P       = Y(:,5);
%% Parameters
Cs = 0.391;                  % carbon content of glucose (ref:enfors)
Cx = 0.488;                  % carbon content of biomas(e coli)  
Ko = 0.0001;
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

%% Calculate and plot 'cellular states' during the fermentation
qS      = (qSmax.*S)./(S+Ks).*exp(-(Kip.*P));

qSox    = (qS-(pAmax.*(qS./(qS+Kap)))).*(DOTm./(DOTm+Ko));

qSof    = qS-qSox;            

pA      = qSof.*Yas;      

qSan    = (qSox-qm).*Yem.*(Cx/Cs);  

qsA     = (qAmax./(1+(qS./Kis))).*(A./(A+Ksa));

qA      = pA-qsA;

mu      = (qSox-qm).*Yem + qsA.*Yxa + (qSof-pA).*Yxsof;

qO      = Yos.*(qSox-qSan)+qsA.*Yoa;

qp      = mu.*Ypx;   



Celstate = [qS,qSof,qSox,pA,qsA,qA,mu,qO,qp];
 %% plot cellular states
tf = {'Sp. subs. uptake rate','Sp. overflow flux','Sp. oxid.flux','Sp. ace. prod. rate',...
      'Sp. ace. re-assim. rate','Sp. ace.excretion rate', 'Sp. growth. rate','Sp. oxy. uptake rate','Sp. prod. formation rate'}; %figure titles
ly = {'q_S [g g^{-1} h^{-1}]','q_{Sof} [g g^{-1} h^{-1}]','q_{Sox} [g g^{-1} h^{-1}]','p_A [g g^{-1} h^{-1}]',...
      'q_{sA} [g g^{-1} h^{-1}]','q_A [g g^{-1} h^{-1}]','\mu [h^{-1}]','q_O [g g^{-1} h^{-1}]','q_P [g g^{-1} h^{-1}]'}; % y labels

  figure
  for k3 = 1:size(Celstate,2)
%        if k3 <=4
%            subplot(3,3,k3)
           [PKS,LOCS] = findpeaks(Celstate(:,k3),'MinPeakWidth',5);
           T_plot = T(LOCS);
% 
%            plot(T_plot,PKS,'LineWidth',1.5)%           
%            figure(4)
           subplot(3,3,k3)
           plot(T,Celstate(:,k3),'LineWidth',1.5)
%        end
        if k3 == 7 || k3 == 8 || k3 == 9
           xlabel('Time [h]') 
        else
            set(gca,'xtick',[])
        end
       ylabel(ly(k3));   title(tf(k3))
       set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold') 
   end

end