function Rank =  Rank_calculation(jacobian,CondNMaxCrit,CollinIndexMaxCrit)

%{
    1.RANK DETERMINATION 

    SVD of J=USV' (Hansen, 1998 ; Liang, 2002; Golub, 1996 & Hansen, 1998)
    

   
    %%%%% Inputs

    jacobian: jacobian of LSQ Objective Function or sensitivity matrix(dy_i/dpar_j) [Ny.Nm x Np]
    par:   estimated parameters using LSQ [1xNp] 
    CondMaxCrit: Condition Number Grah's threshold (2004)  (CondMaxCrit=1000 -default)
    CollinIndexMaxCrit: Collinearity Index Brun's threshold (2002) (CollinIndexMaxCrit=10 -default)


    %%%%%% : Output

    Numerical rank of the sensitivity matrix 
      
%}

    %% -> check if (number of rows in Jacobian--number of datapoints) < (number of parameters)

    [numMeas,numPar] = size(jacobian);
    if  numMeas < numPar
        tempJac = jacobian;
        jacobian = zeros(numPar);
        jacobian(1:numMeas,:) = tempJac;
    end

    %% 1.RANK DETERMINATION 


    J = jacobian;

    %1. SVD of J=USV'
    %V: right and left singular vectors of J are equal to eingenvectors of FIM = J'J, but these are left out of the current svd call.
    [~, S, ~] = svd(J,'econ'); % using economic form of svd to reduce memory usage.


    Sing_values = abs(diag(S)); % absolute values of the diagonal elements of S, which are the singular values
    dimSvalues  = size(Sing_values,1); %number of singular values
    max_Sval    = max(Sing_values); % maximum singular value

    
    CondN_maxCriterion       = CondNMaxCrit;
    CollinIndex_maxCriterion = CollinIndexMaxCrit;

    CondN_CollIdx_selected = [];
    CondN_Sub = zeros(dimSvalues,1);
    CollinIndex_Sub = zeros(dimSvalues,1);

    for i = 1:dimSvalues
        CondN_Sub(i) = max_Sval/Sing_values(i); % condition number is the max singular value/singular value of a parameter
        CollinIndex_Sub(i) = 1/Sing_values(i);  % Collinearity index is the inverse of the singular values
        if abs(CondN_Sub(i)) <= CondN_maxCriterion && (abs(CollinIndex_Sub(i))) <= CollinIndex_maxCriterion % meeting criteria for both sensitivity and linear independence
            CondN_CollIdx_selected  = [CondN_CollIdx_selected; CondN_Sub(i)];
        end
    end

    Rank    = length(CondN_CollIdx_selected); %J's Rank
    
%% Plot identifiability diagnosis results

eps_CondN          = Sing_values(1)/CondN_maxCriterion;
eps_CollinIndex    = 1/CollinIndex_maxCriterion;
eps_SVMin          = max(eps_CondN,eps_CollinIndex);

Num_par = length(Sing_values);
eps_SVMin_vector   = eps_SVMin.*ones(Num_par,1);

SvS1 = figure;
semilogy1 = semilogy(Sing_values,'Color',[0 0 0]);
title('Singular value spectrum')
xlabel('Index'); 
ylabel('Singular value ')    
set(semilogy1,'Marker','s');
set(semilogy1,'MarkerEdgeColor','k');
set(semilogy1,'MarkerFaceColor','g');
set(semilogy1,'LineWidth',2,'Color','b');
grid('on');
xlim([1 Num_par])
hold on
plot1 = semilogy(eps_SVMin_vector,'Color',[0 0 0]);
set(plot1,'DisplayName','Thrsh');
set(plot1,'LineStyle','-.','LineWidth',2,'Color','k');         
hold off
legend('Singular values',' Threshold [\gamma_{max}, \kappa_{max}]');
legend('boxoff')
set(gca,'FontSize', 22, 'FontWeight','bold')

% saveas(SvS1,[pwd '/Results/Singular value spectrum' datestr(now) '.fig']);
% saveas(SvS1,[pwd '/Results/Singular value spectrum' datestr(now) '.tif']);

end