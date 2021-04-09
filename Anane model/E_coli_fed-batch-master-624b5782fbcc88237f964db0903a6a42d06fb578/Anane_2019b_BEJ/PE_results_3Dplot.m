% Eduardo Vieyra Olguin
% Diana C. Lopez C.
% Emmanuel Anane
function PE_results_3Dplot(par_All,par_label,alpha,Ranking)

%{
par_All::         Matrix containing the MonteCarlo results (L executions) of the Parameters (Np x L)
par_label::       Parameters label vector
alpha::           Significance level for the t-student confidence interval
%}
%% Count the Number of Parameters


% Definition of the number of parameters and Monte Carlo repetitions 
    [nPar,repeat]=size(par_All);
    
    parInFig = nPar;

% points:: Number of Points in the Interval for the Distribution Fiting.
    points=40;
    
% Computing the mean (\mu) and standard deviation (\sigma) 
% of the parameter distribution assuming a normal distribution 

    counter=0;
    Nor_mu = zeros(nPar,1);
    Nor_sigma = zeros(nPar,1);
    par_label_vec = cell(nPar,1);

    for i=1:nPar;
         counter=1+counter;
         X=par_All(i,:);
         Nor_mu(counter)=sum(X)/repeat;
         Suma=0;
         for h=1:repeat;
            Suma=(X(h)-Nor_mu(counter))^2 + Suma;
         end
         Nor_sigma(counter)=sqrt(Suma/repeat);
         par_label_vec(counter)=par_label(i);
     end
  
%% Estimator assessment summary    
% Estimator performance assessemnt summary in TableNormalDist, 
    
    faktor=2.5;
    
    [TableNormalDist,Ind,graph_low,graph_up]=tabnorm(par_label_vec,Nor_mu,Nor_sigma,alpha,counter,repeat,Ranking,faktor);

% Code for saving the Table with mean, standard deviation, and the lower and
% upper limits of the confidence interval; the result file will have the
% name TableNormalDist.
   
%     save([pfad ['TableNormalDist_' CaseName]],'TableNormalDist')

%% Creation of the Distribution Fit
%Call for createFit function that do the fitting with normal and
%non-parametric distribution with outputs: NormFit is the Dataset for
%ploting the Normal Distribution and Nonparfit the Dataset for the Non
%parametric distribution.
    [NormFit,NonparFit] = createFit(par_All,points,nPar,graph_low,graph_up);
    
%% Preparation for the Plot
    
    counter=0;
    for i=1:nPar;
        counter=counter+1;
    
%         tetaTrueNormVec(counter)=tetaTrueNorm(i);
        FitPlotY(:,counter)=NormFit.Fit{1,i};
        FitPlotX(:,counter)=NormFit.X{1,i};
        NonparFitPlotY(:,counter)=NonparFit.Fit{1,i};
        NonparFitPlotX(:,counter)=NonparFit.X{1,i};       
    end
%Label Vector for the Plot
    par_labelNew=TableNormalDist(2:end,1);
        
%mu for the plot
    mu_plotNew=TableNormalDist(2:end,3);
    mu_plotNew=cell2mat(mu_plotNew);
    
    k=0;
    for i=1:counter;
        
        if Ranking==1;
            R=Ind(i);
        else
            R=i;    
        end
%Reducing the Data to vectors depending on how many Figures are going
%to be ploted --------------Important just for the plot-----------------
         Plot_Vec_NonParX(:,i-(k*parInFig))=NonparFitPlotX(:,R);
         Plot_Vec_NonParY(:,i-(k*parInFig))=NonparFitPlotY(:,R);
 
         par_label_vec_plot(i-(k*parInFig))=par_labelNew(i,1);
    
         mu_plot(1,i-(k*parInFig))=mu_plotNew(i,1);
        
%          tetaTrueNorm_select(1,i-(k*parInFig))=tetaTrueNormVec(1,R);

         Plot_Vec_FitX(:,i-(k*parInFig))=FitPlotX(:,R);
         Plot_Vec_FitY(:,i-(k*parInFig))=FitPlotY(:,R);
         
    if mod(i,parInFig) == 0;

             k=k+1;
%% Pdf_plot Call
% Call of the function pdf_plot wich has not outputs but the figures and
% has as input: Vectors with the size of variable graf for the
% nonparametric and normal distribution, the points to be ploted as vector
% X, labels vector, the path and experiment number to save the figures, 
% the mean and true value vector.
            pdf_plot(Plot_Vec_NonParX,Plot_Vec_NonParY,Plot_Vec_FitX,Plot_Vec_FitY,par_label_vec_plot,mu_plot,k)

             Plot_Vec_FitX=zeros(points,parInFig);
             Plot_Vec_FitY=zeros(points,parInFig);
             Plot_Vec_NonParX=zeros(points,parInFig);
             Plot_Vec_NonParY=zeros(points,parInFig);

             mu_plot=zeros(1,parInFig);
%              tetaTrueNorm=zeros(1,parInFig);

    else
    
%If the number of parameters can not be divided by the number of Figures--
           if i==counter;
               
               restEnd=counter-(k*parInFig);
               rest=1:restEnd;
               k=k+1;
             pdf_plot(Plot_Vec_NonParX(:,rest),Plot_Vec_NonParY(:,rest),Plot_Vec_FitX(:,rest),Plot_Vec_FitY(:,rest),par_label_vec_plot(1:restEnd),mu_plot(1:restEnd),k)
               
           end
    end
    end

end


%% Tabnorm Function
function [TableNormalDist,Ind,graflow,grafup]=tabnorm(par_label_vec,Nor_mu,Nor_sigma,alpha,counter,repeat,Ranking,faktor)
%{
Function that creates the table TableNormalDist summarizing the estimator performance assessemnt, 
the ordering of the parameters depends on the activation of the ranking in 'Ranking={1,2}', 
graph_low and graph_up:: lower and upper limit for the Plot in order to have a better resolution of the Results in the plot
interval has to be multiply by the factor (faktor) to give a widther graphic.   

Columns in Table TableNormalDist:

TableNormalDist=['label', 'Parameter','Mean','Stand. Dev.',...
    'Low Conf. Int.', 'Up Conf. Int.'];

inputs::
par_label_vec::     Vector with the  Parameters labels
Normu::             Mean of parameter estimates
Norsigma::          Standar Deviation of parameter estimates
alfa::              Significance level for confidence interval
repeat::            Repetitions of the MonteCarlo
Rank::              Varibale for the Ranking Order
faktor::            factor to be multiply 

 %}
% Degrees of Freedom
    DegFree=repeat-1;

% Probability
    alphaup=1-alpha/2;

%t_alpha:: t Critical calculated with the matlab function tinv()
    t_alpha = tinv(alphaup,DegFree);
 
%Calculation of the lower and upper limit of the intervall of confidence
%and also buliding the new vectors for the Table.
    for i=1:counter    
        par(i)=i;
        mu(i)=Nor_mu(i);
        sigma(i)=Nor_sigma(i);
        low(i)=Nor_mu(i)-t_alpha * (Nor_sigma(i));
        graflow(i)=Nor_mu(i)- t_alpha * faktor * (Nor_sigma(i));
        up(i)=Nor_mu(i)+t_alpha * (Nor_sigma(i));
        grafup(i)=Nor_mu(i)+ t_alpha * faktor * (Nor_sigma(i));
        platz(i)=i;
%         platzB(i)=i;
    end
%% Sort
    if Ranking==1
        sigma_relat=100*(sigma./mu);
        [~,Ind] = sort(sigma_relat);
        

        label=par_label_vec(Ind);
        p=par(Ind);
        p=num2cell(p);

        sigma_relat=sigma_relat(Ind);
        sigma_relat=num2cell(sigma_relat);

        m=mu(Ind);
        m=num2cell(m);

        l=low(Ind);
        l=num2cell(l);

        u=up(Ind);
        u=num2cell(u);

        platz=num2cell(platz);

    else
        
        sigma_relat=100*(sigma./mu);
        [~,Ind] = sort(sigma_relat);

        label=par_label_vec;

        p=par;
        p=num2cell(p);

        sigma_relat=num2cell(sigma_relat);

        m=mu;
        m=num2cell(m);

        l=low;
        l=num2cell(l);

        u=up;
        u=num2cell(u);

        ay=0;
        for i=Ind;
            ay=ay+1;
            platz(i)=ay;    
        end
        platz=num2cell(platz);

    end

    %{
    TND=['label', 'Parameter','Mean','Stand. Dev.',...
    'Low Conf. Int.', 'Up Conf. Int.'];
    %}
    TableNormalDist=[['Parameter'; label],...
        ['Orig_Pos'; p'],...
        ['Mean \mu'; m';],...
        ['% StandDev \sigma'; sigma_relat';],...
        ['LCI'; l';],...
        ['UCI'; u';]];
end
%% CreatFit Function
function [NormFit,NonparFit] = createFit(par_All,points,nPar,graph_low,graph_up)


    counter=0;
    for i=1:nPar;
        counter=1+counter;
        x = par_All(i,:);
        %LegHandles = []; LegText = {};

        XGrid = linspace(graph_low(counter),graph_up(counter),points);

%% Create Normal Fit
        pd= fitdist(x', 'normal');

        %save Fit, X, Media, Standard Deviation
        NormFit.Fit{i} = pdf(pd,XGrid);
        NormFit.X{i} = XGrid;
        NormFit.Mu{i}=pd.mu;
        NormFit.Sigma{i}=pd.sigma;

        %% Create Non-Parametic Fit
        npd = fitdist(x','kernel','kernel','normal','support','unbounded');
        NonparFit.Fit{i} = pdf(npd,XGrid);
        NonparFit.X{i}=XGrid;
    end

end
%% Pdf_plot Function
function pdf_plot(valNonpar,pdf_nonpar,parNorm_vals,pdf,par_label,mu_plot,k)
% Function pdf_plot wich has not outputs but the figures and has as input: 
% Vectors with the size of variable graf for the nonparametric and normal 
% distribution, the points to be ploted as vector X, labels vector, the 
% path and experiment number to save the figures, the mean and true 
% value vector.
    
%setlimitz:: Variable that defines the maximal value of the Z-axis
    setlimitz=20;
%Setting the number of points to be ploted
    parNorm_Npts = size(parNorm_vals,1);
    NpActive=size(pdf,2);
    p_i=zeros(parNorm_Npts,NpActive);
    Active = 1:NpActive;

for i=Active
    p_i(:,i)=i;  
end

%Create the actual filled probability distribution
    C=zeros(parNorm_Npts,NpActive);
    for i=Active
        C(:,i)=i;
    end
    
    ThreeD = figure;
    
    h0=fill3(p_i,valNonpar,pdf_nonpar,C,'DisplayName','Actual Distribution');
%     legend(h0,'Non-Parametric Distribution'); 
    hold on

%plot of the line of the NonParametric Distribution
    hnorm=plot3(p_i,parNorm_vals,pdf,...;
        'LineWidth',2,...
        'Color', [0 0 0],'DisplayName','Normal Distribution');
%      legend(hnorm,'Normal Distribution'); 
    hold on
    

     [~, all2]=size(mu_plot);
    h2=plot(1:all2,mu_plot,'d','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',6,'DisplayName','True Value');
            
            legend(h2,'Mean');
            
%             legend([hnorm h2],'Normal Distribution','Mean');
    grid('on');
      set(h0, 'FaceAlpha', 0.8);
      
     
    xlabel('\theta_i','FontWeight','bold','FontSize',28,...
        'FontName','Times New Roman','FontWeight','bold');
    ylabel('Parameter Value','Rotation',-20,'FontWeight','bold','FontSize',28,...
        'FontName','Times New Roman');
    zlabel('pdf','FontWeight','bold','FontSize',28,...
        'FontName','Times New Roman','FontWeight','bold');
    ylim([-0.5 2])

    zlim([0 setlimitz]);

    set(gca,'FontSize',26,'FontName','Times New Roman')  
    set(gca,'XTick',1:NpActive) 
    set(gca,'XTickLabel',par_label,...
        'FontSize',22,'FontName','Times New Roman','FontWeight','bold') 
    set(gca, 'YScale', 'linear','FontWeight','bold');
    
%Saving the Figures 
saveas(ThreeD,[pwd '/Results/MC_3D_plot' datestr(now) '.fig']);
saveas(ThreeD,[pwd '/Results/MC_3D_plot' datestr(now) '.tif']);


end