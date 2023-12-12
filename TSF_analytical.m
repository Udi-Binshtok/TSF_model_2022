function [Pn_full,Pn,Pm] = TSF_analytical(OutputFolder,OutputName,save_output,params)
%TSF_ANALYTICAL Solves and plots the analytical steady state solution of the TSF model equations.  
%   
%   Inputs:
%   OutputFolder    to which folder the data and figures will be saved 
%   OutputName      under what name the data and figures will be saves
%   save_output     if 1 (or 0) the data and figures will (or will not) be saved
%   params          providing parameters (otherwise the function will use its default parameters - in the code)
%
%   Outputs:
%   pn_full         full calculated probabilities for occupancy (n). pn_full(1:N+1) = probabilities for system in inactive state (m=0), pn_full(N+1:2*(N+1)) = probabilities for system in active state (m=1)   
%   pn              the combined probabilies (pn from inactive + active states)
%   pm              the probabilities of the system to be in an inactive (pm(1)) or in an active state (pm(2))  
    
    
    %%% parameters
    if nargin < 4 
        params = get_parameters;
    end
    
    %%% Solve differential equations at steady state
    
    N = params.N; % number of sites
    
    x0Lin = (1/(2*(N+1))).*ones(1,2*(N+1)); % uniform initial pn
    LinConst = @(x)fminconstrLin(x,params);
    x = fmincon(@(x)0,x0Lin,[],[],ones(1,2*(N+1)),1,zeros(1,2*(N+1)),ones(1,2*(N+1)),LinConst);
    
    Pn_full = x;
    Pn = Pn_full(1:(N+1)) + Pn_full((N+2):2*(N+1));
    Pm = [sum(Pn_full(1:(N+1))),sum(Pn_full((N+2):2*(N+1)))];
    
    %%% plot
    F = plot_bar(params,Pn_full);


    %%% saving
    if save_output
        %%% create folder, if output folder does not exist
        if ~exist(OutputFolder, 'dir')
           mkdir(OutputFolder)
        end
    
        %%% save data
        output.pn_full = Pn_full;
        output.pn = Pn;
        output.pm = Pm;
        SaveData = [OutputFolder '/' OutputName ' _Analytical _Data.mat']; % for mac iOS
        save(SaveData,'params','output','-v7.3')
        
        %%% save figures and movie
        SaveFigure_F = [OutputFolder '/' OutputName ' _analytical.fig']; % for mac iOS
        savefig(F,SaveFigure_F)
        
    end

end



%% parametrs

function params = get_parameters

    params.nu = 1;
    params.N = 4;           % total number of binding sites.    
    params.alpha = 0.5;     % alpha = k_on/k_off_1 = [A]/K_d_1              
    params.kappa = 1000;    % kappa = k_off_0/k_off_1            
    params.mu =  1/5;       % mu    = M_0->1/M_1->0_0           
    params.sigma = 0.5;     % sigma = M_1->0_0/k_off_1           


end


%% solver

function [c,ceq] = fminconstrLin(x,params)
    
    nu = params.nu;
    N = params.N;
    alpha = params.alpha;
    mu = params.mu;
    kappa = params.kappa;
    sigma = params.sigma;
    
    c = []; % No nonlinear inequality
    ceq = paramfunLin(x,nu,N,alpha,mu,kappa,sigma); % fsolve objective is fmincon nonlinear equality constraints
end


function F = paramfunLin(x,nu,N,alpha,mu,kappa,sigma)
    
    % variables are:  x(1) = p_0_0 , x(2) = p_1_0 , ... , x(N+1) = p_N_0 , x(N+2) = p_0_1 , x(N+3) = p_1_1 , ... , x(2(N+1)) = p_N_1 
    
    %%% equations for Closed state
    F(1) = -x(1)*(N*alpha + mu*sigma) + x(2)*kappa + x(N+2)*sigma; 
    for i = 2:1:N
        F(i) = -x(i)*((N+1-i)*alpha + (i-1)*kappa + mu*sigma) + x(i-1)*(N+2-i)*alpha + x(i+1)*i*kappa + x(N+1+i)*sigma/(kappa^(nu*(i-1)));
    end
    F(N+1) = -x(N+1)*(N*kappa + mu*sigma) + x(N)*alpha + x(2*(N+1))*sigma/(kappa^(nu*N));
    
    %%% equations for Open state
    F(N+2) = -x(N+2)*(N*alpha + sigma) + x(N+3) + x(1)*mu*sigma;
    for i = (N+3):1:(2*(N+1)-1)
        F(i) = -x(i)*((2*(N+1)-i)*alpha + (i-(N+2)) + sigma/(kappa^(nu*(i-(N+1)-1)))) + x(i-1)*(2*(N+1)-i+1)*alpha + x(i+1)*(i-(N+1)) + x(i-(N+1))*mu*sigma;
    end
    F(2*(N+1)) = -x(2*(N+1))*(N + sigma/(kappa^(nu*N))) + x(2*(N+1)-1)*alpha + x(N+1)*mu*sigma;

end



%% plot

function F = plot_bar(params,pn_full)
    close all
    F = figure;
    N = params.N;
    h = bar(0:N,[pn_full(1:(N+1))' , pn_full((N+2):2*(N+1))']);
    title(['\alpha = ' num2str(params.alpha) '  \kappa = ' num2str(params.kappa) '  \mu = ' num2str(params.mu) '  \sigma = ' num2str(params.sigma)])
%     xlabel('Occupancy (n)')
%     ylabel('Propability at steady state')
    ylim([0 1])
    yticks(0:0.2:1)
%     xlim([0 N])
    xticks(0:1:N)
    set(h, {'DisplayName'}, {' p_n_,_0 ',' p_n_,_1'}')
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',30)
    lgd = legend();
    lgd.FontSize = 30;
    lgd.Orientation = 'horizontal';
    lgd.Location = 'northeast';
    dim = [0.6 0.45 0.3 0.3]; % [x_coordinate y_coordinate length height]
    str = {['S_0 = ' num2str(round(sum(pn_full(1:(N+1))),1))],['S_1 = ' num2str(round(sum(pn_full((N+2):2*(N+1))),1))]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',30,'EdgeColor',[1 1 1]);
end