function [params,output] = TSF_Gillespie(OutputFolder,OutputName,save_output,plot_dynamics,plot_histogram,plot_analysis)
% TSF_GILLESPIE Simulates transcriptional bursting from the TSF model using the Gillespie method. 
% The simulation calculates the dynamics of the state of the system (Closed or Open) and the occupancy (number of bounded transcriptional binding sites by transcriptional factors).  
%
%   Inputs:
%   OutputFolder    to which folder the data and figures will be saved 
%   OutputName      under what name the data and figures will be saved
%   save_output     if 1 (or 0) the data and figures will (or will not) be saved  
%   plot_dynamics   if 1 (or 0) the bursting dynamics will (or will not) be ploted 
%                   (NOTE: to plot the original occupancy, before smoothing, see section "plots of individual iteration" below) 
%   plot_histogram  if 1 (or 0) the probabilities histogram will (or will not) be ploted
%   plot_analysis   if 1 (or 0) the meta-analysis will (or will not) be ploted
%   
%   Outputs:
%   params          data structure containing the parameters and initial conditions of the simulation
%   output          data structure containing the output results from the simulation
    
    
    
    %% settings
    
    % simulation time parameters:
    Total_time = 500*60;        % [min]
    k_off_1 = 1/4;              % [1/min]   (Time is normalized by koff_active. See also kinetic parameters below)
    SmoothingTimeInterval = 10; % [min]     (This is the moving average time interval used in "movmean_time" function)

    % model parameters:
    nu = 1;         % transition M_1->0 dependency on occupancy (M_1->0_n)
    N = 6;          % number of binding sites    
    alpha = 0.5;    % alpha = k_on/k_off_1 = [A]/K_d_1  
    kappa = 1;      % kappa = k_off_0/k_off_1          
    mu = 1/5;       % mu    = M_0->1/M_1->0_0           
    sigma = 0.5;    % sigma = M_1->0,0/k_off_1         
    
    
    % test parameter: 
    parameter2test = 'alpha';  % which parameter is tested (See also switch below);
    
    switch parameter2test
        case 'N'
            factor = 1;
            N = factor*[2 4 6 8 10 12 14 16 18 20];
            Test_parameter = N;                   
            num_of_iter = length(N);
        
        case 'alpha' 
            alpha = 0.02:0.06:2;
%             % Manuscript Figure 4 (Uncomment to plot; See also section "plot analysis" below)  
%             factor = 1;
% %             factor = 0.5; % Manuscript Figure 4C-C' (Uncomment to plot)
%             x_microns = 5:5:35;
%             x_factor = 30;
%             alpha = factor*(1./(1+(x_microns./x_factor).^2));

            Test_parameter = alpha;
            num_of_iter = length(alpha);

        case 'kappa'
            factor = 1;
            kappa = factor*[10 100];

            Test_parameter = kappa;
            num_of_iter = length(kappa);
        
        case 'mu'
            factor = 1;
            mu = factor*[0.1 1];

            Test_parameter = mu;
            num_of_iter = length(mu);
        
        case 'sigma'
            factor = 1;
            sigma = factor*[0.1 1];

            Test_parameter = sigma;
            num_of_iter = length(sigma);
    end
    
    % rearrange model parameters (as vectors 1xnum_of_iter)
    N       = N.*ones(1,num_of_iter);            
    alpha   = alpha.*ones(1,num_of_iter);     
    kappa   = kappa.*ones(1,num_of_iter);   
    mu      = mu.*ones(1,num_of_iter);        
    sigma   = sigma.*ones(1,num_of_iter);  
    
    
    % kinetic parameters (as a function of the model parameters):
    k_on        = alpha.*k_off_1;   % association rate of TFs to a single binding site
    k_off_0     = kappa.*k_off_1;   % disassociation rate of TFs from a single binding site, at the closed system state
    M_1to0_0    = sigma.*k_off_1;   % transition rate of the system from open to closed state, with zero occupancy 
    M_0to1      = mu   .*M_1to0_0;  % transition rate of the system from closed to open state
    
    % transcription threshold:
    BurstOccupancyThreshold = 1.*ones(1,num_of_iter); % number of occupied binding sites considered as transcriptional_ON  



    %% initial conditions
    m0 = 0; % initial state (0 closed , 1 open)
    n0 = 0; % initial occupancy (in [0,N]) 
    t0 = 0; % initial time 


    
    %% other settings
    
    % variables:
    OFF_time = [];           
    OFF_time_SEM = [];  
    ON_time = [];        
    ON_time_SEM = [];
%     % uncomment for amplitues (See section "Calculate outputs for analysis" below):
%     BurstAmplitude = [];
%     BurstAmplitude_SEM = [];
%     BurstAmplitude_scatter_plot_x = [];
%     BurstAmplitude_scatter_plot_y = [];
    
    % plots:
    font_size_titles = 30;
    font_size_ticks = 30;
    

    
    %% loop over the values of test parameter
    for i = 1:num_of_iter
    
            
        %% gillespie
        
        n = n0;         % occupancy (in [0,N])
        occupancy = n;  % occupancy series vector - reccords the occupancy over time
        m = m0;         % state (closed 0, or open 1)
        State = m;      % state series vector - reccords the state over time
        t = t0;         % time series vector - records the time points
        
        while t(end) < Total_time % the simulation runs as long as the time parameter is lower than the total time
            
            % kinetic parameters:
            K_on_n      = (N(i)-n)*k_on(i);                 % binding rate of a single transcriptional factor to multiple binding sites
            K_off_n_0   = n*k_off_0(i);                     % unbinding rate of a single transcriptional factor from multiple binding sites, in the closed system state  
            K_off_n_1   = n*k_off_1;                        % nbinding rate of a single transcriptional factor from multiple binding sites, in the open system state
            M_0to1_     = M_0to1(i);                        % Note: M_0->1 is constant
            M_1to0_n_   = M_1to0_0(i)/(kappa(i)^(nu*n));    % transition rate from open to closed state depends on occupancy
            
            % calculate probabilities:
            rate_bind = K_on_n;
            if m == 0
                rate_unbind = K_off_n_0;
                rate_transition = M_0to1_;
            elseif m == 1   
                rate_unbind = K_off_n_1;
                rate_transition = M_1to0_n_; 
            end
            rate_total = rate_bind + rate_unbind + rate_transition;    % rate of an event
                
            P_bind = rate_bind/rate_total;      % probability of a single TF binding event
            P_unbind = rate_unbind/rate_total;  % probability of a single TF unbinding event
            P_uniform = rand;                   % a random number from a uniform distribution in (0,1)
            
            % choose which event to execute (set a new state or a new occupancy):
            if P_uniform <= P_bind 
                % TF binding
                n = n + 1; 
            elseif P_uniform > P_bind && P_uniform <= (P_bind + P_unbind)
                % TF unbinding
                n = n-1; 
            elseif P_uniform > (P_bind + P_unbind)
                % transition between states
                m = double(~m); 
            end
            
            % calculate event's time point:
            t = [t, t(end) + exprnd(1/rate_total)];

            % set the state and ocuupancy at the event's time point:
            State = [State, m];
            occupancy = [occupancy, n];

        end

        % moving average:
        Smoothed_occupancy = movmean_time(occupancy,t,SmoothingTimeInterval); 
        


        %% plots of individual iteration 
        
        % occupany dynamics:
        if plot_dynamics
%             % Uncomment to plot the original (before smoothing) occupancy plot  
%             figure;
%             yyaxis left
%             plot(t,occupancy,'-b','LineWidth',1)
%             title('Original occupancy plot','FontSize',font_size_titles)
%             xlabel('time [minutes]','FontSize',font_size_titles);
%             xlim([0 Total_time])
% %             xlim([300 900])
%             ylabel('Occupancy','FontSize',font_size_titles)
%             ylim([0 N(i)])
%             yticks(0:1:N(i))
%             hold on
%             yyaxis right
%             plot(t, State,'Color', [255/255 132/255 0],'LineWidth',0.3)
%             ylabel('State','FontSize',font_size_titles)
%             ylim([0 N(i)])
%             yticks([0 1])
%             ax = gca; 
%             ax.FontSize = font_size_ticks;

            % Smoothed occupancy plot using movemean_time function
            F_dynamics(i).fig = figure;
            yyaxis left
            plot(t,Smoothed_occupancy,'-b','LineWidth',1) 
            title(['Smoothed occupancy. Window__dur = ' num2str(SmoothingTimeInterval)],'FontSize',font_size_titles)
            xlabel('time [minutes]','FontSize',font_size_titles);
            xlim([0 Total_time])
%             xlim([300 900])
            ylabel('Smoothed occupancy','FontSize',font_size_titles)
            ylim([0 N(i)])
            yticks(0:1:N(i))
            hold on
            yyaxis right
            plot(t, State,'Color', [255/255 132/255 0],'LineWidth',1)
            ylabel('State','FontSize',font_size_titles)
            ylim([0 N(i)])
            yticks([0 1])
            ax = gca; 
            ax.FontSize = font_size_ticks;
        end
        

        % probabilities histogram:
        if plot_histogram
            dur = diff(t);  % discrete duartaions between events
            
    %         occ_prob = zeros(1,N(i) + 1);       % allocating space
            occ_prob_off = zeros(1,N(i) + 1);   % allocating space
            occ_prob_on = zeros(1,N(i) + 1);    % allocating space
            
            for itr = 0:N(i)
    %             occ_prob(itr + 1) = sum(dur(occupancy(1:end-1) == itr))/t(end);   % caculating probability of occurance by temporal duration of a state
                dur_off_for_occ_prob_off = dur.*(~State(1:end-1));
                occ_prob_off(1,itr + 1) = sum(dur_off_for_occ_prob_off(occupancy(1:end-1) == itr))/t(end);
                dur_on_for_occ_prob_on = dur.*State(1:end-1);
                occ_prob_on(1,itr + 1) = sum(dur_on_for_occ_prob_on(occupancy(1:end-1) == itr))/t(end);
            end
        
%             % Uncomment to plot a combined (without separating open and closed states) probabilities histogram
%             figure;
%             bar(0:N(i),occ_prob)
%             title('Occurence of occupancy')
%             xlabel('Occupancy')
%             ylabel('Occurence')
%             ylim([0 1])
%             yticks(0:0.2:1)

        
            % Separated (closed and open states) probabilities histogram plot 
            F_histogram(i).fig = figure;
            h = bar(0:N(i),[occ_prob_off ; occ_prob_on]');
            title(['N = ' num2str(N(i)) '\alpha = ' num2str(alpha(i)) '  \kappa = ' num2str(kappa(i)) '  \mu = ' num2str(mu(i)) '  \sigma = ' num2str(sigma(i))])
            xlabel('Occupancy','fontsize',font_size_titles)
            ylabel('Occurence','fontsize',font_size_titles)
            ylim([0 1])
            yticks(0:0.2:1)
            set(h, {'DisplayName'}, {' Closed ',' Open'}')
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',font_size_ticks)
            lgd = legend();
            lgd.FontSize = font_size_titles;
            lgd.Orientation = 'horizontal';
            lgd.Location = 'northeast';
            dim = [0.6 0.48 0.3 0.3]; % [x_coordinate y_coordinate length height]
            str = {['S_0 = ' num2str(round(sum(occ_prob_off),1))],['S_1 = ' num2str(round(sum(occ_prob_on),1))]};
            annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',font_size_titles,'EdgeColor',[1 1 1]);
        end

        
        %% Calculate outputs for analysis (of the dynamics) plot 

        TranscriptionON_timeInd = Smoothed_occupancy >= BurstOccupancyThreshold(i);     % transcription-ON time indexes  
        TranscriptionTurn_timeInd = logical([abs(diff([TranscriptionON_timeInd, 2]))]); % transcription initiation/termination (OFF->ON or ON->OFF) time indexes
        TranscriptionTurn_time = t(TranscriptionTurn_timeInd);                          % transcription initiation/termination times
        TranscriptionState_dur = diff([0 TranscriptionTurn_time]);                      % transcription-ON or transcription-OFF durations (this vector does not recognize which period is ON or OFF, but only consecutive peeriods of ON-OFF-ON-OFF-... ) 
        
        TranscriptionState_ind = 1:2:length(TranscriptionState_dur);
        
        if Smoothed_occupancy(1) < BurstOccupancyThreshold(i)
            if TranscriptionState_ind(end) == length(TranscriptionState_dur)
                T_off = TranscriptionState_dur(TranscriptionState_ind);     % pauses durations
                TranscriptionState_ind(1) = [];
                T_on = TranscriptionState_dur(TranscriptionState_ind - 1);  % bursts durations
            else
                T_off = TranscriptionState_dur(TranscriptionState_ind);     % pauses durations
                T_on = TranscriptionState_dur(TranscriptionState_ind + 1);  % bursts durations
            end
        elseif Smoothed_occupancy(1) >= BurstOccupancyThreshold(i)
            if TranscriptionState_ind(end) == length(TranscriptionState_dur)
                T_on = TranscriptionState_dur(TranscriptionState_ind);      % bursts durations
                TranscriptionState_ind(1) = [];                             
                T_off = TranscriptionState_dur(TranscriptionState_ind - 1); % pauses durations
            else
                T_on = TranscriptionState_dur(TranscriptionState_ind);      % bursts durations
                T_off = TranscriptionState_dur(TranscriptionState_ind + 1); % pauses durations
            end
        else
            error(['state ' num2str(m0) ' does not exsit'])
        end
        
        %%% outputs 
        ON_time         = [ON_time, mean(T_on)];                            % average burst duration
        ON_time_SEM     = [ON_time_SEM, std(T_on,1)/sqrt(length(T_on))];    % standard error of the mean
        OFF_time        = [OFF_time, mean(T_off)];                          % average pause
        OFF_time_SEM    = [OFF_time_SEM, std(T_off,1)/sqrt(length(T_off))]; % standard error of the mean
        

%         % uncomment for amplitues:
%         BurstAmplitude = [BurstAmplitude, mean(Smoothed_occupancy(Smoothed_occupancy >= BurstOccupancyThreshold(i)))];    
%         BurstAmplitude_STD = [BurstAmplitude_SEM, std(Smoothed_occupancy(Smoothed_occupancy >= BurstOccupancyThreshold(i)))/length(Smoothed_occupancy(Smoothed_occupancy >= BurstOccupancyThreshold(i)))];
%         BurstAmplitude_scatter_vector = [];
%         t_ind = 1;
%         while t_ind <= length(Transcription_loc)
%             Transcription_state = Transcription_loc(t_ind); % 0 (or 1) transcription is off (or on) 
%             if Transcription_state == 1
%                 n_summation = 0;
%                 c = 0;
%                 while Transcription_state == 1
%                     n_t_ind = Smoothed_occupancy(t_ind);
%                     n_summation = n_summation + n_t_ind;
%                     c = c + 1;
%                     t_ind = t_ind + 1;
%                     if t_ind <= length(Transcription_loc)
%                         Transcription_state = Transcription_loc(t_ind);
%                     else
%                         break
%                     end
%                 end
%                 n_mean = n_summation/c;
%                 BurstAmplitude_scatter_vector = [BurstAmplitude_scatter_vector , n_mean];
%             else
%                 t_ind = t_ind + 1;
%             end
%         end
%         BurstAmplitude_scatter(i).scatter_vector = BurstAmplitude_scatter_vector; 
%         BurstAmplitude_scatter(i).x_value = (factor/Test_parameter(i))*ones(1,length(BurstAmplitude_scatter_vector));
%         BurstAmplitude_scatter_plot_x = [BurstAmplitude_scatter_plot_x BurstAmplitude_scatter(i).x_value];
%         BurstAmplitude_scatter_plot_y = [BurstAmplitude_scatter_plot_y BurstAmplitude_scatter(i).scatter_vector];
        
    end
    


    %% plot analysis
    
    if plot_analysis
        
        % Plot ON-time dependence on test parameter
        F_ON_time = figure('units','normalized','outerposition',[0 0 1 1]);
        errorbar(Test_parameter,ON_time,ON_time_SEM,'.','MarkerSize',60,'LineWidth',5,'CapSize',15)
        xlim([0 1.1*max(Test_parameter)])
        ylim([0 1.1*max(ON_time)+max(ON_time_SEM)])
%         title({['ON-time (n>=' num2str(BurstOccupancyThreshold(1)) ')'] ; 'Mean (SEM)'},'fontsize',font_size_titles);
%         xlabel(parameter2test,'fontsize',font_size_titles);
%         ylabel('time [min]','fontsize',font_size_titles);
        ax = gca; 
        ax.FontSize = font_size_ticks;
%         % Manuscript Figure 4 (Uncomment to plot)
%         errorbar(((factor./Test_parameter-1).^0.5).*x_factor,ON_time,ON_time_SEM,'.','MarkerSize',60,'LineWidth',5,'CapSize',15) % Figure 4 (Uncomment to plot)  
%         xlim([0 60])
%         xticks(5:5:60)
%         ylim([0 100])
%         yticks(0:10:100)
%         xlabel('\mum from distal end','fontsize',font_size_titles);


        % Plot OFF-time dependence on test parameter
        F_OFF_time = figure('units','normalized','outerposition',[0 0 1 1]);
        errorbar(Test_parameter,OFF_time,OFF_time_SEM,'.','MarkerSize',60,'LineWidth',5,'CapSize',15) 
        xlim([0 1.1*max(Test_parameter)])
        ylim([0 1.1*max(OFF_time)+max(OFF_time_SEM)])
%         title({['off period (n<' num2str(BurstOccupancyThreshold(1)) ')'] ; 'Mean (SEM)'},'fontsize',font_size_titles);
%         xlabel([num2str(factor) '/' parameter2test],'fontsize',font_size_titles);
%         ylabel('Time (\tau)','fontsize',font_size_titles);
        ax = gca; 
        ax.FontSize = font_size_ticks; 
%         % Manuscript Figure 4 (Uncomment to plot)
%         errorbar(((factor./Test_parameter-1).^0.5).*x_factor,OFF_time,OFF_time_SEM,'.','MarkerSize',60,'LineWidth',5,'CapSize',15)  
%         xlim([0 60])
%         xticks(5:5:60)
%         ylim([0 100])
%         yticks(0:10:100)
%         xlabel('\mum from distal end','fontsize',font_size_titles);
        

%         % uncomment to plot bursts amplitudes:
%         F_amplitudes = figure('units','normalized','outerposition',[0 0 1 1]);
% %         errorbar(factor./Test_parameter,BurstAmplitude,BurstAmplitude_SEM,'.','MarkerSize',20) 
%         plot(BurstAmplitude_scatter_plot_x,BurstAmplitude_scatter_plot_y,'.')
%         hold on
%         plot(factor./Test_parameter,BurstAmplitude,'-g','LineWidth',3)
%         xlim([0 10])
%         ylim([0 max(N)])
%         ax = gca; 
%         ax.FontSize = font_size_ticks;
%         title({'burst amplitude' ; ['Mean (SEM) from n>=' num2str(BurstOccupancyThreshold)]},'fontsize',font_size_titles); % use STD 
%         xlabel([num2str(factor) '/' parameter2test],'fontsize',font_size_titles);
%         ylabel(['n (from n>=' num2str(BurstOccupancyThreshold) ')'],'fontsize',font_size_titles);



%         % uncomment to display parameters:
%         fprintf(['Parameters:\n nu = ' num2str(nu) '\n N = ' num2str(N(1))'\n alpha = ' num2str(alpha) '\n kappa = ' num2str(kappa(1)) '\n mu = ' num2str(mu(1)) '\n sigma = ' num2str(sigma(1)) '\n kappa_o_f_f__1 = ' num2str(k_off_1(1)) ]);

    end
    

    
%% Save
    
    % arranging parameters and initial conditions:
    params.Total_time = Total_time;
    params.koff_1 = k_off_1;
    params.window_dur = SmoothingTimeInterval; 

    params.nu = nu;
    params.N = N;
    params.alpha = alpha;
    params.kappa = kappa;
    params.mu = mu;
    params.sigma = sigma;

    params.Test_parameter = Test_parameter;

    params.kon = k_on;           
    params.koff_0 = k_off_0; 
    params.M_1to0_0 = M_1to0_0;     
    params.M_0to1 = M_0to1;           
    
    params.BurstOccupancyThreshold = BurstOccupancyThreshold;

    params.m0 = m0;
    params.n0 = n0;         
    params.tau0 = t0;       
    

    % arranging simulation outputs:
    output.BurstDuration = ON_time;
    output.BurstDuration_SEM = ON_time_SEM;
    output.OffPeriod = OFF_time;
    output.OffPeriod_SEM = OFF_time_SEM;
%     % uncomment to save bursts amplitudes data:
%     output.BurstAmplitude = BurstAmplitude;
%     output.BurstAmplitude_SEM = BurstAmplitude_SEM;
%     output.BurstAmplitude_scatter = BurstAmplitude_scatter;

    % saving:
    if save_output
        % create folder, if output folder does not exist
        if ~exist(OutputFolder, 'dir')
           mkdir(OutputFolder)
        end
    
        % save data
        SaveData = [OutputFolder '/' OutputName ' _Data.mat']; % for mac iOS
        save(SaveData,'params','output','-v7.3')
        
        % save figures and movie
        if plot_dynamics
            for i = 1:num_of_iter
                SaveFigure_F_dynamics = [OutputFolder '/' OutputName ' _Dynamics' num2str(i) '.fig']; % for mac iOS
                savefig(F_dynamics(i).fig,SaveFigure_F_dynamics)
            end
        end
        if plot_histogram
            for i = 1:num_of_iter
                SaveFigure_F_histogram= [OutputFolder '/' OutputName ' _Histogram' num2str(i) '.fig']; % for mac iOS
                savefig(F_histogram(i).fig,SaveFigure_F_histogram)
            end
        end
        if plot_analysis
            SaveFigure_F_ON = [OutputFolder '/' OutputName ' _ON_TIME.fig']; % for mac iOS
            savefig(F_ON_time,SaveFigure_F_ON)
            SaveFigure_F_OFF = [OutputFolder '/' OutputName ' _OFF_TIME.fig']; % for mac iOS
            savefig(F_OFF_time,SaveFigure_F_OFF)
%             % uncomment to save bursts amplitudes plot:
%             SaveFigure_F_amplitudes = [OutputFolder '/' OutputName ' _amplitudes.fig']; % for mac iOS
%             savefig(F_amplitudes,SaveFigure_F_amplitudes)
        end
    end

end