%% Generate main simulation results
% Katharine Z Coyte April 2021

% This runs simulations to generate the main results for figures 2 & 3.
% These can then be plotted using corresponding functions.

%% Biological parameters

D = 10;                                                                     % Species number
sigma = 0.015;                                                              % Interspecies interaction strength standard deviation
self_inhibition = 0.1;                                                      % Intra-specific interaction strength

C = 0.7;                                                                    % Community connectivity, probability of interactions

cost_plasmid = 0.005;                                                       % Cost of carrying mobile resistance gene
cost_chromosome = 0.001;                                                    % Cost of carrying chromosomal resistance gene
cost_both = 0.006;                                                          % Cost of both

delta  = 0;%10^-6;                                                          % Rate of segregational loss
mutation_rate = 0;%10^-4;                                                   % Rate of de novo chromosomal mutations

gamma = 0.5*10^-3;                                                          % Average plasmid conjugation rate

res_level =0.1;                                                             % Susceptibility of resistant taxa 
sus_level = ones(D,1)*1;                                                    % Susceptibility of susceptible taxa (vector to allow for variability between taxa)
merc_level = 0.1;                                                           % Magnitude of perturbation

%% Simulation parameters

plot_on=0;                                                                  % Do you want to plot dynamics?

repeats = 100;                                                              % Number of independent communities to simulate

PmSpace =  linspace(0,1,101);                                               % Pm : Proportion of links that are positive, vary between 0 and 1
C_Space =  linspace(0,1,11);                                                % C : probability two given taxa interact, vary between 0 and 1
gamma_space = logspace(-6,-2,101);                                          % gamma : mobility, range between 10^-6 to 10^-2


g_ix = 1;                                                                   % Limitations on conjugation
                                                                            % 1. conjugation occurs between all species (with random variation) 
                                                                            % 2. only between interacting species
                                                                            % 3. between competing species
                                                                            % 4. between cooperating species


%% Initialize results matrices
% Initializing for results - sorry ugly slicing later...

all_resiliences_all = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_resiliences_subset = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_resiliences_resistant = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);

all_robustness_all = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_robustness_subset = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_robustness_resistant = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);

all_is_all = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_is_subset = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_is_resistant = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);

all_abs_all = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_abs_subset = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_abs_resistant = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);

all_pf_all = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_pf_subset = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);
all_pf_resistant = zeros(4, length(gamma_space), length(PmSpace), repeats, 2);



%% Run main simulation
parfor C_ix = 1:length(gamma_space)
    
    % Set current average gamma
    gamma = gamma_space(C_ix)
    
    % Slicing so parfor doesn't get upset
    all_resiliences_all_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_resiliences_subset_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_resiliences_resistant_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    
    all_robustness_all_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_robustness_subset_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_robustness_resistant_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    
    all_is_all_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_is_subset_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_is_resistant_slice = zeros(4, 1, length(PmSpace), repeats, 2);

    all_abs_all_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_abs_subset_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_abs_resistant_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    
    all_pf_all_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_pf_subset_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    all_pf_resistant_slice = zeros(4, 1, length(PmSpace), repeats, 2);
    
    
    for Pm_ix = 1:length(PmSpace)
        
        % Set current Pm value
        Pm = PmSpace(Pm_ix)
        
        for i = 1:repeats
            
            % set current M and Gamma matrices, note we want new matrices for
            % each repeat so they're independent
            
            G_save = gamma*ones(D)+ 0.1*gamma*randn(D);
           
            G_save(G_save<0)=0;
            [M, base_mu, steady_state] = build_M(D, C, sigma, Pm, self_inhibition);
                
            
            for m_ix = 1:2 % Now for each community simulate with and without mercury

                if m_ix == 1 % No mercury selection
                    merc_selection=0;
                elseif m_ix==2 % Prior mercury selection
                    merc_selection=1;
                end
                
                % We might want to change the G matrix such that plasmids only transfer between certain taxa so we create a saved version to not overwrite 
                [G] = alter_G_matrix(g_ix, G_save, M);
                
                % Now actually run the simulations
                [resilience, robustness, braycurtis, final_diff, final_abs, final_pla_freq] = one_comparison_run_oct29(D, M, base_mu, steady_state, cost_plasmid, cost_chromosome, cost_both, G, delta, mutation_rate, sus_level, res_level, merc_level, merc_selection);
                
                % Then save results in oir slice
                all_abs_all_slice(m_ix, 1, Pm_ix, i,:) = mean(final_abs);
                all_abs_resistant_slice(m_ix, 1, Pm_ix, i,:) = final_abs(1,:);
                all_abs_subset_slice(m_ix, 1, Pm_ix, i,:) = mean(final_abs(2:end,:));
                
                all_pf_all_slice(m_ix, 1, Pm_ix, i,:) = mean(final_pla_freq);
                all_pf_resistant_slice(m_ix, 1, Pm_ix, i,:) = final_pla_freq(1,:);
                all_pf_subset_slice(m_ix, 1, Pm_ix, i,:) = mean(final_pla_freq(2:end,:));
                
                all_resiliences_all_slice(m_ix, 1, Pm_ix, i,:) = mean(resilience);
                all_resiliences_resistant_slice(m_ix, 1, Pm_ix, i,:) = resilience(1,:);
                all_resiliences_subset_slice(m_ix, 1, Pm_ix, i,:) = mean(resilience(2:end,:));
                
                all_robustness_all_slice(m_ix, 1, Pm_ix, i,:) = mean(robustness);
                all_robustness_resistant_slice(m_ix, 1, Pm_ix, i,:) = robustness(1,:);
                all_robustness_subset_slice(m_ix, 1, Pm_ix, i,:) = mean(robustness(2:end,:));
                
                all_is_all_slice(m_ix, 1, Pm_ix, i,:) = braycurtis(1,:);
                all_is_resistant_slice(m_ix, 1, Pm_ix, i,:) = 0;
                all_is_subset_slice(m_ix, 1, Pm_ix, i,:) = braycurtis(2,:);
                
            end
            
            
        end
        
    end
    
    % Put slice values into our main results matrices 
    all_resiliences_all(:,C_ix,:,:) = all_resiliences_all_slice(:,1,:,:);
    all_resiliences_subset(:,C_ix,:,:) = all_resiliences_subset_slice(:,1,:,:);
    all_resiliences_resistant(:,C_ix,:,:) = all_resiliences_resistant_slice(:,1,:,:);
    
    all_robustness_all(:,C_ix,:,:) = all_robustness_all_slice(:,1,:,:);
    all_robustness_subset(:,C_ix,:,:) = all_robustness_subset_slice(:,1,:,:);
    all_robustness_resistant(:,C_ix,:,:) = all_robustness_resistant_slice(:,1,:,:);
    
    all_is_all(:,C_ix,:,:) = all_is_all_slice(:,1,:,:);
    all_is_subset(:,C_ix,:,:) = all_is_subset_slice(:,1,:,:);
    all_is_resistant(:,C_ix,:,:) = all_is_resistant_slice(:,1,:,:);
    
    all_abs_all(:,C_ix,:,:) = all_abs_all_slice(:,1,:,:);
    all_abs_subset(:,C_ix,:,:) = all_abs_subset_slice(:,1,:,:);
    all_abs_resistant(:,C_ix,:,:) = all_abs_resistant_slice(:,1,:,:);
    
    all_pf_all(:,C_ix,:,:) = all_pf_all_slice(:,1,:,:);
    all_pf_subset(:,C_ix,:,:) = all_pf_subset_slice(:,1,:,:);
    all_pf_resistant(:,C_ix,:,:) = all_pf_resistant_slice(:,1,:,:);
     
end



