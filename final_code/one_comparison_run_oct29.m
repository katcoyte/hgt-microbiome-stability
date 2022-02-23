%% Run a full set of comparisons for a given community
% Katharine Z Coyte April 2021

function [resilience, robustness, invsimp, final_diff, final_abs, final_pla_freq] = one_comparison_run_oct29(D, M, base_mu, steady_state, cost_plasmid, cost_chromosome, cost_both, G, delta, mutation_rate, sus_level, res_level, merc_level, merc_selection)

% Set intrinsic growth for each strain type (WT - resistance_cost) 
r_free = base_mu;
r_plasmid = r_free - cost_plasmid;
r_chromosome = r_free - cost_chromosome;
r_both = r_free - cost_both;
r = [r_free; r_plasmid; r_chromosome; r_both];


%% Susceptible
% First simulate scenario if all species are susceptible

% Set initial conditions, little_init tells us the proportion of community
% made up of each strain.
little_init = ones(1,D);

% Set initial conditions such that species abundances = community steady state 
init = [little_init,...
    zeros(1,D),...
    zeros(1,D),...
    zeros(1,D)].*repmat(steady_state,4,1)';

% Run one perturbation run
[~, ~, ~, A_pre, t_dur, A_dur, t_post, A_post]  = one_perturbation_run_oct29(D, M, steady_state, r, delta, G, mutation_rate, merc_level, sus_level,  res_level, init, merc_selection);

% Calculate stability, three different types
[robustness_sus, ~, ~, ~] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, 1);
[resilience_sus, ~, ~, ~] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, 3);
[invsimp_sus, final_diff_sus, absolute_sus, prop_plasmid_sus] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, 4);


%% Plasmid
% Now simulate where one species carries a resistance gene
little_init = ones(1,D);
little_init(1,1) = 0;

init = [little_init,...
    1-little_init,...
    zeros(1,D),...
    zeros(1,D)].*repmat(steady_state,4,1)';

% Run one perturbation run
[~, ~, ~, A_pre, t_dur, A_dur, t_post, A_post]  = one_perturbation_run_oct29(D, M, steady_state, r, delta, G, mutation_rate, merc_level, sus_level, res_level, init, merc_selection);

[robustness_pla, ~, ~, ~] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, 1);
[resilience_pla, ~, ~, ~] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, 3);
[invsimp_pla, final_diff_pla, absolute_pla, prop_plasmid_pla] = calculate_stability(D, A_pre, A_dur, A_post, t_post, t_dur, 4);


%% Combine results
robustness = [robustness_sus, robustness_pla];
resilience = [resilience_sus, resilience_pla];
invsimp = [invsimp_sus, invsimp_pla];
final_diff = [final_diff_sus, final_diff_pla];
final_abs = [absolute_sus, absolute_pla];
final_pla_freq = [prop_plasmid_sus, prop_plasmid_pla];


end

