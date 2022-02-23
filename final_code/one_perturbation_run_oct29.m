function [saveT, saveA, t_pre, A_pre, t_dur, A_dur, t_post, A_post]  = one_perturbation_run_oct29(D, M, steady_state, r, delta, G, mutation_rate, merc_level, sus_level, res_level, init, merc_selection)

% Note, each community is composed of 4 different populations, susceptible,
% plasmid, chromosome and both. Functionally we model each as it's own
% node in the interaction network, so we just expland the M matrix

% NOTE: in our analysis we do not allow for de novo mutation, so
% chromosomal and both populations are just zero

M = repmat(M, 4);

r_free = r(1:D);
r_plasmid = r(D+1:2*D);
r_chromosome = r(2*D+1:3*D);
r_both = r(3*D+1:4*D);

% tolerance level to prevent any populations going negative
my_tol=0;

% initialize some things
saveA = [];
saveT = [];

%% No mercury, community acclimatizing 

% Set background mercury level
if merc_selection
    merc = 0.01;
else
    merc = 0;
end

merc_sus = sus_level*merc;
merc_res = res_level*merc;

[t_pre, A_pre] = ode45(@(t,A)plasmid_dynamics5(t, A, M, G, r_free, r_plasmid, r_chromosome, r_both, delta, D, mutation_rate, merc_sus, merc_res, my_tol),[0, 500], init);

saveA = [saveA; A_pre];
saveT = [saveT;t_pre];

%% Mercury, perturbation period

merc = merc_level;

merc_sus = sus_level*merc;
merc_res = res_level*merc;

init = A_pre(end,:); 
t_init = t_pre(end);

[t_dur, A_dur] = ode45(@(t,A)plasmid_dynamics5(t, A, M, G, r_free, r_plasmid, r_chromosome, r_both, delta, D, mutation_rate, merc_sus, merc_res, my_tol),[t_init, t_init+25], init);

saveA = [saveA; A_dur];
saveT = [saveT;t_dur];

%% No mercury

if merc_selection
    merc = 0.01;
else
    merc = 0;
end

merc_sus = sus_level*merc;
merc_res = res_level*merc;

init = A_dur(end,:); 
t_init = t_dur(end);

% Use this if want to check classic stability
%init = A_pre(end,:).*rand(1,length(init)); 

[t_post, A_post] = ode45(@(t,A)plasmid_dynamics5(t, A, M, G, r_free, r_plasmid, r_chromosome, r_both, delta, D, mutation_rate, merc_sus, merc_res, my_tol),[t_init, t_init+25], init);


saveA = [saveA; A_post];
saveT = [saveT;t_post];


end
