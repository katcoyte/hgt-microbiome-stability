% Function describing underlying growth dynamics
% Note, eqns long so have broken into separate terms for clarity. For
% example, we calculate the intrinsic growth and interaction terms
% separately

function [dAdt] = plasmid_dynamics5(t, X, M, G, muS, muP, muC, muQ, delta, D, phi, merc_sus, merc_res, my_tol)

% Abundances of each population
X(X<my_tol) = 0;
S = X(1:D);         % Susceptible
P = X(D+1:2*D);     % Plasmid
C = X(2*D+1:3*D);   % Chromosome
Q = X(3*D+1:4*D);   % Both

% Interaction terms (ie the sum(a_ijXj) component of gLV model) 
compA = M*X;
compS = compA(1:D);
compP = compA(D+1:2*D);
compC = compA(2*D+1:3*D);
compQ = compA(3*D+1:4*D);

% Note, as we're modelling things like segregation rate explicitly for the
% sake of completeness we're going to split the intrinsic (net) growth rate
% into individual birth and death terms. Assuming each species has the same
% underlying death rate. Segregation is so rare none of this really impacts 
% our dynamics, but good to be thorough in case we want to model scenarios 
% where it does 
death_rate = 0.1;

% Removing death rate from "birth" term 
muS = muS + death_rate;
muP = muP + death_rate;
muC = muC + death_rate;
muQ = muQ + death_rate;

igS = S.*(muS); % intrisic growth of S
igP = P.*(muP); % intrisic growth of P
igC = C.*(muC); % intrisic growth of C
igQ = Q.*(muQ); % intrisic growth of Q

mmS = S.*(compS); % MM interactions of S 
mmP = P.*(compP); % MM interactions of P
mmC = C.*(compC); % MM interactions of C
mmQ = Q.*(compQ); % MM interactions of Q

% dxdt of S just dependent on its gLV terms
growthS = igS + mmS; 

% Split dPdt into positive and negative parts, because we only want
% segregation to occur during growth, otherwise cells start being generated
% from nowhere
growthP_pos = heaviside(igP).*igP + heaviside(mmP).*mmP;  
growthP_neg = (1-heaviside(igP)).*igP +  (1-heaviside(mmP)).*mmP;

% C behaves like S
growthC = igC + mmC; 

% Again, split growth of Q so segregation only occurs during growth
growthQ_pos = heaviside(igQ).*igQ + heaviside(mmQ).*mmQ;  
growthQ_neg = (1-heaviside(igQ)).*igQ +  (1-heaviside(mmQ)).*mmQ;

% Plasmid transfer rates
pgainS = S.*(G*P); % S gain from P
pgainS_Q = S.*(G*Q); % S gain from Q
pgainC = C.*(G*P); % C gain from P
pgainC_Q = C.*(G*Q); % C gain from Q

mutS = phi.*S; % rate that S mutates
mutP = phi.*P; % rate that P mutates

% Now build the full system of equations, growth + plasmid gain - mercury -
% death etc:

dS = growthS - pgainS - pgainS_Q + delta*growthP_pos - mutS - merc_sus.*S - death_rate*S;
dP = growthP_pos*(1-delta) + growthP_neg + pgainS + pgainS_Q - mutP - merc_res*P - death_rate*P;
dC = growthC - pgainC - pgainC_Q + delta*growthQ_pos + mutS - merc_res*C - death_rate*C;
dQ = growthQ_pos*(1-delta) + growthQ_neg + pgainC + pgainC_Q + mutP - merc_res*Q - death_rate*Q;


dAdt = [dS;dP;dC;dQ];
dAdt(X<my_tol)=0; % make sure dead things don't grow

end