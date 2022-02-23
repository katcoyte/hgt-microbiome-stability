% Function to alter the G matrix such that plasmids only transfer between
% certain taxa

% G_save : Baseline conjugation matrix
% M : Interaction network
% g_ix conjugation type, where:
    % g_ix = 1 : conjugation occurs between all species (with random variation)
    % g_ix = 2 : only between interacting species
    % g_ix = 3 : only between competing species
    % g_ix = 4 : only between cooperating species

function [G] = alter_G_matrix(g_ix, G_save, M)

D = size(M,1);

if g_ix == 1 % All
    G = G_save;
    
elseif g_ix==2 % Interacting
    G = G_save;
    foo = M;
    foo(abs(M)>0)=1;
    foo = foo+foo';
    foo(foo>0)=1;
    G(foo==0)=0;
    
elseif g_ix==3 % Competing
    G = G_save;
    foo = M;
    foo(M>0) = 0;
    foo(M<0)=1;
    foo = foo+foo';
    foo(foo>0)=1;
    foo(eye(D)==1)=1;
    G(foo==0)=0;
    
elseif g_ix==4 % Cooperating
    G = G_save;
    foo = M;
    foo(M<0) = 0;
    foo(M>0)=1;
    foo = foo+foo';
    foo(foo>0)=1;
    foo(eye(D)==1)=1;
    G(foo==0)=0;
end

end