function D = braycd(X1,X2)
% created by Olivier Zemb
% input : X1 and X2 as two row vectors containing the observations of the
% site 1 and site 2
% output : the bray curtis distance
% Note that this function can be used in combination with pdist.
% braycd= (Si+Sj -2 Cij)/(Si+Sj)

X=[X1; X2];
Cij=sum(min(X));

Si_plus_Sj=sum(X1)+sum(X2);

D=(Si_plus_Sj-2*Cij)/Si_plus_Sj;
end