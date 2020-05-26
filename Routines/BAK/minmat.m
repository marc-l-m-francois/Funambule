function [m,p,q] = minmat(M);
% m valeur min de la matrice
% p numero de ligne du min
% q numero de colonne du min

[mm,P] = min(M);
[m ,q]  = min(mm);
p = P(q);