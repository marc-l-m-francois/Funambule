function [xp1,xp2] = InverseMapApprox(XP1,XP2,X1,X2,x1,x2);
% Attention valeur approximative ! Choisir le plus proche sur 
% la grille X1 X2
% X1 X2 coo de l'image virtuelle
% XP1 XP2 coo des points d'interet
% xp1 xp2 leurs coo dans le repere x1 x2

xp1 = zeros(size(XP1));
xp2 = zeros(size(XP2));

for p=1:size(XP1,1)
    for q=1:size(XP2,2)
        DIS2 = (X1-XP1(p,q)).^2 + (X2-XP2(p,q)).^2;
        [dis,P,Q] = minmat(DIS2);
        if dis<0.707
            xp1(p,q) = x1(Q);
            xp2(p,q) = x2(P);
        else
            % Trous si le point le plus proche est trop loin
            xp1(p,q) = NaN;
            xp2(p,q) = NaN;
        end
    end
end