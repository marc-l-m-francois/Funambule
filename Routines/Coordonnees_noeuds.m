function [xn1,xn2] = Coordonnees_noeuds(li);
% Calcule les coordonnees des noeuds du maillage EF

global XN1 XN2 DN
% XN1 XN2 coo des noeuds dans repere de F
% DN Domaine de definition du maillage (reactualise)

% Init position des noeuds dans l'image g
xn1 = XN1;
xn2 = XN2;
% Calcul des deplacements des noeuds
l = 0;
for p=1:size(XN1,1) 
    for q=1:size(XN2,2)
        if DN(p,q)
            % Les noeuds hors domaine ne bougent pas
            l = l+2;
            % Carte aux noeuds
            xn1(p,q) = XN1(p,q)+li(l-1);
            xn2(p,q) = XN2(p,q)+li(l  );
        end
    end
end