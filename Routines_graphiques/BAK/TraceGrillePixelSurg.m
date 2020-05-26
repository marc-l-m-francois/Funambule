figure(2);

% Raffinement PAIR
raf = 12;

% Courbe fine
xf1 = [0:1/ceil(raf*nraff*longueur):1];
ddxf1 = min(diff(xf1))/1000;
[Xcf1 ,Xcf2 ] = feval(courbe,L,xf1     );
[XcfP1,XcfP2  ] = feval(courbe,L,xf1+ddxf1);

% Calcul de ses derivees numeriques en x1
dXc1fdx1 = (XcfP1-Xcf1)/(ddxf1);
dXc2fdx1 = (XcfP2-Xcf2)/(ddxf1);
ndXcfdx1 = sqrt(dXc1fdx1.*dXc1fdx1 + dXc2fdx1.*dXc2fdx1);
esf1 = dXc1fdx1./ndXcfdx1; % vecteur tangeant
esf2 = dXc2fdx1./ndXcfdx1;
erf1 = esf2;              % vecteur normal
erf2 =-esf1;


% grille pixel tres fine
[XPf1,XPf2] = meshgrid(-0.5:1/raf:size(F,1)+0.5,-0.5:1/raf:size(F,2)+0.5);

% Trace des lignes
for p=1:size(XPf1,1)
    % Point d'avant
    xppf1 = NaN;
    xppf2 = NaN;
    for q=1:size(XPf1,2)
        if XPf2(p,q) == round(XPf2(p,q))
            % Vecteur courbe-point
            Vecf1 = XPf1(p,q)-Xcf1;
            Vecf2 = XPf2(p,q)-Xcf2;
            % Distances a la courbe
            Disf2 = Vecf1.^2 + Vecf2.^2;
            % Le mini 
            [disf2,ind] = min(Disf2);
            % Calcul des coo x1 x2 de g
            xpf2 = (Vecf1(ind)*erf1(ind)+Vecf2(ind)*erf2(ind))/R;
            xpf1 = xf1(ind);
            plot(longueur*[xppf1 xpf1],R*[xppf2 xpf2],[cF,'-'],'LineWidth',trait_fin);
            % plot([xppf1 xpf1],[xppf2 xpf2],[cF,'-'],'LineWidth',trait_fin);
            % Sauvegarde du point d'avant      
            xppf1 = xpf1;
            xppf2 = xpf2;
        end
    end
end

% Trace des colonnes
for q=1:size(XPf1,1)
    % Point d'avant
    xppf1 = NaN;
    xppf2 = NaN;
    for p=1:size(XPf1,2)
        if XPf1(p,q) == round(XPf1(p,q))
            % Vecteur courbe-point
            Vecf1 = XPf1(p,q)-Xcf1;
            Vecf2 = XPf2(p,q)-Xcf2;
            % Distances a la courbe
            Disf2 = Vecf1.^2 + Vecf2.^2;
            % Le mini 
            [disf2,ind] = min(Disf2);
            % Calcul des coo x1 x2 de g
            xpf2 = (Vecf1(ind)*erf1(ind)+Vecf2(ind)*erf2(ind))/R;
            xpf1 = xf1(ind);
            % Eviter les grandes lignes allant de x1<0 a x1>1
            if (xppf1-xpf1)^2+(xppf2-xpf2)^2<1
                plot(longueur*[xppf1 xpf1],R*[xppf2 xpf2],[cF,'-'],'LineWidth',trait_fin);
            end
            % Sauvegarde du point d'avant      
            xppf1 = xpf1;
            xppf2 = xpf2;
        end
    end
end

