function Image_g(li);
% Pour tracer l'image g dans son repere propre x1 x2
global g P1 P2 champ cG trait_gro

% Inits graphiques
figure(2);
clf;
set(2,'Name','Image etat final g');


if ~isempty(li)
    
    % Representation de g, iumage finale dans son repere pixel
    image(round(g*64));colormap gray;
    axis equal tight off;hold on
    
    % Domaine de calcul en repere final
    if ~isempty(champ)
        %  tous les cas sauf EF - on trace le bord deforme
        [p1,p2] = feval(champ,li,P1,P2,0);
        plot(p1,p2,[cG,'-'],'LineWidth',trait_gro);
    else
        % Cas EF - on trace le maillage deforme
        [xn1,xn2] = Coordonnees_noeuds(li);
        global DN
        xn1(~DN) = NaN;
        xn2(~DN) = NaN;
        % Maillage en repere final  
        plot(xn1 ,xn2 ,[cG,'-'],'LineWidth',trait_gro);
        plot(xn1',xn2',[cG,'-'],'LineWidth',trait_gro);
    end
    
end
