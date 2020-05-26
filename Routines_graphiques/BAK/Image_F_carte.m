figure(1);
% Image initiale F dans son repere propre X1 X2

clf;
set(1,'Name',['Carte de ',rep,' sur image etat initial F']);

% Dessin carte
h=imagesc(carte);
colormap hot 
masque = 1-F;
masque(~DG) = F(~DG);
set(h, 'AlphaData',masque);
hold on;axis equal tight off;colorbar

% Dessin du bord du domaine d'etude
if exist('P1')
    plot(P1,P2,[cF,'-'],'LineWidth',trait_gro);
end

% Dessin du maillage
if exist('XNP1');   
    plot(XNP1 ,XNP2 ,[cF,'-'],'LineWidth',trait_gro);
    plot(XNP1',XNP2',[cF,'-'],'LineWidth',trait_gro);
end

% Limitation des valeurs affichees
repp = caxis;
while ~isempty(repp) && size(repp,2)==2
    caxis(repp);
    set(gca,'FontSize',myfontsize);
    commandwindow
    repp = input('        Entrer des valeurs limites [bas haut] ou OK [] : ');
end
commandwindow