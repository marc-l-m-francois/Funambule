figure(1);
% Image initiale F dans son repere propre X1 X2

clf;
set(1,'Name','Image etat initial F');

% Representation de l'image de depart dans son repere pixel
image(round(F*64));colormap gray;
axis equal tight off;hold on

 % Dessin du bord du domaine d'etude
if exist('P1')
    plot(P1,P2,[cF,'-'],'LineWidth',trait_gro);
end

% Dessin du maillage
if exist('XNP1');   
    plot(XNP1 ,XNP2 ,[cF,'-'],'LineWidth',trait_gro);
    plot(XNP1',XNP2',[cF,'-'],'LineWidth',trait_gro);
end