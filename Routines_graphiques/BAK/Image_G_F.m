function Image_G_F(G,DG);

global F

figure(3);
% Carte d'erreur abs(G-F) dans son repere propre X1 X2
gamma_video = 4;
myfontsize = 18;
prof = 12;
clf;
set(3,'Name','Carte d''erreur abs(G-F)');

% Colormap non lineaire pour G-G
nivgris = ([0:(2^prof-1)]/(2^prof-1)).^(1/gamma_video);
gray2 = [nivgris',nivgris',nivgris'];
colormap(gray2);

% Representation de l'ecart
G_F = G-F;
% On annule l'effet des points hors domaine
G_F(~DG) = 0;

imagesc(abs(G_F));hold on
set(gca,'FontSize',myfontsize,'CLim',[0 1]);colorbar;axis equal tight off;

