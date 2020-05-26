function Range_mon_trace(NumFig,l,c,L,C)
% Pour organiser les fenetres graphiques
% NumFig Numero de la figure
% l numero de ligne
% c numero de colonne
% L nombre total des lignes
% C nombre total de colonnes
% rangeoupas a 1 pour ranger a 0 pour laisser
% efface a 1 pour effecer

% REGLAGES

bw = 10;	% Largeur de l'intervalle entre images

% INITS

NbPixEcran = get(0,'ScreenSize');   % Nombre de pixels de l''ecran
L_pix = NbPixEcran(4)-bw;
C_pix = NbPixEcran(3);

% CALCULS

DL_pix = L_pix/L;
DC_pix = C_pix/C;

x = (c-1)*DC_pix;
y = L_pix-l*DL_pix;

dx = DC_pix-bw;
dy = DL_pix - 2*bw;


% GO
figure(NumFig);
set(NumFig,'OuterPosition',[ceil(x+bw),ceil(y+bw),floor(dx),floor(dy)]);