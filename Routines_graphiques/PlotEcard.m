function PlotEcard(G_F,nb,ep_trait,coul,varargin)
% Pour tracer l'ecart |F-G|


% REGLAGES ----------------------------------------------------------------


% Gama reglable pour plus de lisibilite
gama = 1;
% Colormap n/b
colmap = [0:2^8-1]'*[1 1 1]/(2^8-1);
% Derniere couleur en rougepour faire apparaitre les hors limite
colmap (end,:) = [1 0 0];
% Ecart maxi (rouge au dela)
maxc = 0.1;
%maxc = max(abs(GMF(:)));

% STATS -------------------------------------------------------------------


disp('======= Statistiques sur l''ecart F-G =======');

meanFG = mean(G_F(:));
stdFG = std(G_F(:));
disp(['   Moyenne de G-F    = ',num2str(meanFG)]);
disp(['   Ecart-type de G-F = ',num2str(stdFG)]);


% ECART A MOYENNE NULLE ---------------------------------------------------


GMF = ((G_F)-meanFG*ones(size(G_F)));
% mean(GMF(:)) = 0
GMF = abs(GMF);


% TRACES ECART ------------------------------------------------------------


%while ~isempty(maxc)
    % Init
    figure(nb);clf
    set(nb,'Name','Carte d''erreur |(G-F)-<G-F>|','MenuBar','none',...
        'GraphicsSmoothing','off','ToolBar','none');
    % Trace ecard
    imagesc(GMF,[0 maxc]);
    axis equal tight off;hold on
    set(gca,'Position',[0 0 1 1]);
    colormap((colmap).^(1/gama));
    colorbar;
    % Autres traces
    for p=1:2:nargin-4
        plot(varargin{p},varargin{p+1},[coul((p+1)/2),'-'],'LineWidth',ep_trait((p+1)/2));
        if min(size(varargin{p},1))>1
           plot(varargin{p}',varargin{p+1}',[coul((p+1)/2),'-'],'LineWidth',ep_trait((p+1)/2));
        end
    end
    % Choix autre max
    %maxc= input(['    autre limite de niveau de gris (%) [>0 100] ou [] pour sortir ']);
    %maxc = min(1,maxc/100);
%end





