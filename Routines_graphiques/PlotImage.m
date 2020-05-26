function PlotImage(IM,nb,varargin)
% Pour tracer F ou G et les grilles qui vont avec
% nb numero de l'image
% Varagin contient la grille du maillage


% TRACE IMAGE -------------------------------------------------------------


figure(nb);clf
imagesc(IM,[0 1]);
axis equal tight off;hold on
set(gca,'Position',[0 0 1 1]);


% TRACES COMPLEMENATIRES --------------------------------------------------


for p=1:2:nargin-2
    % Option ligne droites (mais c'est faux sur G a cause du champ de
    % deplacement pas forcement lineaire')
    % ,'r-','LineWidth',1);
    plot(varargin{p},varargin{p+1},'r+','LineWidth',1,'MarkerSize',24);
    if min(size(varargin{p},1))>1
       plot(varargin{p}',varargin{p+1}','r+','LineWidth',1,'MarkerSize',24);
    end
end

colormap(gray(256));
