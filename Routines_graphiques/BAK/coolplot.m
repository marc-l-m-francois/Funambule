function coolplot(xwidth, ywidth,name);
% Pour faire des figures correctes pour le papier
h = gcf; % Current figure handle
set(h,'Unit','centimeters')
set(h,'Resize','off');
set(h,'PaperPositionMode','manual');
set(h,'PaperPosition',[0 0 xwidth ywidth]);
set(h,'PaperUnits','centimeters');
set(h,'PaperSize',[xwidth ywidth]); % IEEE columnwidth = 9cm
set(h,'Position',[0 0 xwidth ywidth]);
% 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Pour ecrire du texte en Latex
% xpos, ypos must be set
%  xpos = 2;
%  ypos = 3;
% txlabel = text(xpos,ypos,'$$\mathrm{min}$$','Interpreter','latex','FontSize',9);

switch nargin
    case 3
        % Dump colored encapsulated pdf 
        % -painters pour eviter d'avoir un pdf pourri
        print('-dpdf','-loose', name,'-painters');
        % print('-depsc2','-loose', name);
end