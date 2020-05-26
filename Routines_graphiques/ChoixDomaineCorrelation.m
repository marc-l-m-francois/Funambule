function [P1,P2,DF] = ChoixDomaineCorrelation(F,X1,X2);
% DF Domainde de definition de l'image F en Figure 1
% P1 P2 Coo du contour

% global F X1 X2

P1 = [];
P2 = [];
figure(1);hold on


disp('======= Choix du domaine de correlation =======');
commandwindow
rep = input('   Image complete [] ou selection d''un domaine [1] : ');
if isempty(rep)
    % Domaine de definition de l'image F = image entiere
    b = input('      Largeur de bande d''exclusion au bord (>=0) : ');
    DF = true(size(F));
    P1 = [1+b size(F,2)-b size(F,2)-b 1+b];
    P2 = [1+b 1+b size(F,1)-b size(F,1)-b];
    % Points M dans le domaine
    DF = DetectoDedans(X1,X2,P1,P2);
else
    disp('   Cliquer le contour du domaine. CLIC DROIT ou OPTION-CLIC pour les souris Mac, au dernier');
    button = 1;
    while button==1
        [PP1,PP2,button] = ginput(1);
        P1 = [P1,PP1];
        P2 = [P2,PP2];
        plot(P1,P2,'r+-','LineWidth',1,'MarkerSize',24);
    end
    % Points M dans le domaine
    DF = DetectoDedans(X1,X2,P1,P2);
end
% Completion pour representation graphique
P1 = [P1,P1(1)];
P2 = [P2,P2(1)];
plot(P1,P2,'r-','LineWidth',1);

