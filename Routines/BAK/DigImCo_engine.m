function [li,psi] = DigImCo_engine(li,cptmax,precision);


global F DF g G_F X1 X2 champ P1 P2


% REGLAGES ----------------------------------------------------------------


% Type d'interpollation retenu
type_interp	= 'linear'; 


% INITS GRAPHIQUES --------------------------------------------------------


global cG trait_gro myfontsize


% INITS -------------------------------------------------------------------


% Gradients de l'image F
[gradF10,gradF20]  = gradient(F);


% INIT NEWTON  ------------------------------------------------------------


% Init compteur de boucles et ecart
cpt = 1;
% Champ actuel
[x1,x2] = feval(champ,li,X1,X2,0);
% Image G = image g interpolee
G = interp2(g,x1,x2,type_interp,NaN);
% Son domaine de definition intersection avec celui de F
DG = isfinite(G) & DF;
% Surface = nombre de pixels correles
S = sum(sum(DG));
% Ecart GG-F
G_F = G-F;
% On annule l'effet des points hors domaine
G_F(~DG) = 0;
% Erreur actuelle
Psi(cpt) = sum(sum(G_F.^2))/S;
% Variables internes au Newton
% gradF_dxdl = zeros(size(F,1),size(F,2),size(li,1));
gradF_dxdlp = zeros(size(F,1),size(F,2));
gradF_dxdlq = zeros(size(F,1),size(F,2));
M = zeros(max(size(li)));
V = zeros(max(size(li)),1);



% NEWTON ------------------------------------------------------------------


% Init correcteur pour rentrer dans la boucle
dli = inf*ones(size(li));

while max(abs(dli))>precision && cpt<cptmax
    
    disp(['   iteration ',num2str(cpt)]);
    % Elimination des points hors domaine
    gradF1 = gradF10;
    gradF2 = gradF20;
    gradF1(~DG) = 0;
    gradF2(~DG) = 0;
    % Calcul matrice et vecteur
    for p=1:max(size(li))
         % Calcul de la p-ieme derivee en li
        [dx1dL,dx2dL] = feval(champ,li,X1,X2,p);
        gradF_dxdlp  =  gradF1.*dx1dL +  gradF2.*dx2dL;
        % Vecteur
        V(p) = -sum(sum(gradF_dxdlp.*G_F));
        % Matrice
        for q=1:p
            if q~=p
                [dx1dL,dx2dL] = feval(champ,li,X1,X2,q);
                gradF_dxdlq  =  gradF1.*dx1dL +  gradF2.*dx2dL;
                M(p,q) = sum(sum( gradF_dxdlp.*gradF_dxdlq ));
                M(q,p) = M(p,q);
            else
                M(p,q) = sum(sum( gradF_dxdlp.*gradF_dxdlp));
            end
        end
    end
    
    % Verification de la determination de la matrice
    detM = abs(det(M));
    if detM<1E-6
        warning(['      La matrice a determinant ',num2str(detM),' presque nul)']);
    end       
    
    % Correcteur
    dli = (M\V)';
    li = li + dli;
    
    % Champ actuel
    [x1,x2] = feval(champ,li,X1,X2,0);
    % Image G = image g interpolee
    G = interp2(g,x1,x2,type_interp,NaN);
    % Son domaine de definition intersection avec celui de F
    DG = isfinite(G) & DF;
    % Surface = nombre de pixels correles
    S = sum(sum(DG));
    % Ecart GG-F
    G_F = G-F;
    % On annule l'effet des points hors domaine
    G_F(~DG) = 0;
    % Erreur actuelle
    Psi(cpt) = sum(sum(G_F.^2))/S;
  
    % Test oscillations
    if cpt>=2
        if Psi(cpt)>Psi(cpt-1);% | % Raffiner
            %if norm(dli + dli_0)<=10*epsilon(etape) % Psi(cpt)>Psi(cpt-1)% | % Raffiner
            warning('   Solution oscillante');
            % Bricolage : on sort avec les valeurs precedentes
            %li = li - dli;
            %Psi(cpt) =  Psi(cpt-1); 
            %dli = 0;
        end
    end
    
    % DISPLAYS --------------------------------------------------------
    
    
    % Bords du domaine sur image g
    figure(2);
	[p1,p2] = feval(champ,li,P1,P2,0);
	plot(p1,p2,[cG,'-'],'LineWidth',trait_gro);
    
    % Convergence
    figure(4); 
    semilogy(Psi','o-','MarkerSize',4*trait_gro,'LineWidth',trait_gro,'MarkerFaceColor','w');
    set(gca,'FontSize',myfontsize,'CLim',[0 0.25]);
    xlabel iteration;ylabel \psi;

    % Pour laisser le temps d'afficher
    pause(0.01);
    if cpt >=cptmax;warning ('    Nombre max d''iteration atteint');end
    cpt = cpt+1;

end

psi = Psi(end);
    
    
    
    
    
    
    
  