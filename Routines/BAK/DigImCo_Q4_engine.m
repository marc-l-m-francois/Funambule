function  li = DigImCo_Q4_engine(li,cptmax,precision);


global F DF g LX1 LX2 XN1 XN2 DN TailleEF


% REGLAGES ----------------------------------------------------------------


% Type d'interpollation retenu
type_interp	= 'linear'; 


% INITS GRAPHIQUES --------------------------------------------------------


global cG cF trait_fin trait_moy trait_gro myfontsize gamma_video


% INIT NEWTON  ------------------------------------------------------------


% Init compteur de boucles et ecart
cpt = 1;
% Champ actuel
[x1,x2,xn1,xn2] = Q4(li);
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

disp('   Calcul de la matrice'); % ne change pas en fonction des li
l = 0;  % indice de li
M = zeros(max(size(li)));
% Gradients
[gradF10,gradF20]  = gradient(F);
GF1GF1 = gradF10.*gradF10;
GF1GF2 = gradF10.*gradF20;
GF2GF2 = gradF20.*gradF20;

for p=1:size(XN1,1)
    disp(['   ',num2str(p),'/',num2str(size(XN1,1))]);
    for q=1:size(XN2,2)
        if DN(p,q)
            % Vecteur
            l = l+2;
            % Fonction de base centree au noeud
            LB1 = max(0,1 - abs(LX1-XN1(1,q))/TailleEF);
            LB2 = max(0,1 - abs(LX2-XN2(p,1))/TailleEF);
            m = 0;  % second indice de li
            for r=1:size(XN1,1)
                for s=1:size(XN2,2)
                    if DN(r,s)
                        m = m+2;
                        if m<=l % Symetrie OK
                            if abs(p-r)<=1 & abs(q-s)<=1 % Seult. noeuds voisins  
                                % Fonction de base centree au noeud
                                LC1 = max(0,1 - abs(LX1-XN1(1,s))/TailleEF);
                                LC2 = max(0,1 - abs(LX2-XN2(r,1))/TailleEF);
                                M(l-1,m-1) = (LB2.*LC2)*GF1GF1*(LB1.*LC1)';
                                M(l-1,m  ) = (LB2.*LC2)*GF1GF2*(LB1.*LC1)';
                                M(l  ,m  ) = (LB2.*LC2)*GF2GF2*(LB1.*LC1)';
                                M(l  ,m-1) = M(l-1,m);
                                M(m-1,l-1) = M(l-1,m-1);
                                M(m  ,l-1) = M(l-1,m  );
                                M(m-1,l  ) = M(l  ,m-1);
                                M(m  ,l  ) = M(l  ,m  );
                            end
                        end
                    end
                end
            end
        end
    end
end


% NEWTON ------------------------------------------------------------------


disp('   Calcul du minimum de correlation par methode de Newton');
% Inits
dli = inf*ones(size(li));
% Boucles
while max(abs(dli))>precision && cpt<cptmax
    cpt=cpt+1;disp(['   iteration ',num2str(cpt)]);
    % Elimination des points hors domaine
    gradF1 = gradF10;
    gradF2 = gradF20;
    gradF1(~DG) = 0;
    gradF2(~DG) = 0;
    % Calcul vecteur
    V = zeros(max(size(li)),1);
    l = 0;
    G_FG1 = G_F.*gradF1;
    G_FG2 = G_F.*gradF2;
    for p=1:size(XN1,1)
        for q=1:size(XN2,2)
            if DN(p,q)
                % Vecteur
                l = l+2;
                % Fonction de base centree au noeud
                LB1 = max(0,1 - abs(LX1-XN1(1,q))/TailleEF);
                LB2 = max(0,1 - abs(LX2-XN2(p,1))/TailleEF);
                V(l-1) = -LB2*G_FG1*LB1';
                V(l  ) = -LB2*G_FG2*LB1';
            end
        end
    end
    
    % Verification de la determination de la matrice
    detM = abs(det(M));
    if detM<1E-6
        warning(['      La matrice a determinant ',num2str(detM),' presque nul)']);
    end       

    % Correcteurs
    dli = (M\V)';
    li = li+dli;
    
    % Coo des noeuds dans repere de g
    xnp1 = xn1;xnp1(~DN)=NaN;
    xnp2 = xn2;xnp2(~DN)=NaN; 
    
    % Champ actuel
    [x1,x2,xn1,xn2] = Q4(li);
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
            disp('   Solution oscillante');
            % Bricolage : on calme le jeu
            %dli = (dli+dli_0)/2;
        end
    end
    
    % DISPLAYS --------------------------------------------------------
    
    % Image g
    Image_g(li);
    
    % Ecart F-G
    Image_G_F(G,DG);

    figure(4); 
    % Convergence
    semilogy(Psi','o-','MarkerSize',4*trait_gro,'LineWidth',trait_gro,'MarkerFaceColor','w');
    set(gca,'FontSize',myfontsize,'CLim',[0 0.25]);
    xlabel iteration;ylabel \psi;

    % Pour laisser le temps d'afficher
    pause(0.01);
    if cpt >=cptmax;warning ('    Nombre max d''iteration atteint');end

end
  




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

