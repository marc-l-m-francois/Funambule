function [XN1,XN2,DN,DF] = Mailleur_EF(F,TailleEF,P1,P2);
% F image de travail (on ne se sert que de sa taille)...
% TailleEF taille des elements finis carres
% PA P2 coordonnees des points definissant le domaine de calcul
%
% TXN1,TXN2 coo des noeuds des EF
% DN Domaine de definition des nouds : ceux en dehors de P1, P2 sont false
% DF Domaine de definition de F (correspond a l'interieur du maillage)


% Activer pour les verifications
debuggage = false;

        
    disp('   Creation du maillage');

	% Element centre au centre de l'image (autant que faire se peut)
    xn1 = [ceil((size(F,2)-floor(size(F,2)/TailleEF)*TailleEF)/2):TailleEF:size(F,2)];
    xn2 = [ceil((size(F,1)-floor(size(F,1)/TailleEF)*TailleEF)/2):TailleEF:size(F,1)];

    % Noeud centre au mileu de l'image (souci coo pas entiere)
%     k = floor( (size(F,2)-1)/(2*TailleEF) );
%     XN1 = [ (size(F,2)+1)/2-k*TailleEF:TailleEF:(size(F,2)+1)/2+k*TailleEF ];
%     k = floor( (size(F,1)-1)/(2*TailleEF) );
%     XN2 = [ (size(F,1)+1)/2-k*TailleEF:TailleEF:(size(F,1)+1)/2+k*TailleEF ];
    % Maillage partant de 1,1 (ballot)
%     XN1 = [1:TailleEF:size(F,2)];
%     XN2 = [1:TailleEF:size(F,1)];
    
    % Finalisation
    [XN1,XN2] = meshgrid(xn1,xn2);
    
    % Supression du dernier point double pour l'affichage
    P1 = P1(1:end-1);
    P2 = P2(1:end-1);

    % DN Table d'appartenance des noeuds au domaine de definition (coo noeuds)
    if size(P1,1)==4 && all(P1(1:4)==[1 1 size(F,2) size(F,2)]' & ...
        P2(1:4)==[1 size(F,1)  size(F,1) 1]')
        % Cas ou le domaine de definition est toute l'imagette
        DN = true(size(XN1));
    else
        DN = DetectoDedans(XN1,XN2,P1,P2);
    end

    % DF Table d'appartenance au domaine de def
    DF = false(size(F));

    % Selection des noeuds qui definissent un element carre
    ndn1 = size(DN,1);
    ndn2 = size(DN,2);
    for p=1:ndn1
        for q=1:ndn2
            
            
            % Indicateur de l'appartenance de ce noeud a un element carre
            ok = false;
            % Si ce noeud appartient a priori au domaine
            if DN(p,q)
                % Existence d'un voisin en bas a droite
                if p<ndn1 && q<ndn2 && ...
                    DN(p+1,q  ) && DN(p+1,q+1) && DN(p  ,q+1)
                    DF(xn2(p):xn2(p)+TailleEF,xn1(q):xn1(q)+TailleEF) = true;
                    ok = true;
                end
                % Existence d'un voisin en haut a droite
                if p>1 && q<ndn2 && ...
                    DN(p  ,q+1) && DN(p-1,q+1) && DN(p-1,q  )
                    DF(xn2(p)-TailleEF:xn2(p),xn1(q):xn1(q)+TailleEF) = true;
                    ok = true;
               end
                % Existence d'un voisin en haut a gauche
                if p>1 && q>1 && ...
                    DN(p-1,q  ) && DN(p-1,q-1) && DN(p  ,q-1)
                    DF(xn2(p)-TailleEF:xn2(p),xn1(q)-TailleEF:xn1(q)) = true;
                    ok = true;
                end
                % Existence d'un voisin en bas a gauche
                if p<ndn1 && q>1 && ...
                    DN(p  ,q-1) && DN(p+1,q-1) && DN(p+1,q  )
                    DF(xn2(p):xn2(p)+TailleEF,xn1(q)-TailleEF:xn1(q)) = true;
                    ok = true;
                end
                % Suppression du noeud s'il n'apparait dans aucun element
                DN(p,q) = ok;
            end
        end
    end
    

% Trace de DF pour debug
if debuggage
    [X1,X2]=meshgrid([1:size(F,2)],[1:size(F,1)]);
    X1(~DF) = NaN;X2(~DF) = NaN;
    plot(X1,X2,'w.');
end

% Verification
if sum(sum(DN))==0
    error('   -> Domaine trop petit ou elements trop grand ! ');
end


