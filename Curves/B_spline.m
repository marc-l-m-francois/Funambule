function [Xc1,Xc2,dXc1dL,dXc2dL] = B_spline(L,x1,rien);
% x1 parametre C [0 1]
% L parametres de controle ----> [P1x,P1y,...,Pnx,Pny]

global ordre_B_spline fermee

switch nargin    
    
	case 0 % Init ---------------------------------------------------------
        
        ordre_B_spline = input('    B-spline order (classically 2) ? ');
        fermee = logical(input('    Closed [1] or open [0] ? '));
        disp('    Click control points (enter after the last one)');
        disp('    in case of silhouette, click from x to y (CLOCKWISE !)');
        [x,y] = ginput;
		% Affectation des points en parametres
		for p=1:size(x,1)
				L(2*p-1) = x(p);
				L(2*p  ) = y(p);
		end
		% X=Vk Parametres optimaux depuis les points selectionnes		
		Xc1 = L;

    case 1 % Supplementary drawings ---------------------------------------
        for p=1:size(L,2)/2
            % Drawing control points
            plot(L(2*p-1),L(2*p),'r+');
            text(L(2*p-1),L(2*p),[' P',num2str(p)],'FontSize',18,'Color','r');
        end
   
    case 2 % Calculus -----------------------------------------------------

		% Inits
		nL = size(L,2);
		% Control points
		Q = [L(1:2:nL-1);L(2:2:nL)];
        if ordre_B_spline>size(Q,2);
            error ('    Order too high or to few points');
        end
		% Closing by repetition of last points
		if fermee
			Q = [Q,Q(:,1:ordre_B_spline)];
        else
            %for d=1:ordre_B_spline-1
            %    Q = [Q(:,1),Q,Q(:,end)];
            %end
        end
		nbq = size(Q,2);
		% Nombre de segments de la B-Spline
		nbk = size(x1,2)-1;
		% Vecteur de noeud (regulier pour B-Spline uniforme)
		k = 1:(nbq+ordre_B_spline+1);
		% Parametre
		t = ordre_B_spline+1+x1*(nbq+1-ordre_B_spline-1);
		% Ordre 0
		N = zeros(nbq+ordre_B_spline,nbk+1);
		for l=1:nbq+ordre_B_spline
				if k(l)~=k(l+1)
						N(l, t>=k(l) & t<k(l+1) ) = 1;
				end
        end
		% Ordres >0
		for d=1:ordre_B_spline
				% Init nouveau polynome
				M = zeros(nbq+ordre_B_spline-d,nbk+1);
				for l=1:nbq+ordre_B_spline-d
						% Fonction multiplicatrices
						if k(l+d)~=k(l)
								f = (t-k(l  ))/(k(l+d)-k(l));
						else
								f = zeros(size(t));
						end
						if k(l+1+d)~=k(l+1)
								g = (k(l+1+d)-t)/(k(l+1+d)-k(l+1));
						else
								g = zeros(size(t));
						end
						% Nouveaux polynomes
						M(l,:) = f.*N(l,:) + g.*N(l+1,:);
				end
				N = M;
		end 
		% Courbe
   		x  = Q*N;
		Xc1 = x(1,:);
		Xc2 = x(2,:);
        
        % Derivees
        dXc1dL = zeros(nL,nbk+1);
        dXc2dL = zeros(nL,nbk+1);
        
        for p=1:size(N,1)
            pp = mod(p-1,nL/2);
            dXc1dL(2*pp+1,:) = dXc1dL(2*pp+1,:) + N(p,:);
            dXc2dL(2*pp+2,:) = dXc2dL(2*pp+2,:) + N(p,:);
        end
   
    case 3 % Informations- ------------------------------------------------

		disp(['    B-Spline of order ',num2str(ordre_B_spline)]);
        disp('    Control point coordinates:');
        disp([L(1:2:end-1)',L(2:2:end)']);
        

end


