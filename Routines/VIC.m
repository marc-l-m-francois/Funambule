function [L,Xc1,Xc2,f,psi,X1,X2,ok] = VIC(F,L,R,courbe)
% VIC : Virtual Image Correlation
% Applies for the precise measurement of a a curve or a silhouette
% Requires the parameterized curve in "courbe"
% to initially fit (with given parameters L) at list partially the curve
%
% INPUTS
%
% F image (real)
% L initial shape parameters
% R virtual image radius
% courbe (string) name of the curve
% sil_cont=true : silhouette sil_cont=false: contour (line)
% precision : max variation of L between two iterations
% type_interp = type of interpollation 
% 'nearest' (avoid),'linear' (best choice),'cubic','spline' (slow)
% nraff Mesh refinement. Best precision is reached for 3. Set 1 for faster computation
% cptmax : max number of iterations
% niv_display 
%    0 nothing
%    1 basic
%    2 more
%
% OUTPUTS
%
% L obtained shape parameters
% Xc1 Xc2 curve points
% f physical image in the frame x1 x2 of the virtual one
% psi value of the weight function
% X1 X2 virtual image computation points
% ok = true if everything is OK
%
% AUTHOR
%
% Marc L. M. FRANCOIS, laboratory GeM, Universite de Nantes, France
% marc.francois@univ-nantes.fr


% VIC SETTINGS ------------------------------------------------------------


global niv_display 
global sil_cont precision eps_eig motion_limitation cptmax

% User-chosen display level 0 nothing, 1 text, 2 real-time curve
if isempty(niv_display);niv_display=1;end

% Silhouette or contour
if isempty(sil_cont);
    sil_cont=logical(input('Silhouette [1] or contour [0] ? '));
end

% Precision defines the minimum of the curve parameters
if isempty(precision);precision=1E-3;end

% eps_eig defines the acceptable relative disrepency between eigenvalues
if isempty(eps_eig);eps_eig=1E-12;end

% Activate motion limitation in some problematic cases
if isempty(motion_limitation);motion_limitation=false;end

% Max number of iteration
if isempty(cptmax);cptmax=100;end

% Type of interpolation - linear is proven to be the best (contrary to DIC)
type_interp='linear';

% Degree of refinment - 3 is proven to be enough
nraff = 3;


% INITS -------------------------------------------------------------------


[GF1,GF2] = gradient(F);    % Gradient of the image
DL = inf*ones(size(L));     % Variation of L within each loop
cpt = 0;                    % Newton Iteration counter
psi = [];                   % Absolute error
dpsi= inf;                  % Rel difference of psi between two loops
numind = 1:max(size(DL));   % for displays

while  abs(dpsi) > precision/1000 && cpt<=cptmax % max(abs(DL))>=precision
    cpt = cpt+1;
    ok = true;
    
    
    % VIRTUAL IMAGE -------------------------------------------------------
    
    
    % Virtual image mesh
    if cpt ==1 % Else: computation on previous values
        [Xc1 ,Xc2] = feval(courbe,L,[0:1/1000:1]);
    end
    longueur = sum(sqrt( diff(Xc1).^2+diff(Xc2).^2 ));
    x1 = [0:1/ceil(nraff*longueur):1];
    x2 = [-1:1/ceil(nraff*R):1]';
    % The virtual image
    if sil_cont
        g = (ones(size(x2))+x2)*ones(size(x1))/2;
    else
        g = abs(x2)*ones(size(x1));
    end
    % These are area covered by point, not equal to diff(x1) or diff(x2)
    dx1 = 1/size(x1,2);
    dx2 = 2/size(x2,1);
    % Curvilinear abscissa x1 variationnal term
    ddx1 = min(diff(x1))/1000;
    % Curve and derivatives / parameters L
    [Xc1 ,Xc2 ,dXc1dL ,dXc2dL ] = feval(courbe,L,x1     );
    [XcP1,XcP2,dXc1dLp,dXc2dLp] = feval(courbe,L,x1+ddx1);
    % Numerical derivatives in x1
    dXc1dx1 = (XcP1-Xc1)/(ddx1);
    dXc2dx1 = (XcP2-Xc2)/(ddx1);
    ndXcdx1 = sqrt(dXc1dx1.*dXc1dx1 + dXc2dx1.*dXc2dx1);
    % Numerical second derivatives in x1
    d2Xc1dLdx1  = (dXc1dLp-dXc1dL)/ddx1;
    d2Xc2dLdx1  = (dXc2dLp-dXc2dL)/ddx1;
    % Tangent vector es and normal vector er to the curve
    es1 = dXc1dx1./ndXcdx1;
    es2 = dXc2dx1./ndXcdx1;
    er1 = es2;
    er2 =-es1;
    % Test 
    if min(ndXcdx1)==0
        error('The curve shows a stationnary point !');
    end
    % Mapping : current points X = Xc1 (of the curve) + R x2 er
    X1 = ones(size(x2))*Xc1 + R*x2*er1;
    X2 = ones(size(x2))*Xc2 + R*x2*er2;
    % Image f = image F in the frame of virtual image g
    f = interp2(F,X1,X2,type_interp,1);
    
    
    % DERIVATIVES ---------------------------------------------------------
    
    
    % Eventually set hereafter 'linear' to be faster...
    dFdX1 = interp2(GF1,X1,X2,type_interp,0);
    dFdX2 = interp2(GF2,X1,X2,type_interp,0);
    % Computation of reccurent terms in the loops
    TR = d2Xc1dLdx1.*(ones(size(L'))*er1./ndXcdx1) + ...
         d2Xc2dLdx1.*(ones(size(L'))*er2./ndXcdx1);
    TR1 = R*TR.*(ones(size(L'))*es1);
    TR2 = R*TR.*(ones(size(L'))*es2);
    % Initializations
    dpsidL = zeros(size(L'));
    d2psidL2 = zeros(size(L,2),size(L,2));
    % Vector and matrix computation
    for p=1:size(L,2)
        dX1dLp = ones(size(x2))*dXc1dL(p,:) - x2*TR1(p,:);
        dX2dLp = ones(size(x2))*dXc2dL(p,:) - x2*TR2(p,:);
        % Vector
        dpsidL(p) = sum(sum((dFdX1.*dX1dLp+dFdX2.*dX2dLp).*(f-g)))*dx1*dx2;
        % Erreur globale % sum(phi)*dx1;
        psi(cpt) = sum(sum((f-g).^2))*dx1*dx2/2;
        % Matrix
        for q=p:size(L,2)
            dX1dLq = ones(size(x2))*dXc1dL(q,:) - x2*TR1(q,:);
            dX2dLq = ones(size(x2))*dXc2dL(q,:) - x2*TR2(q,:);
            d2psidL2(p,q) = sum(sum((dFdX1.*dX1dLp+dFdX2.*dX2dLp).*...
                            (dFdX1.*dX1dLq+dFdX2.*dX2dLq)))*dx1*dx2;
            d2psidL2(q,p) = d2psidL2(p,q);
        end
    end
    
	% CORRECTOR DL -------------------------------------------------------- 


    % Inversion with too small eigenvalues removed
    % in place of simple DL = -([p]inv(d2psidL2)*dpsidL)'*R/(nraff);
    % Computation with removing the eigenvalues
    [V,D] = eig(d2psidL2);
    Dm = max(abs(D(D~=0)));
    M = zeros(size(d2psidL2));
    ind = false(size(DL));
    for p=1:numel(L)
        if abs(D(p,p)) > Dm*eps_eig
            M = M+V(:,p)*V(:,p)'/D(p,p);
        else
            if abs(D(p,p)) > 0
                M = M+V(:,p)*V(:,p)'/(sign(D(p,p))*Dm*eps_eig);
            end
            ind(p) = true;
            % ok = false; % Possible to stop with a null eigenvalue
        end
    end
    if niv_display>=2 && any(ind)
        disp(['    ',num2str(sum(ind)),' null eigenvalue'])
    end
    if all(ind)
        error('All the eigenvaluees are null');
    end

    % Corrector of the parameters L *(R/nraff) is the width of the border
    % in f
    DL = -(M*dpsidL)'*(R/(nraff));

    % Possible motion limitation - in case of problems
    if motion_limitation
        ind = abs(DL)>R;
        if any(ind)
            % DL(ind) = R*DL(ind)./max(abs(DL(ind)));
                        DL = R*DL/max(abs(DL));
            ok = false;
            if niv_display>=2
                disp(['    Motion limitation on parameter ',num2str(numind(ind))])
            end
        end
    end
    
    % New shape parameters
    L = L+DL;
    % Relative difference of psi
    if cpt>=2
        dpsi = (psi(cpt-1)-psi(cpt))/psi(cpt-1);
    end

    % Indicator of max number of iterations
    if cpt>=cptmax
        ok = false;
    end
    
    
    % INFORMATIONS --------------------------------------------------------
    

    if niv_display>=1 
        % disp(['    Psi ',num2str(psi(cpt)),', max(abs(DL)) ',num2str(max(abs(DL))),' / ',num2str(precision),'   L= [',num2str(L),']']);
        disp(['    Psi ',num2str(psi(cpt)),'   L= [',num2str(L),']']);
    end
    
    if niv_display>=3
       disp(['----------- internal VIC results loop ',num2str(cpt),' ------------------']);
        disp('Second derivatives, eigenvalues, eigenvectors')
        [V,D] = eig(d2psidL2);
        disp(d2psidL2);
        disp([D*ones(size(L,2),1)]');
        disp(V);
        disp('First derivatives')
        disp(dpsidL)
        disp('Corresponding corrector')
        disp(DL)
        disp('Initial / Final L')
        disp([L-DL;L])
        % Intermediate isplays
        figure(1);
        plot([X1(1  ,:);Xc1;X1(end,:)]',[X2(1  ,:);Xc2;X2(end,:)]',['r--'],'LineWidth',1);
        input('------------------------- go on ? -------------------------');
    end
end
