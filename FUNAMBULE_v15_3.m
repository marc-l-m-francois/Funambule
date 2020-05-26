disp(' ');
disp('-------------------------------------------------------------------------');
disp('|2008                                                               2020|')
disp('|                                                                       |')
disp('|                             Funambule V15.3                           |')
disp('|                                                                       |')
disp('|Marc Francois                                    Labo. GeM Univ. Nantes|')
disp('|marc.francois@univ-nantes.fr                                           |');
disp('-------------------------------------------------------------------------');
disp(' ');
disp('REFERENCES');
disp('[1] Francois, Semin, Auradou, 2010')
disp('    Identification of the shape of curvilinear beams and fibers. App. Mech. Mat. 24-25');
disp('[2] Semin, Auradou, Francois, 2011')
disp('    Accurate measurement of curvilinear shapes by virtualimage correlation, EPJ Ap 56, 2011');
disp('[3] Rethore, Francois, 2013');
disp('    Curves and boudaries measurements using B-splines and virtual images');
disp('    Opt. and Lasers in Engng., 52, 2011');
disp('[4] Francois, refused in Exp. Mech.')
disp('    Virtual Image Correlation uncertainty');


% MATLAB INITS ------------------------------------------------------------


clear all
path(path,'Routines');
path(path,'Routines_graphiques');
path(path,'Curves');
path(path,'Curves/BEAMS');

 
% GENERAL CHOICES ---------------------------------------------------------


% Display level 0 nothing, 1 standard, 2 real-time curve, 3 with checkings
niv_display = 1;
% Automatic mode (auto-sets the width, the luminance, the precision)
auto = true;
% To keep the same image automatically
keepimage = false;


% VIC SETTINGS ------------------------------------------------------------


global sil_cont precision eps_eig motion_limitation cptmax

% sil_cont : silhouette or contour - is asked below

% Precision defines the minimum of the curve parameters
% Initial value below is actualized during computation
precision=0.1;

% eps_eig defines the acceptable relative disrepency between eigenvalues
% Its defaut value can be modified in 'sensitive' curves i.e. Chainette.m
eps_eig = 1E-12;

% Possible limitation of the motion - activate in problematic cases
motion_limitation = false;

% Maximum iteration number for current calculus
cptmax = 20;

% Linear has been tested and proved to be the best for VIC
type_interp = 'linear';
% Degree of mesh refinement. Test showed that 3 is enough
nraff = 3;


% GRAPHICAL SETTINGS ------------------------------------------------------


global cF cG cD trait_gro trait_moy trait_fin myfontsize
cG = 'r';   % Virual image color
cF = 'b';   % Physical image color
cD = 'g';   % mu curve color
trait_gro = 0.70;
trait_moy = 0.35;
trait_fin = 0.18;
myfontsize = 12;


% IMAGE INITS -------------------------------------------------------------


% FIGURE 1 Physical image F
figure(1);clf;Range_mon_trace(1,1,1,2,4);
set(1,'Name','Image F and curve C','MenuBar','none','GraphicsSmoothing','off','ToolBar','none');
colormap gray;

% FIGURE 2 Physical image f in the frame of the virtual image g
figure(2);clf;Range_mon_trace(2,1,2,2,4);
set(2,'Name','Image f (unwrapped) and unconstrained identification mu');

% FIGURE 3 Statistics on local identifier mu
figure(3);clf;Range_mon_trace(3,2,1,2,4);
set(3,'Name','FFT(mu)');

% FIGURE 4 Convergence
if niv_display>=2
    figure(4);clf;Range_mon_trace(4,2,2,2,4);
    set(4,'Name','Convergence');
end


% USER CHOICES ------------------------------------------------------------


if ~keepimage
    % Image selection and preparation
    disptitle('USER CHOICES');
    disp('Select your image');
    [nom_F,ImagePath] = uigetfile('*.*','your choice:');
    IMA = imread(fullfile(ImagePath,nom_F));
    disp(['    Chosen image : ',nom_F]);
    [F,bitdepth] = PrepImage(IMA); clear IMA;
    save Temporaire/F0 F nom_F bitdepth
else
    load Temporaire/F0
end
figure(1);clf;set(gca,'Position',[0 0 1 1]);
imagesc(F,[0 1]);axis equal off;axis([0.5 size(F,2)+0.5 0.5 size(F,1)+0.5]);hold on

% Curve choice
disp('Choose the type of curve between:');
cd Curves/;a = dir;cd ..
cpt=0;
for p=1:size(a,1)
    if strcmp(a(p).name(end),'m') && ~a(p).isdir && size(a(p).name,2)>2
        if strcmp(a(p).name(end-1),'.')
            cpt = cpt+1;
            disp(['    ',num2str(cpt),' ',a(p).name(1:end-2)]);num(cpt) = p;
        end
    end
end
rep = input('Your choice: ');
courbe = a(num(rep)).name(1:end-2);
sil_cont = logical(input('Silhouette [1] or contour (curve) [0] : '));

% CURVE INITIALIZATION ----------------------------------------------------


disptitle('CURVE INITIALIZATION');
L = feval(courbe);
% Drawing the curve
[Xc1,Xc2] = feval(courbe,L,[0:1000]/1000);
plot(Xc1,Xc2,'r-');
% Additional drawings
feval(courbe,L);
% Initialization of the width R
if auto
    R = 2;
else
    % Virtual image initial half width (could be automatized)
    disp('Click a point to define the initial half width of the virtual image');
    disp('   wide enough for the virtual image to contain the silhouette / contour');
    [x1,x2] = ginput(1);
    R = sqrt(min( (Xc1-x1).^2 + (Xc2-x2).^2 ));
    if sil_cont         % Case silhouette
        R = max(R,2);
    end
    disp(['    R = ',num2str(R)]);
end

% Problem of orientation
if sil_cont
    disp('Default orientation is:');
    disp('    For closed curves: black (inside) on white (outside)');
    disp('    For open curves: white on the left and black on the right');
    rep = logical(input('    Is it true [1] or false [0] ? '));
    if ~rep
        F=1-F;
    end
end


% GRADIENT CHECK (Debug) --------------------------------------------------


if niv_display>=2
    disptitle('ANALYTIC GRADIENTS CHECK');
    precisionnum = 1E-3;
    x1 = [0:0.1:1];
    Lc = 100*(rand(size(L))-0.5); % Pixel amplitude
    [Xc1 ,Xc2 ,dXc1dL ,dXc2dL ] = feval(courbe,Lc,x1);
    for p=1:size(Lc,2)
        stop = false;
        dL = zeros(size(Lc)); dL(p) = precisionnum;
        [Xc1p ,Xc2p ,dXc1dLp ,dXc2dLp ] = feval(courbe,Lc+dL,x1);
        dXc1dLpNUM = (Xc1p-Xc1)/precisionnum;  dXc2dLpNUM = (Xc2p-Xc2)/precisionnum;
        if rms( abs(dXc1dLpNUM-dXc1dLp(p,:)))>precisionnum*numel(Lc)
            disp(['    gradient number ',num2str(p),', component 1, is FALSE: see analytic / numeric']);
            disp(dXc1dLp(p,:));
            disp(dXc1dLpNUM);
            stop = true;
        end
        if  rms( abs(dXc2dLpNUM-dXc2dLp(p,:)))>precisionnum*numel(Lc)
            disp(['    gradient number ',num2str(p),', component 2, is FALSE: see analytic / numeric']);
            disp(dXc2dLp(p,:));
            disp(dXc2dLpNUM);
            stop = true;
        end
        if ~stop           
            disp(['    gradient number ',num2str(p),' OK']);
        end
    end
end


% CALCULUS ----------------------------------------------------------------


stop = false;   % Flag for the end of the computation
nstep = 0;
cptmax0 = cptmax;

while ~stop
    stop = true;
    nstep = nstep+1;
    
    % VIC COMPUTATION -----------------------------------------------------

    
    % Call the main routine of the program
    if niv_display>=1
        disptitle (['COMPUTATION STEP ',num2str(nstep)]);
    end
    % First calculus letting the curve uncahnged for width auto setting
    if nstep == 1 && auto
        cptmax=1;
        [LL,Xc1,Xc2,f,psi,X1,X2,stop] = VIC(F,L,R,courbe);
    else
        cptmax = cptmax0;
        [L,Xc1,Xc2,f,psi,X1,X2,stop] = VIC(F,L,R,courbe);    
    end
    % Curve length
    longueur = sum(sqrt( diff(Xc1).^2+diff(Xc2).^2 ));

    
    % ANALYSIS OF f -------------------------------------------------------

    
    % Virtual image coordinates
    x1 = ( 0:1/(size(f,2)-1): 1) ;
    x2 = (-1:2/(size(f,1)-1): 1)';

    % Over-contrasting f (does not modify F !)
    if min(min(f))<0.5 && max(max(f))>0.5
        fwhite = mean(mean(f(f>0.5)));
        fblack = mean(mean(f(f<0.5)));
    else
        fwhite = max(max(f));
        fblack = min(min(f));
    end
        
    if fwhite~=fblack;ff = (f - fblack)/(fwhite - fblack);else ff=f;end
    ff = min(ff,1);   ff = max(ff,0);

    % Local identifier mu
    if sil_cont                                 % silhouette - analytical formula
        dx2 = 2/size(f,1);                      % Area covered by point
        mu = -sum(ff-0.5)*dx2;
    else                                        % contour (curve)
        % f(f>0.1)=1;
        mu = (x2'*(1-ff))./(sum(1-ff)+1E-12);     % center of gravity
    end
    % its FFT 
    [mu_frequency,mu_amplitude] = Spectrum(mu,size(mu,2)-1);

    % Local with
    lwidth = sum(1-ff)/size(x2,1);
                                                % sum(ff<0.5)/size(x2,1); % fmg2<0.01 & 
    
    % Index of detected domains
    ind_det = lwidth>1/size(x2,1);
	rat_det = sum(ind_det/size(f,2));    
    % Index of correlated domains
    ind_cor = ind_det & abs(mu)<lwidth;
	rat_cor = sum(ind_cor/size(f,2));    
    % Mean width
    mwidth = mean(lwidth(ind_det));
    
    % Displays
    if niv_display>=2
        disp(' ');
        disp(['    Detected   ',num2str(rat_det*100),' %']);
        disp(['    Correlated ',num2str(rat_cor*100),' %']);
%        disp(['    Mean width ',num2str(R*mwidth),' pixels']);
    end
    
	% Noise of f (mean of each STD of correlated parts of f along x2)
    sigma_0 = mean(std(ff(:,ind_cor),0,2)); 
    % Quantification noise (classical formula)
    sigma_0q = 1/(2^bitdepth*sqrt(12));
    
    
    % FIGURE 1 physical image F -------------------------------------------
    
    
    if niv_display>=2 || stop     % Refresh image F
        figure(1);clf;set(gca,'Position',[0 0 1 1]);
        imagesc(F,[0 1]);axis equal off;axis([0.5 size(F,2)+0.5 0.5 size(F,1)+0.5]);hold on
        feval(courbe,L);        % Supplemenatary drawings
    end
    % Curve and virtual image borders
    if niv_display>=1
        figure(1);
        plot([X1(1  ,:);Xc1;X1(end,:)]',[X2(1  ,:);Xc2;X2(end,:)]',[cG,'-'],'LineWidth',trait_gro);
        zoom on
    end
    
    
    % FIGURE 2 unwrapped image f and local identifier mu ------------------
    
    
    if niv_display>=2 || stop
        figure(2);clf;
        % Image f - scaled at the pixel scale to help reading
        pcolor(longueur*[0:1/(size(f,2)-1):1],R*[-1:2/(size(f,1)-1):1],ff);shading interp;colormap gray
        xlabel('L x_1 (pixel)');ylabel ('R x_2 (pixel)');
        set(gca,'FontSize',myfontsize);hold on;
        plot([0;longueur],[0;0],[cG,'-'],'LineWidth',trait_gro);
        plot([0 0;longueur longueur],[-R R;-R R],[cG,'-'],'LineWidth',trait_moy);
        % Current virtual image points - for paper onmy (very slow !!!)
        % [xx1,xx2] = meshgrid([0:1/(size(f,2)-1):1],[-1:2/(size(f,1)-1):1]);
        % plot(longueur*xx1,R*xx2,[cG,'.'],'MarkerSize',trait_gro);
        % Pixel grid in the frame (avoid !!!)
        % TraceGrillePixelSurg;
        % Center (curve C) and virtual image borders
        % Local identifier mu
        if rat_det>0
            muplot = plot(longueur*x1(ind_det),R*mu(ind_det),[cD,'.'],'LineWidth',trait_fin,'MarkerSize',1);
            plot([0;longueur],[0;0],[cG,'-'],'LineWidth',trait_gro);
            legend(muplot,'local id. \mu detected')
        elseif rat_cor>0
            muplot = [muplot,plot(longueur*x1(ind_cor),R*mu(ind_cor),['b.'],'LineWidth',6*trait_gro,'MarkerSize',1)];
            legend(muplot,'local id. \mu detected','local id. \mu correlated')
        end
        
    end
    
    
    % THEORETICAL MEASUREMENT PRECISION -----------------------------------
    
    
    % Noise-induced (quadratic average is an assumption)
    sigma_n = sqrt(sigma_0^2+sigma_0q^2)*sqrt(2*size(L,2)/(R*longueur));
    % Discretization-induced
    sigma_d = size(L,2)/(20*longueur);

    
    % FIGURE 3 Spectrum of mu ---------------------------------------------
    
    
    if niv_display>=2 || stop
        figure(3);clf;
        % Trace of mu
        loglog(longueur./mu_frequency,R*mu_amplitude,'ko','MarkerFaceColor',cD,'LineWidth',trait_gro);

        xlabel('wavelength (pixel)'); ylabel('|FFT(R\mu(x_1))| (pixel)');
        myaxis = axis;
        myaxis(2) = 2*longueur;     % To show the infinite wave length
        hold on

        % Continuous term (infinite wave length)
        if mu_amplitude(1)<min(mu_amplitude(2:end))/10
            loglog(myaxis(2),min(R*mu_amplitude(2:end))/10,'kV','MarkerFaceColor',cD,'LineWidth',trait_gro);
            myaxis(3) = min(R*mu_amplitude(2:end))/10;
        else
            loglog(myaxis(2),R*mu_amplitude(1),'k>','MarkerFaceColor',cD,'LineWidth',trait_gro);
            myaxis(4) = max(R*mu_amplitude(1:end));
        end
        text(myaxis(2),myaxis(3),' \infty','rotation',90,'HorizontalAlignment','left','VerticalAlignment','bottom','color','k','FontSize',12)

        % Horizontal lines shows the sigma
        loglog([myaxis(1),myaxis(2)],[sigma_d,sigma_d],'r-','LineWidth',trait_gro);
        text(myaxis(2),sigma_d,'\sigma_d','HorizontalAlignment','right','VerticalAlignment','bottom','color','k','FontSize',12)
        loglog([myaxis(1),myaxis(2)],[sigma_n,sigma_n],'b-','LineWidth',trait_gro);
        text(myaxis(2),sigma_n,'\sigma_n','HorizontalAlignment','right','VerticalAlignment','bottom','color','k','FontSize',12)
        myaxis(3) = min([myaxis(3),sigma_d,sigma_n]);

        % Vertical bars show the representative wavelenght
        loglog(longueur*[1 1],[myaxis(3) myaxis(4)],'k-','LineWidth',trait_gro)
        text(longueur,myaxis(3),' Curve length  ','color','k','rotation',90,'verticalalignment','bottom','horizontalalignment','left')
        loglog(longueur*[1 1]/size(L,2),[myaxis(3) myaxis(4)],'k-','LineWidth',trait_gro);
        text(longueur/size(L,2),myaxis(3),' Curve length / number of parameters  ','color','k','rotation',90,'verticalalignment','bottom','horizontalalignment','left') 
        loglog([1 1],[myaxis(3) myaxis(4)],'k-','LineWidth',trait_moy);
        text(1,myaxis(3),' Pixel size','color','k','rotation',90,'verticalalignment','top','horizontalalignment','left')
        loglog(2*longueur*[1 1]/(size(mu,2)-1),[myaxis(3) myaxis(4)],'k-','LineWidth',trait_moy);
        text(2*longueur/(size(mu,2)-1),myaxis(3),' Sub-pixel','color','k','rotation',90,'verticalalignment','top','horizontalalignment','left')
        
        
        axis(myaxis)
        grid on
    end

    
   % FIGURE 4 Convergence -------------------------------------------------
   
   
    if niv_display>=2
        % Convergence
        figure(4);
        semilogy(psi,[cG,'o-'],'MarkerSize',4*trait_gro,'LineWidth',trait_gro,...
            'color','k','MarkerFaceColor','w');set(gca,'FontSize',myfontsize);
        xlabel iteration; ylabel \psi;grid on;set(gca,'FontSize',myfontsize);drawnow;
    end

    
    % CORRECTION of the half-width R --------------------------------------
    
    
    if auto
        R0 = R;
        % R are rounded at 0.1 pixel in order not to loop for nothing
        if rat_det==0
            R = max(2,ceil(10* (2*R) )/10);
            if R~=R0
                disp(['    No detection at all : increase of the wdth to R=',num2str(R)]);
                stop = false;
            end
        elseif rat_det<(1-1/size(L,2))
            R = max(2,ceil(10* (R*(1-0.1/size(L,2))/rat_det) )/10);
            if R~=R0
                disp(['    Incomplete detection : increase of the wdth to R=',num2str(R)]);
                stop = false;
            end
       elseif rat_cor<(1-0.1/size(L,2))
            R = max(2,ceil(10* (R*(max(abs(mu)) + mwidth)) )/10);
            if R~=R0
                disp(['    Complete detection but incomplete correlation: R=',num2str(R)]);
                stop = false;
            end            
        else
            if sil_cont     % Case silhouette
                R = max(2,ceil(10* (R*(max(abs(mu)-0.5) + mwidth)) )/10);   % 2 Optimal value (neat borders)
            else            % Case contour
                R = max(2,ceil(10* (2*R*mwidth)) /10);
            end
            if R~=R0
                disp(['    Complete correlation: new R=',num2str(R)]);
                stop = false;
            end
        end
    else
        rep = input(['    R=',num2str(R),' New value for R or keep it [] ']);
        if ~isempty(rep);R = rep;stop=false;end
    end
    
    % PRECISION CORRECTION ------------------------------------------------
    
    
    preciso0 = precision;
    % Eveluated prcision
    if rat_cor>(1-0.1/size(L,2))
        precision = max(sigma_n,sigma_d);
        if abs((preciso0 - precision)/preciso0)>0.1
           disp(['    Convergence criterion set at: ',num2str(precision),' pixel']);
           stop = false;
        end
    end

    
    % LUMINANCE CORRECTION ------------------------------------------------

    
    if sil_cont % Silhouette
        fwhite = mean(f(end,ind_cor));
        fblack = mean(f(  1,ind_cor));
    else        % Contour
        fwhite = mean([f(  1,ind_cor),f(end,ind_cor)]);
        fblack = mean( f( (size(f,1)+1)/2,ind_cor ) );
    end
    
    if auto
        if rat_cor>0 % rat_det>0.99
            if ~(abs(fwhite-1)<precision && abs(fblack)<precision)
                % Linear optimal correction
                F = F/(fwhite - fblack) - fblack/(fwhite - fblack);
                disp(['    Linear luminance correction']);
                stop = false;
            end
        end
    else
        rep = logical(input('    Linear correction of the luminance [1] or not [] or [0] ? '));
        if isempty(rep);rep=false;end
        if rep
                 F = F/(fwhite - fblack) - fblack/(fwhite - fblack);
                 stop = false;
        end
    end
    
    if niv_display>=3
        input('Let''s continue ?');
    end
end


% RESULTS -----------------------------------------------------------------


disptitle('RESULTS');
disp(['Parameters found for the ''',courbe,''' to fit the image ''',nom_F,''':']);
feval(courbe,L,[0:1/size(Xc1,2):1],1);


% ANALYSIS OF THE VALIDITY OF THE RESULTS FROM MU -------------------------


disptitle('UNCERTAINTY')
disp(['Theoretical uncertainty sigma_n associated to image noise: ',num2str(sigma_n),' pixel']);
disp(['Theoretical uncertainty sigma_d due to discretisation    : ',num2str(sigma_d), ' pixel']);
disp(' ');

% ANALYSIS OF MU FOR USER --------------------------------------------------


ecart_the = max(sigma_d,sigma_n);

% Infinite wavelength - bias
ecart_mes = R*mu_amplitude(1);
if ecart_mes>ecart_the
    fprintf(2,'Warning: for inifinite wavelenghts (bias) - see Figures 2 and 3.\n');
    disp(['    Distance from mu to VIC is ',num2str(100*(ecart_mes-ecart_the)/ecart_the),' % above the theoretical value']);
    disp('    In case of a silhouette, this may be due to an ill-contrasted image');
end

% Wavelengths between lenght of the inflence domain to length
ecart_mes = R*max(mu_amplitude(mu_frequency>0 & mu_frequency<=size(L,2)));
if ecart_mes>ecart_the
    fprintf(2,'Warning: for long wavelenghts (between the curve lenght and the length par parameter) - see Figures 2 and 3.\n');
    disp(['    Distance from mu to VIC is ',num2str(100*(ecart_mes-ecart_the)/ecart_the),' % above the theoretical value']);
    disp('    This reveals an unadapted curve');
end

% Wavelenghts between the pixel size (1) to lenght of the inflence domain
ecart_mes = R*max(mu_amplitude(mu_frequency>size(L,2) & mu_frequency<=longueur));
if ecart_mes > ecart_the
    fprintf(2,'Warning: for medium wavelenghts (between 1, the pixel size, and the length par parameter) - see Figures 2 and 3.\n');
    disp(['    Distance from mu to VIC is ',num2str(100*(ecart_mes-ecart_the)/ecart_the),' % above the theoretical value']);
    disp('    This reveals an ill border imaging');
end

ecart_mes = R*max(mu_amplitude(mu_frequency>longueur));
if ecart_mes > ecart_the
    fprintf(2,'Warning: for short wavelenghts (below the pixel size) - see Figures 2 and 3.\n');
    disp(['    Distance from mu to VIC is ',num2str(100*(ecart_mes-ecart_the)/ecart_the),' % above the theoretical value']);
    disp('    This should not occur');
end


% COMPLEMENTARY INFORMATIONS ----------------------------------------------


disptitle('COMPLEMENTARY INFORMATION:')
disp(['Estimated image noise: ',num2str(sigma_0*100),' %']);
disp(['Quantification noise for ',num2str(bitdepth),' bits: ',num2str(sigma_0q*100),' %']);
disp(['The curve length is: ',num2str(longueur),' pixel']);
disp(['The half-width of the virtual image was ',num2str(R),' pixels']);
disp('Optimal parameters are in the variable L');
disp('Curve points are in variables Xc1 and Xc2');
disp('Identified curve parameters L =');
disp(L);

disptitle('THANK YOU')


