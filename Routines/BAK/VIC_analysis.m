function [F1,R1] = VIC_analysis(F,f,Xc1,Xc2,L,R,sil_cont);


% SETTINGS ----------------------------------------------------------------


% Acceptable distance from local measurement mu to global VIC measurement
MAXIPEAK = 0.1;

global cF cG cD trait_gro trait_moy trait_fin myfontsize


% MU CURVE (UNCONSTAINED FIT)  --------------------------------------------


% Retrieving internal values
dx1 = 1/size(f,2);
dx2 = 2/size(f,1);
x1=[0:size(f,2)-1]/(size(f,2)-1);
x2=[-(size(f,1)-1):2:(size(f,1)-1)]/(size(f,1)-1);
longueur = sum(sqrt( diff(Xc1).^2+diff(Xc2).^2 ));

% Mu curve (local identification)
if sil_cont
    % Exact
    mu = -R*sum(f-0.5)*dx2;
else
    % Dirty
    mu = -R*sum(f-0.5)*dx2/2;
    mu = mu-mean(mu);
end
% Drawing mu
figure(2);
plot(longueur*x1,mu,[cD,'-'],'LineWidth',trait_moy);
plot([0;longueur],[0;0],[cG,'-'],'LineWidth',trait_gro);


% COMPARIZON TO UNCONSTRAINED CURVE MU  -----------------------------------


% FFT(mu)
fftmu = abs(fft(mu));
indm = floor(size(fftmu,2)/2);
fftmu = [fftmu(1),2*fftmu(2:indm)];
[maxfft,indmaxfft] = max(fftmu/size(x1,2));
freqmu = x1(1:indm)*size(x1,2)/longueur;
% longdondmu = 1./freqmu;
% Statistics on mu
% disp(['The avarage of mu, the local identifier, is ',num2str(mean(mu)),...
%     ' and its RMS is ',num2str(rms(mu)),' pixels']);

if maxfft>MAXIPEAK
    fit = false;
    disp('    This curve cannot fit the contour or silhouette of interest');
    ratio = x1(indmaxfft)*size(x1,2);
    disp('    The mu curve (figure 3) shows:');
    disp(['    a peak of wavelength = curve length / ',num2str(ratio),...
          ' (frequency : ',num2str(x1(indmaxfft)*size(x1,2)/longueur),' pixel-1),']);
    disp(['    of power ',num2str(maxfft),' pixels.']);
    disp('    Suggestion : choose a richer, more flexible, curve.');
else
    fit = true;
    disp('    The curve seems to fit perfectly the contour / silhouette.');
end

% FIGURE 3 FFT(mu)
figure(3);clf;
%semilogx(longdondmu,fftmu/size(x1,2),'ko','MarkerFaceCOlor',cD,'LineWidth',trait_gro);
plot(freqmu,fftmu/size(x1,2),'ko','MarkerFaceCOlor',cD,'LineWidth',trait_gro);
%axis([0 longueur 0 max(MAXIPEAK,maxfft)]);
set(gca,'FontSize',myfontsize);
%xlabel('wavelength (pixel)');
xlabel('frequency (pixel)-1');
ylabel('|FFT(\mu)| (pixel)');
set(gca,'FontSize',myfontsize);hold on
plot([freqmu(1) freqmu(end)],MAXIPEAK*[1 1],'r-');
axis on
grid on
% coolplot(15,3,'FFTmu'); % Bouding box adaptee [38 0 390 82] % ICICI


% LIGHT ANALYSIS AND PROPOSAL ---------------------------------------------


F1 = [];
% Light analysis
if sil_cont
    fblack = mean(f(1  ,:));    % bord int
    fwhite = mean(f(end,:));    % bord ext
else
    fwhite = mean([f(1  ,:),f(end,:)]);     % bords int et ext
    fblack = mean(f(ceil(size(f,1)/2),:));  % milieu
end

if maxfft<R && (abs(fblack)>1E-3 | abs(fwhite-1)>1E-3)
    disp('    Image F brightness appears not ideal. You can use the proposed optimal correction.');
    % Linear optimal correction
    F1 = F/(fwhite - fblack) - fblack/(fwhite - fblack);
    % Possible limitation to [0 1] ==> loss of information
    % F = min(F,1); F = max(0,F);
end


% VIRTUAL IMAGE WIDTH ANALYSIS AND PROPOSAL -------------------------------


R1 =[];
if fit 
    if ~sil_cont
        if mean(mean(f))>0.55 |  mean(mean(f))<0.45
            R1 = R*0.5/mean(mean(f));
            disp(['    Suggestion : choisir un R final = ',num2str(R1),'.']);
        else
            disp(['    The retained width R = ',num2str(R),' is perfect.']);
        end
    else
        if R~=1.5
            R1 = 1.5;
            disp(['    Optimal value for neat image = ',num2str(R1)]);
        end
    end
end







