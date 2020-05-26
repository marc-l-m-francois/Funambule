%%% VIC_analysis2

MAXIPEAK = 0.1;


% MU THE UNCONSTRAINED FITTING --------------------------------------------


dx2 = 2/size(f,1);
if sil_cont
    mu = -R*sum(f-0.5)*dx2;
else                        % Dirty calculus !
    mu = -R*sum(f-0.5)*dx2/2;
    mu = mu-mean(mu);
end
figure(2);
plot(longueur*x1,mu,[cD,'-'],'LineWidth',trait_moy);
plot([0;longueur],[0;0],[cG,'-'],'LineWidth',trait_gro);


% ANALYSIS OF MU ----------------------------------------------------------



% Analysis of mu
fs = size(x1,2)-1;
[freqmu,fftmu] = Spectrum(mu,fs);

figure(3);clf;
% plot(freqmu,fftmu,'ko','MarkerFaceCOlor',cD,'LineWidth',trait_gro);
semilogx(freqmu,fftmu,'ko','MarkerFaceCOlor',cD,'LineWidth',trait_gro);

set(gca,'FontSize',myfontsize);
xlabel('frequency x_1^{-1}');
ylabel('|FFT(\mu(x_1))| (pixel)');
set(gca,'FontSize',myfontsize);hold on
% plot([freqmu(1) freqmu(end)],MAXIPEAK*[1 1],'r-');
axis on
grid on
% coolplot(15,3,'FFTmu'); % Bouding box adaptee [38 0 390 82] % ICICI
disp('Fig.3 spectral analysis of mu.');
disp('Feq. 0 : bias. Freq. 1 to 100 can be corrected by richer curve.');
disp('High frequency are related to pixelization.');


% ADVISE ON THE NATURE OF THE CURVE ---------------------------------------


maxfft = max(fftmu);
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
