% Created by Justin Cooke
% Purpose of this script is to conduct a spectral analysis on the smooth
% portion of the dune field so that we may determine a proper cutoff
% frequency

clc;
clear;
close all;

%% Load the Data -- Original Data

myDir = dir('./Smooth/May1_LogSpaced/z*');

N = length(myDir);

cr = 0;
cx = 0;

for i = 1:N
   
    tempFile = myDir(i).name;
    tempFolder = myDir.folder;
    
    tempName = strcat(tempFolder,'/',tempFile);
    
    if contains(tempName,'README')
        xyz = load(tempName);
        xyz = [xyz(301,2:end),xyz(468,2:end),xyz(635,2:end),xyz(801,2:end)];
        cr = cr + 1;
        loc_data{cr} = xyz;
    elseif contains(tempName,'comp')
        u_x = load(tempName);
        u_x = [u_x(:,301),u_x(:,468),u_x(:,635),u_x(:,801)];
        cx = cx + 1;
        ux_data{cx} = u_x;
    end
    
    
end

clear temp* myDir cx cr urms_x xyz u_x

%% Now load in new smooth data

myDir = dir('./Smooth/Combined/z*comp*');

N = length(myDir);

for i = 1:N
   
    tempFile = myDir(i).name;
    tempFolder = myDir.folder;
    
    tempName = strcat(tempFolder,'/',tempFile);
    
    tempux = load(tempName);
    tempux = tempux(:,4:end);
    ux_data{i} = [ux_data{i} ; tempux];

end

clear temp* myDir 

%% Choose an x location (x = 1000, 1250, 1500, 1750)

store_ux_data = ux_data;

count = 0;
for i = 1:N
    if i ~= 2 && i ~= 3 && i ~= 5
        count = count + 1;
        ux_data_trunc{count} = store_ux_data{i}(:,4);
    end
end

%% FOR XY PLANE


dt = 0.65;

[Nt,Nx] = size(ux_data_trunc{1});


count = 0;
z = zeros(N-3,1);
for i = 1:N
    if i ~= 2 && i ~= 3 && i ~= 5
        count = count + 1;
        z(count) = loc_data{i}(1,3);
    end
end

Nz = length(z);


Lt = dt*Nt;
dw  = 2*pi/Lt;
n_w  = -Nt/2 : 1 : Nt/2 - 1;
omega = n_w.*dw;

u_tau = 0.1313; % Val at 1805
nu = 0.32e-4;
deltanu = nu/u_tau;
delta = 300;

clear u_bar

% Creating mean and fluctuating parts
u_bar{Nz,1}   = [];
u_rms{Nz,1}   = [];
lam_x{Nz,Nx}   = [];
lam_xp{Nz,Nx}  = [];
omega_p{Nz,1} = [];
kx{Nz,Nx}      = [];

for i = 1:Nz
    u_bar{i} = mean(ux_data_trunc{i});
    u_rms{i} = ux_data_trunc{i} - u_bar{i};

    % Now create streamwise wavelength
    for j = 1:Nx
        lam_x{i,j} = u_bar{i}(j)./omega;
        lam_xp{i,j} = (u_bar{i}(j)/u_tau)./(omega.*(deltanu/u_tau));
        kx{i,j} = (2*pi)./lam_x{i,j};
    end
    
    omega_p{i} = omega.*(deltanu/u_tau);
    
end

u_bar = cell2mat(u_bar);

t = linspace(0,Lt,Nt);
T = t.*(u_bar(end)/delta);
total_T = T(end);

% %% Conduct FFT
% 
% for i = 1:Nz
%       urms_temp = u_rms{i};
% 
%       u_hat_w = (1/Nt) * fftshift(fft(urms_temp));
%       psi_w = abs(u_hat_w).^2;
%       E_hat_w = psi_w ./ (dw);
%          
%     % Pre-multiplied Energy Spectra
%       E_kx = E_hat_w.*u_bar(i); % Taylor's hypothesis to go E(w) -> E(kx)
%       kxE_uu = E_kx.*kx{i}';
%      
%       E_w{i} = E_hat_w;
%       kxE_kx{i} = kxE_uu; 
%       
% end
% 
% clear urms_temp psi_w E_hat_w E_uu u_hat_w



%% Plot Spectra

% kxEw = cat(2,kxE_kx{:});
% 
% close all;
% 
% zz = z./delta;
% zplus = z./deltanu;
% parfor i = 1:Nz
%     ll{i} = lam_x{i}./delta;
% end
% 
% for i = 1:Nz
%     [Z,L{i}] = meshgrid(zz,ll{i});
%     [ZP,LP{i}] = meshgrid(zplus,lam_xp{i});
% end
%     
% close all;
% 
% figure()
% plot3((Z),(L{end}),((kxEw./(u_tau^2)))); hold on;
% % zlim([0 20]);
% % ylim([0.0 5]);
% zlabel('kx\cdotE_{uu}/u^2_\tau');
% ylabel('\lambda/\delta');
% xlabel('z^+');
% grid on;
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% figure()
% contour(ZP,LP{end},kxEw./(u_tau^2),20); hold on;
% 
% grid on;
% grid minor;
% ax1 = gca;
% set(ax1,'xscale','log');
% set(ax1,'yscale','log');
% xlabel('z/\delta');
% ylabel('\lambda/\delta');
% colorbar;
% caxis([0,20]);



%% Create Large Scales

lxp_cutoff = 70000;
lx_cutoff = 340; %lxp_cutoff*deltanu;

half = Nt/2;

for i = 1:Nz
    
    omega_cutoff = u_bar(i)/lx_cutoff;
    omega_pcutoff = omega_cutoff*(deltanu/u_tau);
    
    up_rms{i} = u_rms{i}./u_tau;
    
    uhat_w{i} = (1/Nt) * fftshift(fft(up_rms{i}));
    
    ii_l(i) = find(abs(omega_p{i}) <= omega_pcutoff,1);
    ii_u(i) = half + (half - ii_l(i)) + 2;
    
    temp = uhat_w{i};
    temp(1:ii_l(i) - 1) = 0;
    temp(ii_u(i) + 1:end) = 0;
    
    uhat_wcutoff = real(temp);
    
    uhat_S = (uhat_w{i}) - uhat_wcutoff;
    uL{i} = Nt * ifft(ifftshift(uhat_wcutoff));
    uS{i} = Nt * ifft(ifftshift(uhat_S));
    
end

%% Now do Hilbert Transform and Stuff

for i = 1:Nz
    
    tempsmall = real(uS{i});

    ush{i} = hilbert(tempsmall);
    ushi = imag(ush{i});
    ushr = real(ush{i});
    Env_uS{i} = sqrt((ushi.^2) + (ushr.^2));
    
    % Now filter the envelope
    
    store_me = (1/Nt) * fftshift(fft(Env_uS{i}));
    store_me(1:ii_l - 1) = 0;
    store_me(ii_u + 1:end) = 0;
    
    EnvL_uS{i} = Nt * ifft(ifftshift(store_me));
    
    
end

%%

close all;

zl_ref = 5;
zs_ref = 5;

figure()
tiledlayout(3,1)
p1 = nexttile;
plot(T,up_rms{zs_ref});
ylabel('u^+','Interpreter','tex','FontName','SansSerif','FontSize',16);
% myTitle = strcat('Z_L= ',num2str(round(z(zl_ref))),' & Z_S= ',...
%     num2str(round(z(zs_ref))),' & \lambda_c =',num2str(lxp_cutoff));
% title(myTitle);
ylim([-10,10]);
yline(0);
text(0.25,7.5,sprintf('A'),'FontName','SansSerif','FontSize',16);
set(gca,'FontSize',16);

p2 = nexttile;
ap1 = plot(T,uL{zl_ref},T,EnvL_uS{zs_ref}-mean(EnvL_uS{zs_ref}));
ap2 = yline(0);
ylabel('u^+_L','Interpreter','tex','FontName','SansSerif','FontSize',16);
legend([ap1(1),ap1(2)],'u^+_L','E_L(u^+_S)','','Interpreter','tex',...
    'FontName','SansSerif','FontSize',16,'Box','off','Orientation','Horizontal');
ylim([-1,1]);
set(gca,'FontSize',16);
text(0.25,.75,sprintf('B'),'FontName','SansSerif','FontSize',16);

p3 = nexttile;
ap3 = plot(T,uS{zs_ref},T,Env_uS{zs_ref});
ap4 = yline(0);
ylabel('u^+_S','Interpreter','tex','FontName','SansSerif','FontSize',16);
legend([ap3(1),ap3(2)],'u^+_S','E(u^+_S)','','Interpreter','tex',...
    'FontName','SansSerif','FontSize',16,'Box','off','Orientation',...
    'Horizontal','Location','SouthEast');
ylim([-10,10]);
text(0.25,7.5,sprintf('C'),'FontName','SansSerif','FontSize',16);

linkaxes([p1,p2,p3],'x')
xticklabels([p1,p2],{});
xlabel('tU_\infty / \delta','Interpreter','tex','FontName','SansSerif','FontSize',16);
xlim([0,30]);
set(gca,'FontSize',16);

clear ap* ax1 p1 p2 p3




%% AM Coeff Calculation - Single Point

% zref = zl_ref; % 50m approx.

clear uis* Env*uis uoL*

for i = 1:Nz
    
    uoL = real(uL{i});
    uis = real(uS{i});
  
    uish = hilbert(uis);
    uishi = imag(uish);
    uishr = real(uish);
    
    Env_uis = sqrt((uis.^2) + (uishi.^2));
    
    store_me = (1/Nt) * fftshift(fft(Env_uis));
    store_me(1:ii_l(i)-1) = 0;
    store_me(ii_u(i)+1:end) = 0;
    
    store_me = real(store_me);
    
    EnvL_uiS = Nt * ifft(ifftshift(store_me));
    EnvL_uis{i} = EnvL_uiS - mean(EnvL_uiS);
    
    top_prod = uoL.*real(EnvL_uis{i});
    top = mean(top_prod)';
    
    bot = rms(uoL)*rms(real(EnvL_uis{i}));
    
    R(i) = top/bot;
%     R(i) = real(R(i));
    
end

clear temp* clear bot1 bot 2

%%
delta=300;
centerLL = 3.9*sqrt(delta/deltanu);
zplus = z./deltanu;
zdelta = z./300;

close all;


% open('/home/jpcooke/Downloads/twopoint_wall/Images/SinglePointCoeff.fig');
% semilogx(zplus,R); hold on;
% xline(centerLL,'--');
% myTitle = strcat('One-Point AM Coefficient: ',num2str(round(total_T)),'T');
% title(myTitle);        
% xlabel('z/\delta');
% ylabel('R');
% grid on;
% grid minor;
% legend('Expt: Re_\tau = 15,000','3.9\surd(15,000)',...
%     '3.9\surd(1.23\times10^6)','Numerical Data');


%----------- NOTE: SCALING IS 30 not 300 so as to see how it looks 

open('/home/jpcooke/Downloads/twopoint_wall/Images/SinglePointCoeff_zdelta.fig'); hold on;
% open('/home/jpcooke/Downloads/twopoint_wall/Images/SinglePointCoeff.fig');
semilogx(z/70,R); 
% xline(centerLL,'--');
myTitle = strcat('One-Point AM Coefficient: ',num2str(round(total_T)),'T');
title(myTitle);        
xlabel('z/\delta');
ylabel('R');
grid on;
grid minor;


%% Two Point Correlation

% Find time difference
% 
% myTime = zeros(Nz,1);
% II = zeros(Nz,1);
% sig3{Nz,1} = [];
% 
% for i = 1:Nz 
%    
%     sig1 = fft(ux_data_trunc{i});
%     sig2 = fft((ux_data_trunc{1}));
%     sig2 = conj(sig2);
%     
%     sig3{i} = ifft(sig1.*sig2);
%     
%     myMax = max(sig3{i});
%     II(i) = find(sig3{i} == myMax,1);
%     
%     myTime(i) = T(II(i));
%     
% 
% end
% 
% close all;
% 
% figure();
% plot(T,sig3{6})
% 
% clear myMax 
% 
% %% Two-Point R_AM
% 
% clear uis* Env*uis uoL*
% 
% zref = 18;
% for i = 1:Nz
%     
%     uoL = uL{i};
%     uis = uS{1};
%     
% %     dtau = (z(i) - z(4)) / (u_bar(end)*tand(14.5));
% %     
% %     if i == 3
%         ind_tc(i) = 0;
% %     else
% %         %ind_tc(i) = find(t >= abs(dtau),1);
% %         ind_tc(i) = find(T >= myTime(3),1);
% %     end
% %     
%     if i <= 3
%         ul = uoL(1:Nt-ind_tc(i));
%         us = uis(ind_tc(i)+1:end);
%     elseif i > zref
%         ul = uoL(ind_tc(i)+1:end);
%         us = uis(1:Nt-ind_tc(i));
%     end
%         
%     uish = hilbert(us);
%     uishi = imag(uish);
%     uishr = real(uish);
%     
%     Env_uis = sqrt((us.^2) + (uishi.^2));
%     
%     store_me = (1/(Nt-ind_tc(i))) * fftshift(fft(Env_uis));
%     store_me(1:ii_l(i)-1) = 0;
%     store_me(ii_u(i)+1:end) = 0;
%     
%     EnvL_uiS = (Nt-ind_tc(i)) * ifft(ifftshift(store_me));
%     EnvL_uis{i} = EnvL_uiS - mean(EnvL_uiS);
%     
%     top_prod = ul.*real(EnvL_uis{i});
%     top = mean(top_prod)';
%     
%     bot = rms(ul)*rms(EnvL_uis{i});
%     
%     R_tp(i) = top/bot;
%     
% end
% 
% clear temp* clear bot1 bot 2
% 
% 
% 
% centerLL = 3.9*sqrt(1000/deltanu);
% zplus = z./deltanu;
% 
% close all;
% 
% figure();
% semilogx(zplus,R_tp,'-o'); hold on;
% % xline(centerLL);
% myTitle = strcat('Two-Point AM Coefficient: \lambda^+_c= ',num2str(lxp_cutoff));
% title(myTitle);     
% xlabel('z^+');
% ylabel('R');


%% Ask to continue the analysis

prompt = 'Do you want to continue the analysis? (Y/N) \n';

yn = input(prompt,'s');

if contains(yn,'Y') || contains(yn,'y')
    ...
elseif contains(yn,'n') || contains(yn,'N')
    return
else
    disp('Invalid choice');
    return
end

%% Load Dune Data

myDir = dir('./DuneField/Combined/x_*');

df_N = length(myDir);

cx = 0;
cr = 0;

for i = 1:df_N
    
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        df_xyz = load(loadMe);
        df_xyz = df_xyz(:,2:end);
        cx = cx + 1;
        df_loc{cx} = df_xyz;
    else
        df_ux = load(loadMe);
        df_ux = df_ux(:,4:end);
        cr = cr + 1;
        df_uxdata{cr} = df_ux;
    end
    
    
end

clear myName myFolder loadMe df_xyz df_ux cx cr

myDir = dir('./DuneField/WSS/WallStress*');

surf_N = length(myDir);

for i = 1:surf_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        surf_xyz = load(loadMe);
        surf_xyz = surf_xyz(:,2:end);
    else
        surf_tau = load(loadMe);
        surf_tau = surf_tau(:,4:end);
    end
    
    
end

clear myName myFolder loadMe

myDir = dir('./DuneField/WSS/SWSS*');

surf_N = length(myDir);

cx = 0;
cr = 0;

for i = 1:surf_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        swss_xyz = load(loadMe);
        swss_xyz = swss_xyz(:,3:end);
        cx = cx + 1;
        swss_y{cx} = swss_xyz;
    else
        swss_tau = load(loadMe);
        swss_tau = swss_tau(:,4:end);
        cr = cr + 1;
        swss_tauw{cr} = swss_tau;
    end
    
    
end

clear myName myFolder loadMe cr cx swss_xyz swss_tau

% Note that df_uxdata is set up as the cell # represents the x location,
% the size of each cell is Nt x Nz

%% Wall Shear Stress
 
% Streamwise wall stress averaging

surf_tau_bar = mean(surf_tau);

for i = 1:df_N/2
   
    x(i) = df_loc{i}(1,1);
    loc_z{i} = df_loc{i}(:,end);
    
    II = find(surf_xyz(:,1) > x(i),1);
    
    tau(i) = surf_tau_bar(II);
    df_utau(i) = sqrt(tau(i)/1.23);
    df_dnu(i) = nu./df_utau(i);
    
    df_Retau(i) = delta/df_dnu(i);

end

% Spanwise station based averaging

for i = 1:df_N/2
   
    swss_taubar(i) = mean(swss_tauw{i},'All');
        
end

df_taubar = mean(surf_tau);

close all;

figure();
plot(surf_xyz(:,1),df_taubar,'k'); hold on
plot(x,swss_taubar,'rdiamond','MarkerFaceColor','r');
% plot(surf_xyz(:,1),...
%     (2.*mean(surf_tau))./(1.23*u_bar(end)*u_bar(end)),'k'); hold on;
% plot(x,(2.*tau)./(1.23*u_bar(end)*u_bar(end)),'rdiamond');
xlabel('x');
ylabel('\tau_w');
xlim([0 8000]);

x_hat = x - 1850;

figure();
tiledlayout(2,1);
ax1 = nexttile;
plot(surf_xyz(:,1)-1850,mean(surf_tau),'k'); hold on;
plot(x_hat,swss_taubar,'rdiamond','MarkerFaceColor','r');
ylabel('\tau_{wall}');
ax2 = nexttile;
plot(surf_xyz(:,1)-1850,surf_xyz(:,end),'k');
ylabel('z');
% for i = 1:df_N/2
%     myText = strcat(num2str(round(x(i))),' [m]');
%     lt = xline(x(i),'--k',{myText});
% end
linkaxes([ax1,ax2],'x');
xticklabels(ax1,{})
myplot.TileSpacing = 'compact';
xlabel('x-x_0');
xlim([-1846 6300]);



figure();
yyaxis left
plot(surf_xyz(:,1)-1850,mean(surf_tau),'Color',"#D95319"); hold on;
plot(x_hat,swss_taubar,'kdiamond','MarkerFaceColor',"blue");
ylabel('\tau_{wall} [Pa]');
ax1 = gca;
ax1.YColor = "#D95319";
yyaxis right
area(surf_xyz(:,1)-1850,surf_xyz(:,end),'FaceColor',"#808080",...
    'FaceAlpha',0.2,'EdgeColor',"#808080");
ax2 = gca;
ax2.YColor = "#808080";
ylim([0 20]);
ylabel('z [m]');
xlabel('x-x_0');
xlim([-1846 6300]);

%% Dune specific data

[df_Nt,df_Nz] = size(df_uxdata{1});
df_Nx = df_N/2;

df_z = zeros(N,1);
for i = 1:N
   df_z(i) = loc_data{i}(1,3);
end

df_n = -df_Nt/2 : 1 : df_Nt/2 - 1;
df_Lt = df_Nt*dt;
df_t = linspace(0,df_Lt,df_Nt);
df_dw = (2*pi)/df_Lt;
df_omega = df_n.*df_dw;

for i = 1:df_Nx

    df_ubar{i} = mean(df_uxdata{i});
    df_uprime{i} = df_uxdata{i} - df_ubar{i};
    
    df_upplus{i} = df_uprime{i}./df_utau(i);
    
end



%% Create Large Scales

df_half = df_Nt/2;

for i = 1:df_Nx

    df_utz = df_upplus{i}; % u'(z,t) at single x-loc
    df_ubarz = df_ubar{i}; % u_bar(z) at single x-loc
    df_lx_cutoff = lx_cutoff; %delta_iblLocal(i); %lxp_cutoff*df_dnu(i); % lambda_x cutoff val at x-loc
    df_omega_p = df_omega.*(df_dnu(i)/df_utau(i)); % create omega_p at each x-loc
    
    for j = 1:df_Nz % Now we are isolating z at each x-location
       
        df_ut = df_utz(:,j); % u'(t) 
        df_ubz = df_ubarz(j); % u_bar is a constant at each z
        df_uhat = (1/df_Nt) * fftshift(fft(df_ut)); % uhat(omega)
        
        omega_cutoff = df_ubz/df_lx_cutoff; % using Taylor's hypothesis to get cutoff value for omega
        omegap_cutoff = omega_cutoff*(df_dnu(i)/df_utau(i)); % converting to omega^+
                
        df_ii_l{i}(j) = find(abs(df_omega_p) <= omegap_cutoff,1);
        df_ii_u{i}(j) = df_half + (df_half - df_ii_l{i}(j)) + 2;

        temp = df_uhat;
        temp(1:df_ii_l{i}(j) - 1) = 0;
        temp(df_ii_u{i}(j) + 1:end) = 0;

        df_uhat_L = temp;
        
        temp_uL{j} = df_Nt * ifft(ifftshift(df_uhat_L));
        temp_uS{j} = df_ut - temp_uL{j};
        
        
    end
    
    mat_uL = cell2mat(temp_uL);
    mat_uS = cell2mat(temp_uS);
    
    df_uL{i} = mat_uL;
    df_uS{i} = mat_uS;
    
end

clear mat_u* temp* df_utz df_ubarz df_ut df_ubz df_uhat 


%% Now we will conduct amplitude modulation analysis 

for i = 1:df_Nx
    
    ulzt = df_uL{i};
    uszt = df_uS{i};
    ix_l = df_ii_l{i};
    ix_u = df_ii_u{i};
    
    for j = 1:df_Nz
       
        ult = ulzt(:,j);
        ust = uszt(:,j);
        
        ust_h = hilbert(ust);
        ust_hi = imag(ust_h);
        ust_hr = real(ust_h);
        
        Eust = sqrt((ust_hr.^2) + (ust_hi).^2);
        
        Eusw = (1/df_Nt) * fftshift(fft(Eust));
        
        temp = Eusw;
        temp(1:ix_l(j)-1) = 0;
        temp(ix_u(j)+1:end) = 0;
        
        ELusw = temp;
        ELust_temp = real(df_Nt * ifft(ifftshift(ELusw)));
        ELust = ELust_temp - mean(ELust_temp);
        
        topprod = ult.*ELust;
        top = mean(real(topprod));
        bot = rms(ult)*rms(ELust);
        
        R_temp(j) = top/bot;
        
    end
    
    df_R{i} = R_temp;
    
end

clear temp* clear ELust_temp %topprod top bot

%% DF Single-Point Plot

close all;

figure();

for i = 2:df_Nx
    
    nexttile;
    semilogx(zdelta,R); hold on;
    semilogx(df_z/delta,df_R{i}); hold on
    yline(0,'--k');
    ylim([-0.8,0.5]);
    xlabel('z/\delta');
    ylabel('R(z/\delta)');
    myTitle = strcat(num2str(round(x(i))),' [m]');
    title(myTitle);
    legend('Smooth Playa','Dunes','Location','NorthEast');
end

% open('/home/jpcooke/Downloads/twopoint_wall/Images/SinglePointCoeff_zdelta.fig');
% % open('/home/jpcooke/Downloads/twopoint_wall/Images/SinglePointCoeff_zdelta_150T.fig');
% semilogx(zdelta,R); hold on;
% myTitle = strcat('One-Point AM Coefficient');
% title(myTitle);        
% xlabel('z/\delta');
% ylabel('R');
% grid on;
% grid minor;

%% Two Point Analysis


for i = 1:df_Nx
    
    ulzt = df_uL{i};
    uszt = df_uS{i};
    ix_l = df_ii_l{i};
    ix_u = df_ii_u{i};
    
    for j = 1:df_Nz
       
        ult = ulzt(:,15);
        ust = uszt(:,j);
        
        ust_h = hilbert(ust);
        ust_hi = imag(ust_h);
        ust_hr = real(ust_h);
        
        Eust = sqrt((ust_hr.^2) + (ust_hi).^2);
        
        Eusw = (1/df_Nt) * fftshift(fft(Eust));
        
        temp = Eusw;
        temp(1:ix_l(j)-1) = 0;
        temp(ix_u(j)+1:end) = 0;
        
        ELusw = temp;
        ELust_temp = df_Nt * ifft(ifftshift(ELusw));
        ELust = ELust_temp - mean(ELust_temp);
        
        topprod = ult.*ELust;
        top = mean(topprod);
        bot = rms(ult)*rms(ELust);
        
        R_temp(j) = top/bot;
        
    end
    
    df_R2{i} = R_temp;
    
end

clear temp* clear ELust_temp %topprod top bot

%% DF Two-Point Plot

% close all;
% 
% figure();
% 
% for i = 1:df_Nx
%     
%     nexttile;
% %    semilogx(zdelta,R_tp); hold on;
% %     semilogx(df_z/delta,df_R2{i}); hold on
%     yline(0,'--k');
%     ylim([-0.5,0.5]);
%     xlabel('z/\delta');
%     ylabel('R(z/\delta)');
%     myTitle = strcat(num2str(round(x(i))),' [m]');
%     title(myTitle);
%     legend('Smooth Playa','Dunes','Location','NorthEast');
% end

%% Dune Field Mean Velocity Profiles

close all;

for i = 1:df_Nx
    df_zplus{i} = df_loc{i}(:,end)./df_dnu(i);
end

figure();

for i = 1:df_Nx
    
    nexttile;
    semilogx(df_zplus{i},df_ubar{i}./df_utau(i)); hold on;
    semilogx(zplus,u_bar./u_tau);
    xlabel('z^+');
    ylabel('u^+');
    myTitle = strcat(num2str(round(x(i))),' [m]');
    title(myTitle);
    ylim([0,60]);
    xlim([10^3,5*10^6]);
    legend('Dunes','Smooth Playa','Location','NorthWest');
end

figure();

for i = 1:df_Nx
    
    nexttile;
    plot(df_ubar{i},df_z); hold on;
    plot(u_bar,z);
    ylabel('z');
    xlabel('u');
    myTitle = strcat(num2str(round(x(i))),' [m]');
    title(myTitle);
%     ylim([0,60]);
    xlim([3,6.5]);
    legend('Dunes','Smooth Playa','Location','NorthWest');
end

%% Plot <u'u'>

for i = 1:df_Nx
    
    for j = 1:df_Nz
   
        temp = df_upplus{i}(:,j);
        store_dfuu{j} = temp.*temp;
        
    end
    
    mat_dfuu = cell2mat(store_dfuu);
    df_upup{i} = mat_dfuu;
    df_upupbar{i} = mean(df_upup{i});
    
end

figure();

for i = 1:df_Nx
    
    nexttile;
    semilogx(df_zplus{i},df_upupbar{i}); hold on;
    xlabel('z^+');
    ylabel('\langle u^\prime u^\prime \rangle^+');
    myTitle = strcat(num2str(round(x(i))),' [m]');
    title(myTitle);
    ylim([0,15]);
    xlim([10^3,5*10^6]);
%     legend('Smooth Playa','Dunes','Location','NorthEast');
end

%% Teresa Saxton-Fox IBL Calc Method

myDir = dir('./DuneField/TSF_IBL/x*');

ntemp = length(myDir);

c1 = 0;
c2 = 0;
c3 = 0;

for i = 1:ntemp
   
    datatemp = load(strcat(myDir(i).folder,'/',myDir(i).name));
    
    if contains(myDir(i).name,'u_avg')
        c1 = c1 + 1;
        uxbar_tsf{c1} = mean(datatemp(:,4:end));
    elseif contains(myDir(i).name,'README')
        c2 = c2 + 1;
        x_tsf{c2} = datatemp(:,2);
        z_tsf{c2} = datatemp(:,end);
    else
        c3 = c3 + 1;
        ux_tsf{c3} = mean(datatemp(:,4:end));
    end
    
end


clear dpyp df_dpyp max_dpyp idp ibl_tsf 

low = 1;
high = 100;

for i = 1:10
   
    mean_vel = (uxbar_tsf{i});
    local_z = z_tsf{i};
    global_z = z_tsf{i} - z_tsf{i}(1);
    U_infty = uxbar_tsf{end}(end);
    delta_0 = 340;
    
    for j = low+1:high-1
        yp = local_z(j);
        yi = local_z(low);
        yf = local_z(high);
        
        Uyp = mean_vel(j);
        Uyi = mean_vel(low);
        Uyf = mean_vel(high);
        
        pt_behind = (Uyp - Uyi)/sqrt(yp - yi); 
        pt_forward = (Uyf - Uyp)/sqrt(yf - yp);
        
        dpyp(j-(low)) = (pt_behind - pt_forward) * (sqrt(delta_0)/U_infty);
        
    end
    
    
    df_dpyp{i} = dpyp;
    
    max_dpyp(i) = max(dpyp);
    idx = find(df_dpyp{i} == max_dpyp(i),1);
    myY = local_z(idx);
    ibl_tsf(i) = myY;
    
end

close all;

figure();
for i = 1:10
nexttile;
plot(sqrt(z_tsf{i}/delta),ux_tsf{i}/U_infty,'o');
xlabel('\surd(y/\delta_0)');
ylabel('U/U_0');
xlim([0.25 0.98]);
end
% figure();
% plot(sqrt(z_tsf{1}(low:high)/delta),df_dpyp{1},'.');
% xlabel('\surd(y/\delta_0)');
% ylabel('\Delta_p');

for i = 1:10
   df_xVals(i) = x_tsf{i}(1) - x_tsf{1}(1); 
end

% based on eyeing it 
ibl_tsf_est = [0,46.74,57.3,63.6,71,98.05,109.7,123.9,159.2,165.3];
% ibl_tsf(5) = ibl_tsf_est(5);
ibl_tsf(1) = 0;

ft = fittype('a*x^b');
f3 = fit(df_xVals',ibl_tsf',ft,'StartPoint',[0,0]);

figure();
plot(df_xVals,ibl_tsf,'k*'); hold on;
plot(f3,'k--');
text(4000,30,sprintf('a_{loc} = %f\nb_{loc} = %f',f3.a,f3.b));
xlabel('x [m]');
ylabel('\delta_{ibl}');

%% Li method for delta ibl calculation



for i = 1:df_Nx
   
    % Calculate u'u'/u_0^2 first
    temp = df_uprime{i};
    uu_u0u0{i} = (temp.*temp)./(U_infty^2);
    uu_u0u0bar{i} = mean(uu_u0u0{i});
    
    delta_0 = 300;
    
    x_hat(i) = x(i) - 1850;
    xh_delta0(i) = x_hat(i)/delta_0;
    log_xhd(i) = log10(xh_delta0(i));
    
end

for i = 2:df_Nx
   
    duu = uu_u0u0bar{i} - uu_u0u0bar{i-1};
    
    dlogxhd = log_xhd(i) - log_xhd(i-1);
    
    ratio{i-1} = duu/dlogxhd;
    
    
end

%% Plot IBL

close all;

figure();

xhatdelta = [1.22,2.45,3.89,5.56,7.52,9.79,12.45,15.55,19.17];

for i = 1:2:df_Nx-2
    
    nexttile;
    semilogx(df_z./delta,ratio{i},'ko','MarkerFace','red'); hold on;
    %yline(5*10^-5);
    %semilogx(df_loc{i}(:,end)./delta,ratio{i},'o');
    xlabel('z/\delta_{ABL}','Interpreter','tex','FontSize',16,'FontName','SanSerif');
    ylabel('a(z)','Interpreter','tex','FontSize',16,'FontName','SanSerif');
       
    set(gca,'FontSize',16,'FontName','SanSerif');
end

delta_ibl = [0,0.109,0.1943,0.1943,0.1602,0.1602,0.2355,0.2856,...
    0.4198,0.5089].*delta;
delta_iblLocal = [0.0074266667,0.1644,0.3654,0.3208,0.3565,0.3377,0.3817,0.52,...
    0.6128,0.7403].*delta;

ft = fittype('a*x^b');
f = fit(x_hat',delta_ibl'./300,ft,'StartPoint',[0,0]);
f2 = fit(x_hat',delta_iblLocal',ft,'StartPoint',[0,0]);
f3 = fit(x_hat',ibl_tsf_est',ft,'StartPoint',[0,0]);
 

figure();
plot(x_hat,delta_ibl,'k^','MarkerFace','blue'); hold on;
%plot(x_hat,delta_iblLocal,'ko','MarkerFace','red');
plot(x_hat,ibl_tsf,'ksquare','MarkerFace','red'); 
plot(f,'b');
%plot(f2,'r--');
plot(f3,'r-.');
text(-1750,200,sprintf('Li Method:\na = %f\nb = %f',f.a,f.b),'Color','blue',...
    'FontName','SansSerif','FontSize',14);
%text(250,0.55,sprintf('Li Local: a = %f\nb = %f',f2.a,f2.b),'Color','red',...
%    'FontName','SansSerif','FontSize',16);
text(-1750,125,sprintf('Parthasarathy Method:\na = %f\nb = %f',f3.a,f3.b),...
    'FontName','SansSerif','FontSize',14);
xlabel('x - x_0 [m]','FontName','SansSerif','FontSize',16);
ylabel('\delta_{ibl}','FontName','SansSerif','FontSize',16);
legend('Li \it{et al} Method','Parthasarathy and Saxton-Fox Method',...
    'Li \it{et al} Fit','Parthasarathy and Saxton-Fox Fit',...
    'location','northwest','Box','off','NumColumns',2,'Orientation','Horizontal');
set(gca,'FontSize',16,'FontName','SansSerif');
yticks([0 100 200 300]);
xticks([-2000 0 2000 4000 6000]);
xticklabels({'-2000','0','2000','4000','6000'});
xlim([-2000 6000]);

%% Single-Point Dune AM Normalized with delta_ibl

close all;

delta_ibl = [.1,0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*300;

% figure();
% 
% for i = 2:df_Nx
%     
%     nexttile;
%     semilogx(zdelta,R); hold on;
%     p1 = semilogx(df_z/delta_ibl(i),df_R{i}); hold on
%     p1.Color = "#2F4858";
%     p2 = semilogx(df_z/delta,df_R{i}); hold on
%     p2.Color = "#F6AE2D";
%     yline(0,'--k');
%     ylim([-0.8,0.8]);
%     xlabel('z/\delta');
%     ylabel('R(z/\delta)');
%     myTitle = strcat(num2str(round(x(i))),' [m]');
%     title(myTitle);
%     legend('Smooth Playa','Dunes with \delta_{IBL}',...
%         'Dunes with \delta_{ASL}','Location','NorthEast');
% end



% mySymbol = ["square",">","<","x","*","diamond","^",...
%     "o","v","*"];
% myColors = ["#04e4ff","#00dfdf","#00d7b5","#0acd84","#58c14e",...
%     "#81b100","#a59d00","#c78200","#e65d00","#ff0000"];

% myColors = ["#FFD700", "#DA7800", "#EB2427", "#DF0069", "#7A0245", ...
%     "#C876B5", "#BB00FF", "#6E00CF", "#0000FD", "#2C196A"];

myColors = [ "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000", "	#000000" ];

figure();
ax1 = subplot(1,2,1);
for i = 2:df_Nx
semilogx(df_z,df_R{i},'Color',myColors(i),'LineWidth',1.25); hold on;
end
semilogx(z,R,'k-.','LineWidth',1.5);
yline(0,'--k');
ylim([-0.8,0.5]);
xlim([0 400]);
xlabel('$z$','Interpreter','latex');
ylabel('$R(z)$','Interpreter','latex');
%legend('off')
% legend(strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(2)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(3)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(4)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(5)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(6)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(7)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(8)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(9)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(10)-1850)/300),2))),...
%     'Smooth Playa','NumColumns',2,'Location','NorthEast',...
%     'interpreter','latex');
ax1.FontSize = 16;

ax2 = subplot(1,2,2);
% open('/home/jpcooke/Downloads/twopoint_wall/Images/SinglePointCoeff_zdelta.fig'); hold on;
for i = 2:df_Nx
semilogx(df_z/delta_ibl(i),df_R{i},'Color',myColors(i),'LineWidth',1.25); hold on;
end
semilogx(z/70,R,'k-.','LineWidth',1.5);
% title('Collapse of R_{AM} with \delta_{ibl} in Dune Field');
yline(0,'--k');
ylim([-0.8,0.5]);
xlim([5*10^-3 10^1]);
xlabel('$z/\delta_{IBL}$','Interpreter','latex');
ylabel('$R(z/\delta_{IBL})$','Interpreter','latex');
% legend('off')
% legend(strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(2)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(3)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(4)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(5)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(6)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(7)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(8)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(9)-1850)/300),2))),...
%     strcat('$\hat{x}/\delta_{ASL}$ = ',num2str(round(((x(10)-1850)/300),2))),...
%     'Smooth Playa','NumColumns',1,'Location','EastOutside',...
%     'interpreter','latex');
legend(strcat('$\hat{x}/\delta$ = ',num2str(round(((x(2)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(3)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(4)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(5)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(6)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(7)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(8)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(9)-1850)/300),2))),...
    strcat('$\hat{x}/\delta$ = ',num2str(round(((x(10)-1850)/300),2))),...
    'Smooth Playa','NumColumns',5,'Location','NorthOutside',...
    'interpreter','latex');

ax2.FontSize = 16;

%% FIG 4 Updated

close all;

figure();
semilogx(z/30,R,'k','LineWidth',2); hold on;
for i = 2:df_Nx
    if i == 2 
        semilogx(df_z/30,df_R{i},'Color',myColors(i),'LineWidth',1.25);
    else
        semilogx(df_z/delta_ibl(i),df_R{i},'Color',myColors(i),'LineWidth',1.25); 
    end
end
yline(0,'--k');
ylim([-0.8,0.7]);
xlim([5*10^-3 10^1]);
xlabel('$z/\delta_{IBL}$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
ylabel('$R(z/\delta_{IBL})$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
% legend('$\hat{x}_{-1}$','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
%     '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$','$\hat{x}_9$',...
%     'NumColumns',5,'Location','NorthEast','Orientation','Horizontal',...
%     'interpreter','latex','FontName','SansSerif','FontSize',16,'Box','off');
set(gca,'FontName','SansSerif','FontSize',16);

%%
figure();
for i = 2:df_Nx
semilogx(df_z,df_R{i},'Color',myColors(i),'LineWidth',1.25); hold on;
end
semilogx(z,R,'k-.','LineWidth',1.5);
yline(0,'--k');
ylim([-0.8,0.7]);
xlim([0 400]);
xlabel('$z$','Interpreter','latex');
ylabel('$R(z)$','Interpreter','latex');
set(gca,'FontName','SansSerif','FontSize',16);

%% Fig 4 with inset

open('../Images/PNAS Figures/Fig4/AM_Collapse.fig');

axes('position',[3*10^(-2) 3*10^(-2) 0.45 0.45])
box on;
semilogx(z,R,'k','LineWidth',2); hold on;
for i = 2:df_Nx
semilogx(df_z,df_R{i},'Color',myColors(i),'LineWidth',1.25);
end
xlabel('$z$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
ylabel('$R(z/)$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);

%%

close all;

h1 = openfig('../Images/PNAS Figures/Fig4/RSS_Collapse_Updated10.4.23.fig');
ax1 = gca;
h2 = openfig('../Images/PNAS Figures/Fig4/AM_Collapse_Updated10.4.23.fig');
ax2 = gca;
% h3 = openfig('../Images/PNAS Figures/Fig4/AM_NoCollapse.fig');
% ax3 = gca;
% 
% h4 = openfig('../Images/PNAS Figures/Fig4/CombinedCollapse.fig');

h4 = figure();
tcl = tiledlayout(1,2);
ax1.Parent=tcl;
ax1.Layout.Tile=1;
% ax1.Layout.TileSpan=[1,2];
ax2.Parent=tcl;
ax2.Layout.Tile=2;
% ax2.Layout.TileSpan=[1,2];
% ax3.Parent=tcl;
% ax3.Layout.Tile=5;



%% FIGURE 2 IN MANUSCRIPT

close all;

lim1 = surf_xyz(1,1)-1850;
lim2 = surf_xyz(end,1)-1850;

delta_ibl = [0.006701,0.109,0.1943,0.1943,0.1602,0.1602,0.2355,0.2856,...
    0.4198,0.5089].*delta;

figure();

plot(surf_xyz(:,1)-1850,surf_xyz(:,end)./delta,'k','LineWidth',2); hold on;
fp1 = plot(f,'b');
fp1.LineWidth = 2;
for i = 1:9
    plot(x_hat(i),delta_ibl(i)./delta,'kdiamond','MarkerFaceColor',...
        myColors(i),'MarkerSize',10);
end
plot(-100,0.003,'kdiamond','MarkerFace','white','MarkerSize',10);
fp1 = plot(f,'b');
fp1.LineWidth = 2.5;
ylabel('$z/\delta_{ABL}$','FontSize',16,'FontName','SansSerif',...
    'interpreter','latex');
% for i = 1:df_N/2
%     myText = strcat(num2str(round(x(i))),' [m]');
%     lt = xline(x(i),'--k',{myText});
% end
xlabel('$\hat{x} = x - x_0 [m]$','FontSize',16,'FontName','SansSerif',...
    'interpreter','latex');
% text(-1000,200,sprintf('IBL Height = ax^b \n\na = %f\nb = %f',f.a,f.b));
legend('off')
xlim([lim1 lim2]);
ax = gca;
ax.FontSize = 16;

%%

% close all;
% 
% df_T = df_t.*(u_bar(end)/delta);
% iT = find(df_T >= 30,1);
% 
% [DF_T,DF_Z] = meshgrid(df_T(1:iT),df_z);
% DF_T = DF_T';
% DF_Z = DF_Z';
% 
% %map = [1 0 0
% %    1 1 1
%  %   0 0 1];
% 
% 
% figure();
% %contourf(DF_T,DF_Z,df_upplus{1}(1:iT,:));
% contourf(DF_T,DF_Z./delta,df_upplus{10}(1:iT,:),'LineColor','none');
% xlabel('t*U_\infty/\delta');
% ylabel('z/\delta');
% title(strcat('x = ',num2str(round(x(10),0)),'m'));
% xlim([0 30]);
% colormap(jet);
% colorbar;
% % set(gca,'yscale','log')


%%

% open('../Images/DF500_MedV4/Misc Plots/Time Height Contours/x1900.fig');
% 
% open('../Images/DF500_MedV4/Misc Plots/Time Height Contours/x2586.fig');
% 
% open('../Images/DF500_MedV4/Misc Plots/Time Height Contours/x3518.fig');
% 
% open('../Images/DF500_MedV4/Misc Plots/Time Height Contours/x4788.fig');
% 
% open('../Images/DF500_MedV4/Misc Plots/Time Height Contours/x7601.fig');
%% EXIT