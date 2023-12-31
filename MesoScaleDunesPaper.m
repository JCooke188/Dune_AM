% Created by Justin Cooke
% Purpose of this script is to create a centralized script containing all
% of the code used to create the figues and analyze the data in Cooke et
% al. 2024 "Mesoscale Structure of the Atmospheric Boundary Layer Across a
% Natural Roughness Transition".

clc;
clear;
close all;

%% Load the Data -- Original Data

% Loading data for Alkali Flat, this part is !slow! because there is a
% lot of data and it's being run serially in MATLAB. Sorry, I didn't know
% how to use python at time of writing.

myDir = dir('./AlkaliFlat/OrigData/z*');

N = length(myDir);

cr = 0;
cx = 0;

for i = 1:N
   
    tempFile = myDir(i).name;
    tempFolder = myDir.folder;
    
    tempName = strcat(tempFolder,'/',tempFile);
    
    if contains(tempName,'README')
        xyz = load(tempName); % Load all xyz information
        % Only keep x = 1000m, 1250m, 1500m, 1750m
        xyz = [xyz(301,2:end),xyz(468,2:end),xyz(635,2:end),xyz(801,2:end)];
        cr = cr + 1;
        loc_data{cr} = xyz;
    elseif contains(tempName,'comp')
        u_x = load(tempName); % Load all streamwise velocity
        % Same as above
        u_x = [u_x(:,301),u_x(:,468),u_x(:,635),u_x(:,801)];
        cx = cx + 1;
        ux_data{cx} = u_x;
    end
    
    
end

clear temp* myDir cx cr urms_x xyz u_x

%% Now load in additional Alkali Flat data

% This data was only probed at four x locations (x = 1000m, 1250m, 1500m,
% and 1750m) hence why the previous data which was probed at many more was
% truncated to only these four locations. 

myDir = dir('./AlkaliFlat/Data/z*comp*');

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

% Can choose to run the Alkali Flat analysis at any of the four x-stations,
% in the paper we chose x = 1750 due to its proximity to the dune field
% (~100m upstream).

store_ux_data = ux_data;

count = 0;
for i = 1:N
    % Ignore these data points because they are repeats due to probe 
    % spacing being smaller than grid spacing
    if i ~= 2 && i ~= 3 && i ~= 5 
        count = count + 1;
        ux_data_trunc{count} = store_ux_data{i}(:,4); % Change store_ux_data{i}(:,X) for x = 1,2,3,4
    end
end

%% Necessary information needed for creating cutoff filter plus constants 

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

u_tau = 0.1313; % Value at 1805
nu = 0.32e-4; % Kin. viscosity of air
deltanu = nu/u_tau; 
delta = 300; % ABL Height
zplus = z./deltanu;
zdelta = z./delta;
delta_asl = 30;

% Because we convert u_bar from cell to matrix and if you re-run this
% and don't do this, you'll have a bad time
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
        % Using Taylor's Frozen Hypothesis to go from t -> x
        lam_x{i,j} = u_bar{i}(j)./omega;
        lam_xp{i,j} = (u_bar{i}(j)/u_tau)./(omega.*(deltanu/u_tau));
        kx{i,j} = (2*pi)./lam_x{i,j};
    end
    
    omega_p{i} = omega.*(deltanu/u_tau);
    
end

u_bar = cell2mat(u_bar);

t = linspace(0,Lt,Nt);
T = t.*(u_bar(end)/delta); % Large-eddy turnover time
total_T = T(end);


%% Create Large Scales

% Lambda_x cutoff chosen as lambda_x,c = delta_abl
lx_cutoff = 300; 

half = Nt/2;

for i = 1:Nz
    
    % Convert to omega using Taylor's Frozen Hypothesis u_bar changes with
    % z
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

%% Now do Hilbert Transform
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

%% Plot Sanity Check

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

%% Check point

% The above was just for the alkali flat, now we have the data from the
% dune field, which was captured at 10 log-spaced x location within the
% dune field. Overall flow is the same as above. Now we've added wall
% stress data though, and we plan to calculate internal boundary layer
% heights

%% Load Dune Data

myDir = dir('./DuneField/DuneData/x_*');

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
 
% Streamwise wall stress averaging in time

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

% Spanwise station based averaging so now we have time and space average

for i = 1:df_N/2
   
    swss_taubar(i) = mean(swss_tauw{i},'All');
        
end

df_taubar = mean(surf_tau);

close all;

x_hat = x - 1850;

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

% Same procedure as before and same cutoff as before.

df_half = df_Nt/2;

for i = 1:df_Nx

    df_utz = df_upplus{i}; % u'(z,t) at single x-loc
    df_ubarz = df_ubar{i}; % u_bar(z) at single x-loc
    df_lx_cutoff = lx_cutoff; 
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

% Didn't use this plot, but shows R_AM at each individual x-location 

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



%% Create <u'u'>

for i = 1:df_Nx
    
    for j = 1:df_Nz
   
        temp = df_upplus{i}(:,j);
        store_dfuu{j} = temp.*temp;
        
    end
    
    mat_dfuu = cell2mat(store_dfuu);
    df_upup{i} = mat_dfuu;
    df_upupbar{i} = mean(df_upup{i});
    
end



%% Calculate IBL Height using the method of Li et al.
% Ref: M Li, et al., Experimental study of a turbulent boundary layer with
% a rough-to-smooth change in surface conditions at high Reynolds numbers.
% J. Fluid Mech. 923, A18 (2021)

% Highest avg velocity data at last station
U_infty = df_ubar{end}(end);

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

%% Determine IBL Location

% Need to determine where the value of the difference being calculated goes
% to zero. Threshold is ~10^-4 for where values begin to flatten out.

close all;

figure();

xhatdelta = [1.22,2.45,3.89,5.56,7.52,9.79,12.45,15.55,19.17];

for i = 1:2:df_Nx-2
    
    nexttile;
    semilogx(df_z./delta,ratio{i},'ko','MarkerFace','red'); hold on;
    xlabel('z/\delta_{ABL}','Interpreter','tex','FontSize',16,'FontName','SanSerif');
    ylabel('a(z)','Interpreter','tex','FontSize',16,'FontName','SanSerif');
       
    set(gca,'FontSize',16,'FontName','SanSerif');
end

%% Now plot IBL heights and find the data fit 


delta_ibl = [0,0.109,0.1943,0.1943,0.1602,0.1602,0.2355,0.2856,...
    0.4198,0.5089].*delta;

% Fit type -> delta_ibl/z02 = a*(x/z02)^b
ft = fittype('a*x^b');
f = fit(x_hat',delta_ibl'./delta,ft,'StartPoint',[0,0]);
 

figure();
plot(x_hat,delta_ibl./delta,'k^','MarkerFace','blue'); hold on;
plot(f,'b');
text(-1750,200,sprintf('Li Method:\na = %f\nb = %f',f.a,f.b),'Color','blue',...
    'FontName','SansSerif','FontSize',14);

xlabel('x - x_0 [m]','FontName','SansSerif','FontSize',16);
ylabel('\delta_{ibl}','FontName','SansSerif','FontSize',16);
legend('Li \it{et al} Method',...
    'location','northwest','Box','off','NumColumns',2,'Orientation','Horizontal');
set(gca,'FontSize',16,'FontName','SansSerif');
yticks([0 100 200 300]);
xticks([-2000 0 2000 4000 6000]);
xticklabels({'-2000','0','2000','4000','6000'});
xlim([-2000 6000]);

%% Single-Point Dune AM Normalized with delta_ibl

% Initial plot but opted to create two separate plots then stitched
% together.

close all;


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
ax1.FontSize = 16;

ax2 = subplot(1,2,2);
for i = 2:df_Nx
semilogx(df_z/delta_ibl(i),df_R{i},'Color',myColors(i),'LineWidth',1.25); hold on;
end
semilogx(z/30,R,'k-.','LineWidth',1.5);
% title('Collapse of R_{AM} with \delta_{ibl} in Dune Field');
yline(0,'--k');
ylim([-0.8,0.5]);
xlim([5*10^-3 10^1]);
xlabel('$z/\delta_{IBL}$','Interpreter','latex');
ylabel('$R(z/\delta_{IBL})$','Interpreter','latex');

ax2.FontSize = 16;

%% FIG 4B and 4C

close all;

% 4B
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

% 4C
figure();
semilogx(z/30,R,'k--','LineWidth',2); hold on;
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
set(gca,'FontName','SansSerif','FontSize',16);


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
xlabel('$\hat{x} = x - x_0 [m]$','FontSize',16,'FontName','SansSerif',...
    'interpreter','latex');
legend('off')
xlim([lim1 lim2]);
ax = gca;
ax.FontSize = 16;

%% Figure 3 in Manuscript

% First load the data

myDir = dir('./DuneField/WSS/Slice/SpanwiseWSS*');

surf_N = length(myDir);

cx = 0;
cr = 0;

swss_x = zeros();

tic;

for i = 1:surf_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        swss_xyz = load(loadMe);
        cx = cx + 1;
        swss_x(cx) = round(swss_xyz(1,2));
    else
        swss_temp = load(loadMe);
        cr = cr + 1;
        swss_tau{cr} = swss_temp(end-149:end,4:end);
    end
    
    
end

toc;

clear myName myFolder loadMe cr cx swss_xyz swss_temp


%% Average the Data

% Average the full field data in time only
surf_tau_bar = mean(surf_tau);

% Average the spanwise data in time and space
swss_Nx = surf_N/2;

for i = 1:swss_Nx
   
    temp = swss_tau{i};
    swss_tau_bar{i} = mean(temp); % time average
    swss_tau_da(i) = mean(swss_tau_bar{i}); % average across y
    
end


%% Create Figure 3

close all;

flat_wss_x = [715 1100 1385 1615];
flat_wss = [0.02098, 0.02075, 0.02093, 0.02089];
smooth_wss = mean(flat_wss);

figure();
yyaxis left
plot(surf_xyz(:,1)-1850,surf_tau_bar./smooth_wss,...
    'Color',"#CA1F00",'LineWidth',1.25); hold on;
plot(swss_x-1850,swss_tau_da./smooth_wss,...
    'k^','MarkerFaceColor',"#0037CA",'MarkerSize',7.5);
%plot(swss_x(1:2:end)-1850,y_fit,'--','LineWidth',1.5,'Color','#0037CA');
ylabel('$\tau_{b}/\tau_0$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
ax1 = gca;
ax1.YColor = "#CA1F00";
yyaxis right
area(surf_xyz(:,1)-1850,surf_xyz(:,end)./300,'FaceColor',"#808080",...
    'FaceAlpha',0.2,'EdgeColor',"#808080");
ax2 = gca;
ax2.YColor = "#808080";
ylim([0 0.075]);
ylabel('$z/\delta_{ABL}$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
xlabel('$\hat{x} = x - x_0$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
xlim([-1846 6300]);
ax1.FontSize = 16;



%% EXIT