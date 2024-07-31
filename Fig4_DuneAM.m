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
rho = 1.23; % Density of air 
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

%% Dune Wall Stress Data

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

%% Load data for RSS

myDir = dir('./DuneField/DuneData/RSSData/x*');

Ndir = length(myDir);

iuw = 0;


for i = 1:Ndir
    

    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    myPath = strcat(myFolder,'/',myName);
    
    data = load(myPath);
    
    if contains(myPath,'comp(u_rey,1)')
        iuw = iuw + 1;
        uw{iuw} = data(:,4:end);
    end
        
end

clear my* Ndir i* data


%% Time-averaged RSS


N_uw = length(uw);

for i = 1:N_uw
    
    tempuw = uw{i};
    UW{i} = mean(tempuw);
    
end
    
clear temp*


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





%% Plot RSS

close all;

set(0,'defaultTextInterpreter','latex');

myColors = [ "Black", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000", "#000000" ];

delta_ibl = [0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta;

z_plot = linspace(1.5,400,100);
uw_plot = cell2mat(UW');

figure();
tiledlayout(1,2);
p1 = nexttile;
plot(uw_plot(1,:),z_plot,'k--','LineWidth',2); hold on;
for i = 2:10
        plot(uw_plot(i,:),z_plot,'Color',myColors(i),'LineWidth',2); hold on;
end
set(gca,'FontName','SansSerif','FontSize',20);
ylabel('$z$','FontName','SansSerif','FontSize',36);
xlabel('$\langle u^\prime w^\prime \rangle$','FontName','SansSerif','FontSize',36);
ylim([0 200]);
legend({'Alkali Flat','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$',...
    '$\hat{x}_4$','$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$',...
    '$\hat{x}_8$','$\hat{x}_9$'},'Interpreter','Latex',...
    'Location','NorthWest','NumColumns',3,'FontSize',30);


uw_plot2 = uw_plot(2:end,:);

p2 = nexttile;
% figure();
plot(uw_plot(1,:)./(0.12^2),z_plot./30,'k--','LineWidth',2); hold on;
for i = 2:10
    if i == 2
        plot(uw_plot2(i,:)./(0.12^2),z_plot./30,'Color',myColors(i),'LineWidth',2); hold on;
    else
        plot(uw_plot2(i,:)./(0.12^2),z_plot./delta_ibl(i-1),'Color',myColors(i),'LineWidth',2); hold on;
    end
end
set(gca,'FontName','SansSerif','FontSize',20);
ylabel('$z/\hat{\delta}$','FontName','SansSerif','FontSize',36);
xlabel('$\langle u^\prime w^\prime \rangle/u^2_{\tau,0}$','FontName','SansSerif','FontSize',36);
ylim([0 4]);


%% FIG 4B and 4C

close all;

% 4C Inset
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
        semilogx(df_z./30,df_R{i},'Color',myColors(i),'LineWidth',1.25);
    else
        semilogx(df_z./delta_ibl(i-1),df_R{i},'Color',myColors(i),'LineWidth',1.25); 
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




%% END

