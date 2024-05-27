%% START

% Created by Justin Cooke
% Purpose of the script is to reproduce Fig. 2 in 
% https://www.pnas.org/doi/10.1073/pnas.2320216121

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Data

% Local Directory
myDir = dir('./DuneField/Combined/x_*');

% On GitHub
% myDir = dir('./DuneField/DuneData/x_*');

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

% Local and GitHub have same Directory names
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

% Local and GitHub have same Directory names
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

%% Params we need

nu = 0.32e-4; %[[m^2/s]
delta_abl = 300; % [m]
U_infty = 6.0758; % [m/s] found at highest probe point

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
    
    df_Retau(i) = delta_abl/df_dnu(i);

end

% Spanwise station based averaging

for i = 1:df_N/2
   
    swss_taubar(i) = mean(swss_tauw{i},'All');
        
end

df_taubar = mean(surf_tau);

%% Dune specific data

[df_Nt,df_Nz] = size(df_uxdata{1});
df_Nx = df_N/2;

myLL = log(1.5);
myUL = log(400);

thisLogZ = linspace(myLL,myUL,30);
df_z = exp(thisLogZ);

for i = 1:df_Nx

    df_ubar{i} = mean(df_uxdata{i});
    df_uprime{i} = df_uxdata{i} - df_ubar{i};
    
    df_upplus{i} = df_uprime{i}./df_utau(i);
    
end

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

%% Determine IBL Location

% Need to determine where the value of the difference being calculated goes
% to zero. Threshold is ~10^-4 for where values begin to flatten out.

close all;

figure();

xhatdelta = [1.22,2.45,3.89,5.56,7.52,9.79,12.45,15.55,19.17];

for i = 1:df_Nx-2
    
    nexttile;
    semilogx(df_z./delta_abl,ratio{i},'ko','MarkerFace','red'); hold on;
    xlabel('$z/\delta_{ABL}$','FontSize',16,'FontName','SanSerif');
    ylabel('$a(z)$','FontSize',16,'FontName','SanSerif');
       
    set(gca,'FontSize',16,'FontName','SanSerif');
end

%% Find a Fit

% Found by inspection
delta_ibl = [0,0.109,0.1943,0.1943,0.1602,0.1602,0.2355,0.2856,...
   0.4198,0.5089].*delta_abl;

% Want a fit of the type y = a*x^b
ft = fittype('a*x^b');
f = fit(x_hat',delta_ibl'./delta_abl,ft,'StartPoint',[0,0]);

%% Now Plot

myColors = [ "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000", "	#000000" ];

close all;

lim1 = surf_xyz(1,1)-1850;
lim2 = surf_xyz(end,1)-1850;

figure();

plot(surf_xyz(:,1)-1850,surf_xyz(:,end)./delta_abl,'k','LineWidth',2); hold on;
fp1 = plot(f,'b');
fp1.LineWidth = 2;
for i = 1:9
    plot(x_hat(i),delta_ibl(i)./delta_abl,'kdiamond','MarkerFaceColor',...
        myColors(i),'MarkerSize',10);
end
plot(-100,0.003,'kdiamond','MarkerFace','white','MarkerSize',10);
fp1 = plot(f,'b');
fp1.LineWidth = 2.5;
asl_line = yline(0.1); %Added a line to signify delta_asl
asl_line.LineWidth = 2.5;
asl_line.LineStyle = '--';
ylabel('$z/\delta_{ABL}$','FontSize',16,'FontName','SansSerif');
xlabel('$\hat{x} = x - x_0 [m]$','FontSize',16,'FontName','SansSerif');
legend('off')
xlim([lim1 lim2]);
ax = gca;
ax.FontSize = 16;
