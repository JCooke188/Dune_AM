%% START

% Created by Justin Cooke
% Purpose of the script is to reproduce Fig. 3 in 
% https://www.pnas.org/doi/10.1073/pnas.2320216121

%% START

clc;
clear;
close all;

set(0,'defaultTextInterpreter','latex');

%% Load Data

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

%% Load Spanwise Slice Data  

myDir = dir('./DuneField/WSS/Slice/SpanwiseWSS*');

surf_N = length(myDir);

cx = 0;
cr = 0;

swss_x = zeros();

for i = 1:surf_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        swss_xyz = load(loadMe);
        cx = cx + 1;
        swss_x(cx) = round(swss_xyz(1,2));
        swss_z{cx} = swss_xyz(:,end); % added just now
    else
        swss_temp = load(loadMe);
        cr = cr + 1;
        swss_tau{cr} = swss_temp(end-149:end,4:end);
    end
    
    
end


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


%% Get Average dune height k_a at each spanwise slice

k_a = zeros(40,1);
k_rms = zeros(40,1);

for i = 1:40
   
    temp = swss_z{i};
    k_a(i) = mean(temp);
    k_rms(i) = rms(temp-k_a(i));
    
end

clear temp

%% Now Plot Everything

close all;

% Data points in the Alkali Flat
flat_wss_x = [715 1100 1385 1615]; 
flat_wss = [0.02098, 0.02075, 0.02093, 0.02089];
smooth_wss = mean(flat_wss);

% Add moving average
swss_tau_movMean = movmean((swss_tau_da./smooth_wss),5);


figure();
yyaxis left
plot(surf_xyz(:,1)-1850,surf_tau_bar./smooth_wss,...
    'Color',"#CA1F00",'LineWidth',1.25); hold on;
plot(swss_x-1850,swss_tau_da./smooth_wss,...
    'k^','MarkerFaceColor',"#0037CA",'MarkerSize',7.5);
plot(swss_x-1850,swss_tau_movMean,'-','LineWidth',1.5,'Color',...
    'black');
ylabel('$\tau_{b}/\tau_0$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
ax1 = gca;
ax1.YColor = "#CA1F00";
yyaxis right
area(surf_xyz(:,1)-1850,surf_xyz(:,end)./300,'FaceColor',"#808080",...
    'FaceAlpha',0.2,'EdgeColor',"#808080"); hold on;
errorbar(swss_x-1850,k_a./300,k_rms./300,'o','MarkerSize',7.5,...
    'Color','black','MarkerFaceColor','#808080');
ax2 = gca;
ax2.YColor = "#808080";
ylim([0 0.25]);
ylabel('$z/\delta_{ABL}$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
xlabel('$\hat{x} = x - x_0 [m]$','Interpreter','latex',...
    'FontName','SansSerif','FontSize',16);
xlim([-1846 6300]);
ax1.FontSize = 16;




%% END
