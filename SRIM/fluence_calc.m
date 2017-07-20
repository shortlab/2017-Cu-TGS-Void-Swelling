%Calculate the necessary fluences to get a series of displacements for a
%given vacancy.txt file from SRIM

target_peak_doses=[0.5 1 5 10 15 25 50]; %for Ni WRONG
%target_peak_doses=[5 10 30 50 75 100]; %for Cu

%Make this same damn things do a peak damage calc for me based on actual
%fluences measured
% fluence_meas[]*10^16;

Ni=0;
Cu=1;

energy=35; %in MeV

if Ni
    cd(strcat('Ni_Ni/',num2str(energy),'MeV/'))
    n=9.14*10^28*10^-6; %in cm^-3
elseif Cu

    
    cd(strcat('Cu_Cu/',num2str(energy),'MeV/'))
    n=8.491*10^28*10^-6; %in cm^-3
else
    display('you done fucked up');
end

vac_data=dlmread('VACANCY.txt','',[27 0 126 2]);
depth=vac_data(:,1)/(10^4); %depth in um
total_vac=vac_data(:,2)+vac_data(:,3); %in #/angstrom*ion
total_vac_cm=total_vac*10^8; %in #/cm*ion

[peak_vac,peak_location]=max(total_vac_cm);

fluences=(target_peak_doses*n)/peak_vac; %in #/cm^2

display(fluences)
display(depth(peak_location))

figure('Position',[100 100 150 360])
% plot(depth,total_vac_cm*fluences(end)/n,'r-')
plot(total_vac_cm*(6.3184*10^16)/n,depth,'r-','linewidth',1.5)
hold on
plot([0 95],[4.8 4.8],'k--');
xlim([0 95])
ylim([0 5.25])
set(gca,...
'Ydir','reverse',...
'Xdir','reverse',...
'xaxisLocation','top',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,...
'FontName','Times')
ylabel({'Depth ($\mu$m)'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',16,...
'FontName','Times')
xlabel({'Dose (dpa)'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',16,...
'FontName','Times')
% title({strcat(num2str(energy),' MeV Cu$^{6+}$ into Cu')},...
% 'FontUnits','points',...
% 'interpreter','latex',...
% 'FontSize',20,...
% 'FontName','Times');


cd('../..')

    