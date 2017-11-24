%Data dump of Irradiated {111} SC copper

%define fit for max finding
ST=[1850 50 0];
OPS=fitoptions('Method','NonLinearLeastSquares','Start',ST);
TYPE=fittype('A+B*cos(x*(pi/30)+p)','options',OPS,'coefficients',{'A','B','p'});
%Need different parameters as a starting point for min finding
ST2=[1750 50 0];
OPS2=fitoptions('Method','NonLinearLeastSquares','Start',ST2);
TYPE2=fittype('A+B*sin(x*(pi/30)+p)','options',OPS2,'coefficients',{'A','B','p'});

%0dpa data from PRB paper
% % SAW_angles_0dpa=0:5:60;
% % SAW_speeds_0dpa=[1.7346 1.8170 1.9033 1.9402 1.9009 1.8117 1.7284 1.6534 1.6095 1.5941 1.6126 1.6640 1.7408]*10^3;
% %
% % f00=fit(SAW_angles_0dpa(1:floor(end/2))',SAW_speeds_0dpa(1:1:floor(end/2))','fourier1');
% % f00_neg=fit(SAW_angles_0dpa(floor(end/2):end)',SAW_speeds_0dpa(floor(end/2):end)','fourier1');
% %
% % PSAW_angles_0dpa=[0 5 25 30 35 40 45 50 55 60];
% % PSAW_speeds_0dpa=[2.5391 2.4274 2.4208 2.5064 2.6082 2.6737 2.6980 2.6719 2.6070 2.5183]*10^3;

%0dpa data taken on new dual-detector setup
SAW_angles_0dpa_0_base=0:5:60;
SAW_angle_shift_0dpa_0=6.51; %from 7
SAW_speeds_0dpa_0=[1.8607 1.9283 1.9444 1.8897 1.8046 1.7173 1.6479 1.6096 1.6018 1.6287 1.6838 1.7652 1.8536]*10^3;

SAW_angles_0dpa_0_base=SAW_angles_0dpa_0_base(1:end-1);
SAW_speeds_0dpa_0=SAW_speeds_0dpa_0(1:end-1);

[SAW_angles_0dpa_0,SAW_speeds_0dpa_0]=phase_shift(SAW_angles_0dpa_0_base,SAW_speeds_0dpa_0,60,SAW_angle_shift_0dpa_0);

f00=fit(SAW_angles_0dpa_0(1:floor(end/2))',SAW_speeds_0dpa_0(1:1:floor(end/2))',TYPE);
f00_neg=fit(SAW_angles_0dpa_0(floor(end/2):end)',SAW_speeds_0dpa_0(floor(end/2):end)',TYPE2);

PSAW_speeds_0dpa_0=[2.4209 2.4072 2.4847 2.5814 2.6701 2.7128 2.7153 2.6588 2.6081 2.4884 2.3968]*10^3;
PSAW_angles_0dpa_0=[0 15 20 25 30 35 40 45 50 55 60];

%error in freqs
SAW_freq_err_0dpa_0=[1.1658 1.1658 1.7625 2.0190 1.1658 2.0194 1.1656 1.1657 0.4407 2.4533 1.1657 0.7633 1.7623]*10^5;
PSAW_freq_err_0dpa_0=[3.6792 1.2709 1.6123 1.5154 1.7172 1.7251 2.2587 2.2565 4.9612 4.8700 4.2300]*10^6; %yes, 10^6

%error in SAW w/o grat error
% % SAW_speed_err_0dpa_0=[0.6363 0.6363 0.9619 1.1020 0.6363 1.1021 0.6362 0.6363 0.2405 1.3390 0.6363 0.4166 0.9618];
% % PSAW_speed_err_0dpa_0=[20.0809 6.9363 8.7996 8.2711 9.3720 9.4154 12.3277 12.3158 27.0778 26.5799 23.0869];

%error in SAW speed w/ grat error from 2 meas - real
SAW_speed_err_0dpa_0=[5.8305 6.0399 6.1323 5.9883 5.6569 5.4613 5.1720 5.0536 4.9949 5.2468 5.2832 5.5139 5.8532];
PSAW_speed_err_0dpa_0=[21.4500 10.2141 11.7186 11.5352 12.5302 12.6509 14.9499 14.8412 28.2701 27.6869 24.2639];

%error in SAW speed w/ grat error representative of 3 - more like others,
%not real
% % SAW_speed_err_0dpa_0=[0.8607 0.8750 1.1367 1.2493 0.8490 1.2251 0.8174 0.8100 0.5539 1.4319 0.8246 0.6898 1.1218];
% % PSAW_speed_err_0dpa_0=[20.0950 6.9767 8.8335 8.3101 9.4089 9.4533 12.3567 12.3436 27.0900 26.5912 23.0989];

SAW_speed_err_0dpa_0=SAW_speed_err_0dpa_0(1:end-1);

[~,SAW_speed_err_0dpa_0]=phase_shift(SAW_angles_0dpa_0_base,SAW_speed_err_0dpa_0,60,SAW_angle_shift_0dpa_0);

%92dpa data from Sandia [first data processing run, no peak sliding]
SAW_angles_92dpa_0=0:5:60;
SAW_speeds_92dpa_0=[1.7064 1.7777 1.8377 1.8600 1.8229 1.7523 1.6964 1.6362 1.5970 1.5901 1.6075 1.6549 1.7207]*10^3;

PSAW_angles_92dpa_0=30:5:60;
PSAW_speeds_92dpa_0=[2.2838 2.4896 2.5610 2.6108 2.5897 2.4957 2.3629]*10^3;

%92dpa data from Sandia, the set to use

SAW_angles_92dpa_base=0:5:60;
SAW_angle_shift_92dpa=0.76; %from 1
SAW_speeds_92dpa=[1.7098  1.7785 1.8395 1.8615 1.8243 1.7556 1.6983 1.6382 1.5989 1.5903 1.6098 1.6538 1.7204]*10^3;
SAW_angles_92dpa_base=SAW_angles_92dpa_base(1:end-1);
SAW_speeds_92dpa=SAW_speeds_92dpa(1:end-1);

[SAW_angles_92dpa,SAW_speeds_92dpa]=phase_shift(SAW_angles_92dpa_base,SAW_speeds_92dpa,60,SAW_angle_shift_92dpa);

f92=fit(SAW_angles_92dpa(1:floor(end/2))',SAW_speeds_92dpa(1:1:floor(end/2))',TYPE);
f92_neg=fit(SAW_angles_92dpa(floor(end/2):end)',SAW_speeds_92dpa(floor(end/2):end)',TYPE2);

PSAW_speeds_92dpa=[2.4524 2.5115 2.5817 2.6201 2.5920 2.5250 2.4158]*10^3;
PSAW_angles_92dpa=[30 35 40 45 50 55 60];

%error in freqs
SAW_freq_err_92dpa=[0.8811 1.7622 3.9161 4.4059 4.5999 6.6676 3.0843 0.4410 0.8812 3.0842 4.6630 2.6801 2.3315]*10^5;
PSAW_freq_err_92dpa=[0.3093 0.4704 0.0764 0.1152 0.3024 0.2683 1.0255]*10^7;

%error in SAW w/o grat error
% % SAW_speed_err_92dpa=[0.4196 0.8392 1.8650 2.0983 2.1907 3.1755 1.4689 0.2100 0.4197 1.4689 2.2208 1.2764 1.1104];
% % PSAW_speed_err_92dpa=[14.7290 22.4050 3.6405 5.4880 14.4024 12.7794 48.8388];

%error in SAW speed w/ grat error from 3 meas
SAW_speed_err_92dpa=[1.9485 2.1498 2.7693 2.9486 2.9867 3.7283 2.3937 1.8351 1.8282 2.3000 2.8533 2.2397 2.2133];
PSAW_speed_err_92dpa=[14.9798 22.5787 4.6377 6.2145 14.6884 13.0847 48.9127];

SAW_speed_err_92dpa=SAW_speed_err_92dpa(1:end-1);

[~,SAW_speed_err_92dpa]=phase_shift(SAW_angles_92dpa_base,SAW_speed_err_92dpa,60,SAW_angle_shift_92dpa);

%50dpa data from Sandia

SAW_angles_50dpa_base=0:5:60;
SAW_angle_shift_50dpa=22.95; %from 22.5
SAW_speeds_50dpa=[1.8348 1.7501 1.6821 1.6218 1.5962 1.6024 1.6343 1.6961 1.7726 1.8471 1.9064 1.9113 1.8586]*10^3;
SAW_angles_50dpa_base=SAW_angles_50dpa_base(1:end-1);
SAW_speeds_50dpa=SAW_speeds_50dpa(1:end-1);

[SAW_angles_50dpa,SAW_speeds_50dpa]=phase_shift(SAW_angles_50dpa_base,SAW_speeds_50dpa,60,SAW_angle_shift_50dpa);

f50=fit(SAW_angles_50dpa(1:floor(end/2))',SAW_speeds_50dpa(1:1:floor(end/2))',TYPE);
f50_neg=fit(SAW_angles_50dpa(floor(end/2):end)',SAW_speeds_50dpa(floor(end/2):end)',TYPE2);

PSAW_speeds_50dpa=[2.5020 2.6109 2.6394 2.6227 2.5575]*10^3; %Could maybe pick up more of these with the derivative on to clean up the DC end
PSAW_angles_50dpa=[15 20 25 30 35];

%error in freqs
SAW_freq_err_50dpa=[3.9163 8.8126 9.5620 6.4002 8.4063 6.4906 5.9600 3.6062 5.2872 1.5256 6.6675 4.0381 1.9206]*10^5;
PSAW_freq_err_50dpa=[2.6831 0.2442 0.2653 0.1946 1.5143]*10^7;

%error in SAW w/o grat error
% % SAW_speed_err_50dpa=[1.8729 4.2144 4.5728 3.0608 4.0202 3.1040 2.8503 1.7246 2.5285 0.7296 3.1886 1.9311 0.9185];
% % PSAW_speed_err_50dpa=[128.3159 11.6787 12.6898 9.3072 72.4184];

%error in SAW speed w/ grat error from 3 meas
SAW_speed_err_50dpa=[2.2824 4.3943 4.7266 3.2707 4.1773 3.3065 3.0780 2.1043 2.8252 1.5022 3.4647 2.3613 1.6093];
PSAW_speed_err_50dpa=[128.3282 11.8253 12.8278 9.4921 72.4412];

SAW_speed_err_50dpa=SAW_speed_err_50dpa(1:end-1);

[~,SAW_speed_err_50dpa]=phase_shift(SAW_angles_50dpa_base,SAW_speed_err_50dpa,60,SAW_angle_shift_50dpa);

%30dpa data from Sandia

SAW_angles_30dpa_base=0:5:60;
SAW_angle_shift_30dpa=11.3; %from 11.5
SAW_speeds_30dpa=[1.9036 1.9213 1.8703 1.7911 1.7104 1.6456 1.6073 1.6002 1.6244 1.6788 1.7525 1.8345 1.9028]*10^3;
SAW_angles_30dpa_base=SAW_angles_30dpa_base(1:end-1);
SAW_speeds_30dpa=SAW_speeds_30dpa(1:end-1);

[SAW_angles_30dpa,SAW_speeds_30dpa]=phase_shift(SAW_angles_30dpa_base,SAW_speeds_30dpa,60,SAW_angle_shift_30dpa);

f30=fit(SAW_angles_30dpa(1:floor(end/2))',SAW_speeds_30dpa(1:1:floor(end/2))',TYPE);
f30_neg=fit(SAW_angles_30dpa(floor(end/2):end)',SAW_speeds_30dpa(floor(end/2):end)',TYPE2);

PSAW_speeds_30dpa=[2.3464 2.4320 2.4747 2.5884 2.6545 2.6725 2.6323 2.5711 2.4856 2.3951 2.3069]*10^3;
PSAW_angles_30dpa=[10 15 20 25 30 35 40 45 50 55 60];

%error in freqs
SAW_freq_err_30dpa=[0.7633 1.9206 2.4530 1.1657 1.1657 0.7631 0.6026 0.8813 1.5264 0.7633 1.3220 1.1657 1.5886]*10^5; %note that value in arrary index 7 here came out chance way too low. Took error for this point (0.6) from the SAW_err from peak fitting in process script
PSAW_freq_err_30dpa=[1.6899 1.8337 1.9160 2.5383 2.0904 1.4666 1.0962 3.5377 2.8882 2.7264 2.4599]*10^6;

%error in SAW w/o grat error
% % SAW_speed_err_30dpa=[0.3650 0.9185 1.1731 0.5575 0.5575 0.3650 0.2882 0.4214 0.7300 0.3650 0.6322 0.5575 0.7597];
% % PSAW_speed_err_30dpa=[8.0815 8.7691 9.1628 12.1391 9.9971 7.0137 5.2423 16.9186 13.8124 13.0384 11.7639];

%error in SAW speed w/ grat error from 3 meas
SAW_speed_err_30dpa=[1.4017 1.6460 1.7732 1.3901 1.3377 1.2255 1.1785 1.2132 1.3663 1.2481 1.3972 1.4184 1.5515];
PSAW_speed_err_30dpa=[8.2519 8.9380 9.3302 12.2778 10.1737 7.2665 5.5664 17.0170 13.9250 13.1491 11.8777];

SAW_speed_err_30dpa=SAW_speed_err_30dpa(1:end-1);

[~,SAW_speed_err_30dpa]=phase_shift(SAW_angles_30dpa_base,SAW_speed_err_30dpa,60,SAW_angle_shift_30dpa);

%10dpa data from Sandia, first process through used derivative filter, but
%found a small systematic difference, so I turned it off for the second
%time I ran it, those number reported here

%This is the data from the first run with the derivative filter turned on
% SAW_speeds_10dpa=[1.6605 1.7327 1.8200 1.8943 1.9248 1.8769 1.8015 1.7180 1.6541 1.6163 1.5918 1.6148 1.6614]*10^3;
% PSAW_angles_10dpa=0:5:60; % Haven't down selected the real ones
% PSAW_speeds_10dpa=[2.5206 2.4202 2.3106 6.7401 8.9692 8.9687 2.4238 2.5551 2.6748 2.7489 2.7644 2.7482 2.6950]*10^3;

SAW_angles_10dpa_base=0:5:60;
SAW_angle_shift_10dpa=55.64; %from 56
SAW_speeds_10dpa=[1.6762 1.7419 1.8298 1.9051 1.9325 1.8839 1.8099 1.7225 1.6593 1.6163 1.6017 1.6195 1.6719]*10^3;
SAW_angles_10dpa_base=SAW_angles_10dpa_base(1:end-1);
SAW_speeds_10dpa=SAW_speeds_10dpa(1:end-1);

[SAW_angles_10dpa,SAW_speeds_10dpa]=phase_shift(SAW_angles_10dpa_base,SAW_speeds_10dpa,60,SAW_angle_shift_10dpa);

f10=fit(SAW_angles_10dpa(1:floor(end/2))',SAW_speeds_10dpa(1:1:floor(end/2))',TYPE);
f10_neg=fit(SAW_angles_10dpa(floor(end/2):end)',SAW_speeds_10dpa(floor(end/2):end)',TYPE2);

PSAW_speeds_10dpa=[2.5360 2.4571 2.3168 2.5558 2.6786 2.7461 2.7628 2.7416 2.6734]*10^3;
PSAW_angles_10dpa=[0 5 10 35 40 45 50 55 60];

%error in freqs
SAW_freq_err_10dpa=[0.2203 0.8088 0.6882 0.5728 0.5194 1.1252 0.5342 0.5082 0.1588 0.2330 0.4965 0.7644 1.2407]*10^6;
PSAW_freq_err_10dpa=[0.0919 0.1164 0.1346 0.0452 0.0448 0.0458 0.0131 0.0226 0.0457]*10^8;

%error in SAW w/o grat error
% % SAW_speed_err_10dpa=[1.0536 3.8683 3.2914 2.7396 2.4842 5.3814 2.5549 2.4305 0.7595 1.1145 2.3746 3.6559 5.9340];
% % PSAW_speed_err_10dpa=[43.9711 55.6676 64.3907 21.6062 21.4036 21.9016 6.2829 10.7921 21.8694];

%error in SAW speed w/ grat error from 3 meas
SAW_speed_err_10dpa=[1.1479 3.8972 3.3288 2.7882 2.5391 5.4057 2.6019 2.4752 0.8834 1.1979 2.4142 3.6823 5.9514];
PSAW_speed_err_10dpa=[43.9765 55.6716 64.3938 21.6174 21.4160 21.9144 6.3277 10.8178 21.8815];

SAW_speed_err_10dpa=SAW_speed_err_10dpa(1:end-1);

[~,SAW_speed_err_10dpa]=phase_shift(SAW_angles_10dpa_base,SAW_speed_err_10dpa,60,SAW_angle_shift_10dpa);

%5dpa data from Sandia, don't need the filted on this one becuase it came
%out pretty nice

SAW_angles_05dpa=0:5:60;
SAW_angle_shift_05dpa=15;
SAW_speeds_05dpa=[1.9634 1.9505 1.8778 1.7934 1.7149 1.6469 1.6136 1.6131 1.6458 1.7107 1.7867 1.8813 1.9364]*10^3;

[SAW_angles_05dpa,SAW_speeds_05dpa]=phase_shift(SAW_angles_05dpa,SAW_speeds_05dpa,60,SAW_angle_shift_05dpa);

PSAW_angles_05dpa=0:5:60;
PSAW_speeds_05dpa=[9.0211 8.9747 2.4751 2.5617 2.6399 2.6965 2.6966 2.6376 2.5929 2.5196 2.4470 2.3360 8.6425]*10^3;

%Second try of 5 dpa data, something got messed up with the angles on the
%first one as the 0 and 60 wound up not lining up very well.

SAW_angles_05dpa_1_base=0:5:60;
SAW_angle_shift_05dpa_1=8.48;
SAW_speeds_05dpa_1=[1.9231 1.9690 1.9466 1.8762 1.7886 1.7039 1.6430 1.6154 1.6203 1.6550 1.7220 1.8077 1.8998]*10^3;
SAW_angles_05dpa_1_base=SAW_angles_05dpa_1_base(1:end-1);
SAW_speeds_05dpa_1=SAW_speeds_05dpa_1(1:end-1);

[SAW_angles_05dpa_1,SAW_speeds_05dpa_1]=phase_shift(SAW_angles_05dpa_1_base,SAW_speeds_05dpa_1,60,SAW_angle_shift_05dpa_1);

f05=fit(SAW_angles_05dpa_1(1:floor(end/2))',SAW_speeds_05dpa_1(1:1:floor(end/2))',TYPE);
f05_neg=fit(SAW_angles_05dpa_1(floor(end/2):end)',SAW_speeds_05dpa_1(floor(end/2):end)',TYPE2);

PSAW_speeds_05dpa_1=[2.4265 2.4155 2.4819 2.5770 2.6584 2.7028 2.6957 2.6209 2.5641 2.5043 2.4324 2.3193]*10^3;
PSAW_angles_05dpa_1=[0 10 15 20 25 30 35 40 45 50 55 60];

%error in freqs
SAW_freq_err_05dpa_1=[2.4531 0.7633 0.4405 2.2893 0.7630 2.0191 2.3314 2.6803 1.3219 2.6436 2.6799 1.7623 0.4404]*10^5;
PSAW_freq_err_05dpa_1=[3.2649 3.8943 0.9562 0.9231 1.9847 0.4766 0.7333 3.4038 9.5837 7.9507 1.1929 3.9221]*10^6;

%error in SAW w/o grat error
% % SAW_speed_err_05dpa_1=[1.1736 0.3652 0.2108 1.0952 0.3650 0.9659 1.1154 1.2823 0.6324 1.2647 1.2821 0.8431 0.2107];
% % PSAW_speed_err_05dpa_1=[15.6195 18.6306 4.5746 4.4163 9.4948 2.2800 3.5082 16.2840 45.8495 38.0370 5.7068 18.7636];

%error in SAW speed w/ grat error from 3 meas
SAW_speed_err_05dpa_1=[4.5748 4.5419 4.4807 4.4508 4.1287 4.0351 3.9390 3.9294 3.7788 4.0099 4.1618 4.2411 4.3733];
PSAW_speed_err_05dpa_1=[16.5860 7.3139 7.3900 11.2921 6.6195 7.1222 17.3632 46.2270 38.4703 7.9904 19.5067];

SAW_speed_err_05dpa_1=SAW_speed_err_05dpa_1(1:end-1);

[~,SAW_speed_err_05dpa_1]=phase_shift(SAW_angles_05dpa_1_base,SAW_speed_err_05dpa_1,60,SAW_angle_shift_05dpa_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use errors at the closest points to the shifted line as representative of
%the error in the percent deviation

max_00dpa=max(f00(0:.01:30));
max_05dpa=max(f05(0:.01:30));
max_10dpa=max(f10(0:.01:30));
max_30dpa=max(f30(0:.01:30));
max_50dpa=max(f50(0:.01:30));
max_92dpa=max(f92(0:.01:30));

max_vec=[max_00dpa max_05dpa max_10dpa max_30dpa max_50dpa max_92dpa];

max_sigma_00dpa=SAW_speed_err_0dpa_0(4);
max_sigma_05dpa=SAW_speed_err_05dpa_1(3);
max_sigma_10dpa=SAW_speed_err_10dpa(4);
max_sigma_30dpa=SAW_speed_err_30dpa(4);
max_sigma_50dpa=(SAW_speed_err_50dpa(4)+SAW_speed_err_50dpa(3))/2;
max_sigma_92dpa=SAW_speed_err_92dpa(4);

max_sigma_vec=[max_sigma_00dpa max_sigma_05dpa max_sigma_10dpa max_sigma_30dpa max_sigma_50dpa max_sigma_92dpa];

% max_error_vec=sqrt((max_sigma_vec/max_vec(1)).^2+(max_vec*max_sigma_vec(1)/max_vec(1)^2).^2);
max_error_vec=max_sigma_vec/max_vec(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting below here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!!!!
% Want to check on the shift angle by comparing the curves of the fitted
% functions, it doesn't matter since I'm using the computed max from each
% amd I'm pretty close, but better to be better on the plot I'm showing
% !!!!!

plot_errs=0;

figure('Position', [100, 100, 1000, 500])
plot([15 15],[1550 2000],'k--');
hold on
plot([45 45],[1550 2000],'k--');
hold on
p0=plot(SAW_angles_0dpa_0,SAW_speeds_0dpa_0,'ko','MarkerSize',6,'MarkerFaceColor','k','DisplayName','0 dpa');
hold on
if plot_errs
    errorbar(SAW_angles_0dpa_0,SAW_speeds_0dpa_0,SAW_speed_err_0dpa_0,'k-');
else
    plot(SAW_angles_0dpa_0,SAW_speeds_0dpa_0,'k-');
end
hold on
p5=plot(SAW_angles_05dpa_1,SAW_speeds_05dpa_1,'d','color',[1 0.5 0],'MarkerSize',6,'MarkerFaceColor',[1 0.5 0],'DisplayName','5 dpa');
hold on
if plot_errs
    errorbar(SAW_angles_05dpa_1,SAW_speeds_05dpa_1,SAW_speed_err_05dpa_1,'-','color',[1 0.5 0]);
else
    plot(SAW_angles_05dpa_1,SAW_speeds_05dpa_1,'-','color',[1 0.5 0]);
end
hold on
p10=plot(SAW_angles_10dpa,SAW_speeds_10dpa,'mp','MarkerSize',8,'MarkerFaceColor','m','DisplayName','10 dpa');
hold on
if plot_errs
    errorbar(SAW_angles_10dpa,SAW_speeds_10dpa,SAW_speed_err_10dpa,'m-');
else
    plot(SAW_angles_10dpa,SAW_speeds_10dpa,'m-');
end
hold on
p30=plot(SAW_angles_30dpa,SAW_speeds_30dpa,'b^','MarkerSize',6,'MarkerFaceColor','b','DisplayName','30 dpa');
hold on
if plot_errs
    errorbar(SAW_angles_30dpa,SAW_speeds_30dpa,SAW_speed_err_30dpa,'b-');
else
    plot(SAW_angles_30dpa,SAW_speeds_30dpa,'b-');
end
hold on
p50=plot(SAW_angles_50dpa,SAW_speeds_50dpa,'s','color',[0 0.5 0],'MarkerSize',6,'MarkerFaceColor',[0 0.5 0],'DisplayName','50 dpa');
hold on
if plot_errs
    errorbar(SAW_angles_50dpa,SAW_speeds_50dpa,SAW_speed_err_50dpa,'-','color',[0 0.5 0]);
else
    plot(SAW_angles_50dpa,SAW_speeds_50dpa,'-','color',[0 0.5 0]);
end
hold on
p90=plot(SAW_angles_92dpa,SAW_speeds_92dpa,'rv','MarkerSize',6,'MarkerFaceColor','r','DisplayName','90 dpa');
hold on
if plot_errs
    errorbar(SAW_angles_92dpa,SAW_speeds_92dpa,SAW_speed_err_92dpa,'r-');
else
    plot(SAW_angles_92dpa,SAW_speeds_92dpa,'r-'); 
end
hold on
legend([p0 p5 p10 p30 p50 p90])
% legend('Baseline','5 DPA','10 DPA','30 DPA','50 DPA','92 DPA');
% legend('Baseline','Baseline','10 DPA','10 DPA','30 DPA','30 DPA','50 DPA','50 DPA','92 DPA','92 DPA');

set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'LineWidth',1.25)
ylabel({'SAW speed (m/s)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')
xlabel({'Relative surface angle (degrees)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')

figure()
errorbar([0 5 10 30 50 90],max_vec/max_vec(1),max_error_vec,'k--','LineWidth',1.25);
hold on
% plot([0 5 10 30 50 90],max_vec/max_vec(1),'k--');
% hold on
plot([0 5 10 30 50 90],max_vec/max_vec(1),'ko','MarkerSize',6,'MarkerFaceColor','k');
xlim([0 100]);
set(gca,...
    'FontUnits','points',...
    'FontSize',16,...
    'FontName','Times',...
    'LineWidth',1.25)
ylabel({'Speed along $\langle 11\bar{2}\rangle$($c/c_0$)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Times')
xlabel({'Peak dose (dpa)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Times')

% figure()
% errorbar([0 5 10 30 50 90],max_vec/max_vec(1),max_error_vec);

figure()
plot([0 5 10 30 50 90],[min(f00_neg(30:.01:60)) min(f05_neg(30:0.01:60)) min(f10_neg(30:0.01:60)) min(f30_neg(30:.01:60)) min(f50_neg(30:.01:60)) min(f92_neg(30:.01:60))]/min(f00_neg(30:.01:60)),'k-');

display(max_vec(end))
display(min(f92_neg(30:.01:60)))

anisotropy=max_vec./[min(f00_neg(30:.01:60)) min(f05_neg(30:0.01:60)) min(f10_neg(30:0.01:60)) min(f30_neg(30:.01:60)) min(f50_neg(30:.01:60)) min(f92_neg(30:.01:60))];

% sample=0:0.01:30;
% figure()
% plot(sample,f00(sample),sample,f10(sample),sample,f30(sample),sample,f50(sample),sample,f92(sample));
