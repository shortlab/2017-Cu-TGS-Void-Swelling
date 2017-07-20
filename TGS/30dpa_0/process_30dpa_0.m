%Process the 30 dpa copper SC data that was irradiated at Sandia

grat=4.7823; %Calibrated grating spacing in um +/- 3.4e-9m from 3 batches
grat_err=3.4e-9;
grat2=4.7781; %Calibrated grating spacing in um +/- 5.2e-9 from 3 batches
grat2_err=5.2e-9;
angles=0:5:60;
SAW_freqs=zeros(3,length(angles));
PSAW_freqs=zeros(3,length(angles));
SAW_err=zeros(3,length(angles));
PSAW_err=zeros(3,length(angles));

temp_SAW_freq=zeros(3,1);
temp_PSAW_freq=zeros(3,1);

temp_SAW_err=zeros(3,1);
temp_PSAW_err=zeros(3,1);

str_base='Cu_Sandia_2016-30dpa_0-04.80um-';

for ii=1:length(angles)
    display(angles(ii))
    for jj=1:3
        if angles(ii)<10
            angle_str=strcat('0',num2str(angles(ii)));
        else
            angle_str=num2str(angles(ii));
        end
        [peaks,~,errors]=param_extract_time(jj,strcat(str_base,angle_str,'_deg-POS-',num2str(jj),'.txt'),strcat(str_base,angle_str,'_deg-NEG-',num2str(jj),'.txt'),grat,1);
        temp_SAW_freq(jj)=peaks(2);
        temp_PSAW_freq(jj)=peaks(1);
        temp_SAW_err(jj)=errors(1,:,2);
        temp_PSAW_err(jj)=errors(1,:,1);
    end
    SAW_freqs(:,ii)=temp_SAW_freq;
    PSAW_freqs(:,ii)=temp_PSAW_freq;
    
    SAW_err(:,ii)=temp_SAW_err;
    PSAW_err(:,ii)=temp_PSAW_err;
end

SAW_freq_means=zeros(1,length(angles));
SAW_freq_std=zeros(1,length(angles));

PSAW_freq_means=zeros(1,length(angles));
PSAW_freq_std=zeros(1,length(angles));

for kk=1:length(SAW_freq_means)
    SAW_freq_means(kk)=mean(SAW_freqs(:,kk));
    SAW_freq_std(kk)=std(SAW_freqs(:,kk));
    
    PSAW_freq_means(kk)=mean(PSAW_freqs(:,kk));
    PSAW_freq_std(kk)=std(PSAW_freqs(:,kk));
end

SAW_speeds=SAW_freq_means*grat*10^(-6);
PSAW_speeds=PSAW_freq_means*grat*10^(-6);

SAW_speed_err=sqrt((grat*10^(-6)*SAW_freq_std).^2+(SAW_freq_means*grat_err).^2);
PSAW_speed_err=sqrt((grat*10^(-6)*PSAW_freq_std).^2+(PSAW_freq_means*grat_err).^2);

figure()
errorbar(angles,SAW_speeds,SAW_speed_err);

% figure()
% plot(angles,SAW_speeds,'b')
% hold on
% plot(angles,PSAW_speeds,'r')