dir_base='/Users/Cody/Documents/MIT/Short_Lab/LSAW/Measurements/080_Data/LSAW/Tungsten_Calibration/';

day1='2015-11-03/';
day2='2015-11-05/';

freq=zeros(1,length(1:50));
therm=zeros(1,length(1:50));
freq_fit_error=zeros(2,length(1:50));

go_here=strcat(dir_base,day1);
cd(go_here)

base1='Tungsten_Calibration-2015-11-03-06.40um-spot';

for j=1:30
    if j<10
        pos_str=strcat(base1,'0',num2str(j),'-POS-1.txt');
        neg_str=strcat(base1,'0',num2str(j),'-NEG-1.txt');
    else
        pos_str=strcat(base1,num2str(j),'-POS-1.txt');
        neg_str=strcat(base1,num2str(j),'-NEG-1.txt');
    end
    [freq(j),therm(j),freq_fit_error(:,j)]=param_extract(3,pos_str,neg_str,4.80);
end

go_here=strcat(dir_base,day2);
cd(go_here)

base2='Tungsten_Calibration-2015-11-05-04.80um-spot';

for j=31:50
    pos_str=strcat(base2,num2str(j),'-POS-1.txt');
    neg_str=strcat(base2,num2str(j),'-NEG-1.txt');
    [freq(j),therm(j),freq_fit_error(:,j)]=param_extract(3,pos_str,neg_str,4.80);
end

freq_ave=mean(freq);
freq_std=std(freq);

figure()
plot(1:50,freq,'k*');