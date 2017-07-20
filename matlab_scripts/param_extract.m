function [freq,therm,errs]=param_extract(num,pos_file,neg_file,grat,two_peak)

peak_fit_plot=0;

if nargin<5
    two_peak=0;
end

[fft_1,therm_1]=make_fft_embed(num,pos_file,neg_file,'therm',6,grat,0);

[peak1,errs1,ft_1]=fit_spectra_peaks(fft_1,peak_fit_plot);

%Include case of meaningfully double-peaked spectra. This is the case where
%both a SAW and a PSAW are detectable
if two_peak
    minus_peak_1=[fft_1(:,1) fft_1(:,2)-ft_1(fft_1(:,1))];
%     figure()
%     plot(fft_1(:,1),fft_1(:,2),'b-');
%     hold on
%     plot(minus_peak_1(:,1),minus_peak_1(:,2),'r-');
    [peak2,errs2,ft_2]=fit_spectra_peaks(minus_peak_1,peak_fit_plot);
    
    freq=[peak1 peak2];
    if freq(1)>freq(2)
        errs(:,:,1)=errs1;
        errs(:,:,2)=errs2;
    else
        freq=[freq(2) freq(1)];
        errs(:,:,1)=errs2;
        errs(:,:,2)=errs1;
    end
    
else
    freq=peak1;
    errs=errs1;
end

if two_peak
    [fft_2,therm_2]=make_fft_embed(num,pos_file,neg_file,'thermsaw2',6,grat,[peak1 peak2]);
else
    [fft_2,therm_2]=make_fft_embed(num,pos_file,neg_file,'thermsaw',6,grat,peak1);
end

therm=therm_2;

end