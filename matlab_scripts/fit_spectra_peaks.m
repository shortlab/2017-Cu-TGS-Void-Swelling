function [peak,errs,ft]=fit_spectra_peaks(fft,plotty)
%Feed function a power spectral density and it will return the position of
%the peak and the error in the gaussian fitting

[mx_val,ind]=max(fft(:,2));

min_pct_fit=0.50;

pos_ind=ind;
go_p=1;
while go_p
    pos_ind=pos_ind+1;
    if fft(pos_ind,2)<=min_pct_fit*mx_val
        go_p=0;
    end
end

neg_ind=ind;
go_n=1;
while go_n
    neg_ind=neg_ind-1;
    if fft(neg_ind,2)<=min_pct_fit*mx_val
        go_n=0;
    end
end

pk_trace=fft(neg_ind:pos_ind,:);

% This is the gaussian fitting for the peaks, which does well enough
ft=fit(pk_trace(:,1),pk_trace(:,2),'gauss1');
peak=ft.b1;
error_mat=confint(ft,0.95);
errs=[peak-error_mat(1,2); error_mat(2,2)-peak];

if plotty
    figure()
    plot(pk_trace(:,1),ft(pk_trace(:,1)),pk_trace(:,1),pk_trace(:,2));
end

%Below is a record of some previous peak fitting procedures that I tried at
%some point.

% Write a custom fit, or get one, for a lorentzian
% LB=[0 0 0 0];
% UB=[10^-6 10^9 10^9 10^-6];
% ST=[10^-8 510*10^6 10^7 10^-7];
% OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Start',ST,'Upper',UB);
% TYPE=fittype('A./(1+((x-x0)/(w/2)).^2)+B;','options',OPS,'coefficients',{'A','x0','w','B'});
% [f0 gof]=fit(pk_trace(:,1),pk_trace(:,2),TYPE);

%Fit a scaled version of the peak so the fitting algorithm doesn't have to
%handle 16 orders of magnitude
% trace_scaled=[pk_trace(:,1)/(10^8) pk_trace(:,2)/100];
% LB=[0 4.5 0 0];
% UB=[10 6 10 10];
% ST=[1 5.15 .1 1];
% OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Start',ST,'Upper',UB);
% TYPE=fittype('A./(1+((x-x0)/(w/2)).^2)+B;','options',OPS,'coefficients',{'A','x0','w','B'});
% [f0 gof]=fit(trace_scaled(:,1),trace_scaled(:,2),TYPE);
% 
% peak=(f0.x0)*10^8;
% 
% display(f0);

% %Contains the one sigma confidence interval on the fitted parameters. 
% errs=confint(f0,0.95);
% pk_err=((errs(2,2)-errs(1,2))/2)*10^8;

% display(errs(1,2));
% display(errs(2,2));
% display(pk_err);

% fft_scaled=[fft_1(:,1)/(10^8) fft_1(:,2)/100];
% 
% figure()
% plot(trace_scaled(:,1),f0(trace_scaled(:,1)),'r-');
% hold on
% plot(fft_scaled(:,1),fft_scaled(:,2),'g-');

% Try a center of mass weighting method to find the peak position
% peak=(sum(pk_trace(:,1).*pk_trace(:,2)))/(sum(pk_trace(:,2)));
end