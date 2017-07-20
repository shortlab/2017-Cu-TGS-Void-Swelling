function [out] = make_fft(num_in,pos_file,neg_file,mod,strt,pass)
%   NOTE: File command to edit tab autocomplete file is as follows
%       edit(fullfile(matlabroot,'toolbox/local/TC.xml'))

nelson=1;

derivative=0;
pulse=0;
grat=5.50; %in um
plotty=1;
plotfft=1;
saveout=1;

custom=0;
therm=0;
thermsaw=0;

if nargin==3
    %     mod='qexp';
    mod='therm';
    strt=5;
end

if nargin==4
    strt=5;
end

if nargin<6
    pass=cell(1,num_in);
end


pos_base=pos_file(1:end-5);
neg_base=neg_file(1:end-5);

if nelson
    
    % Old Nelson has different header length
    % pos=dlmread(strcat(pos_base,num2str(num_in),'.txt'),'',21,0);
    % neg=dlmread(strcat(neg_base,num2str(num_in),'.txt'),'',21,0);
    
    hdr_len=22;
else
    hdr_len=15;
end

pos=dlmread(strcat(pos_base,num2str(num_in),'.txt'),'',hdr_len,0);
neg=dlmread(strcat(neg_base,num2str(num_in),'.txt'),'',hdr_len,0);

fixed=[pos(:,1) pos(:,2)-neg(:,2)];

pass{num_in}=fixed;

if num_in>1
    out=make_fft(num_in-1,pos_file,neg_file,mod,strt,pass);
else
    x_tot=pass{1}(:,1);
    y_tot=pass{1}(:,2);
    
    for j=2:length(pass)
        x_tot=x_tot+pass{j}(:,1);
        y_tot=y_tot+pass{j}(:,2);
    end
    
    dat=[x_tot/length(pass) y_tot/length(pass)];
    
    %     %For data saved directily off the scope on 4/28/15 only
    %     dat=dlmread('Al_SC_111.csv',',',8,0);
    %     flr=mean(dat(end-500:end,2));
    %     dat(:,2)=dat(:,2)-flr;
    
    %     %For Nb data saved directly from the scope only
    %     pos=dlmread('C3pes_oft00000.csv',',',8,0);
    %     neg=dlmread('C3pes_oft00001.csv',',',8,0);
    %     dat=[pos(:,1) pos(:,2)-neg(:,2)];
    
    if saveout
        dlmwrite('dat_temp.txt',dat);
    end
    
    
    if plotty
        figure()
        plot(dat(:,1),dat(:,2),'r')
    end
    
    [mx ind]=max(dat(:,2));
    
         out=dat;
    
    fixed_short=dat(ind:end,:);
    fixed_short(:,2)=fixed_short(:,2)-mean(fixed_short(end-50:end,2));
    %     fixed_short(:,2)=fixed_short(:,2)-mean(dat(1:10,2));
    
    
    %     out=fixed_short;    %This saves temporal profile out, not fft
    %     display(fixed_short(:,1))
    %     figure()
    %     plot(fixed_short(:,1),fixed_short(:,2))
    
    custom=0;
    therm=0;
    
    if strcmp(mod,'exp1') %Single Exponential
        mod_str=mod;
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
    elseif strcmp(mod,'exp2') %Sum of Exponentials
        mod_str=mod;
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
    elseif strcmp(mod,'power2') % Two term power fit
        mod_str=mod;
    elseif strcmp(mod,'qexp') %Forced two term decaying exponential fit
        custom=1;
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
    elseif strcmp(mod, 'therm') %fit thermal decay to find diffusivity
        therm=1;
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
        %         fixed_short(:,1)=fixed_short(:,1)-195*(fixed_short(2,1)-fixed_short(1,1));
    elseif strcmp(mod, 'thermsaw')
        thermsaw=1;
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
    else
        mod_str='power1'; %single term power fit (goes something like 1/x)
    end
    
    %     out=fixed_short;
    
    %Fit data to one of a couple of models determined by the mod argument in
    %the function
    
    if custom
        q=2*pi/(grat*10^(-6));
        %         LB=[0 0 0 0];
        %         UB=[.7 .7 1000 10000];
        %         ST=[0.5 0.5 100 300];
        LB=[0 0];
        UB=[1 1];
        ST=[0.05 10^-5];
        OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
        %         TYPE=fittype('A.*exp(-x./(t1*(10^-9)))+B.*exp(-x./(t2*(10^-9)));','options',OPS,'coefficients',{'A','B','t1','t2'});
        TYPE=fittype('A.*x^(-1/2)*exp(-x*k*q^2)','options',OPS,'problem','q','coefficients',{'A','k'});
        [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',q);
        %         display(f0)
        out=f0.k;
    elseif therm
        q=2*pi/(grat*10^(-6));
        LB=[0 0 0];
        UB=[10 5*10^-4 0.1];
        ST=[0.05 1*10^-4 0];
        OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
        TYPE=fittype('A.*erfc(q*sqrt(k*x))+c;','options',OPS,'problem','q','coefficients',{'A','k','c'});
        [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',q);
        
        display(f0.k);
        
        
        %         [f0 gof]=fit(fixed_short(strt:end,1),sgolayfilt(fixed_short(strt:end,2),1,101),TYPE,'problem',q);
        %         out=f0.k;
    elseif thermsaw
        q=2*pi/(grat*10^(-6));
        LB=[0 0 0 1.027*10^9 0 0 0];
        UB=[10 10^-3 10 1.028*10^9 2*pi 10^-6 1];
        ST=[0.05 1*10^-4 0.01 1.026*10^9 0 10^-5 0];
        OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
        TYPE=fittype('A.*erfc(q*sqrt(k*x))+B.*sin(2*pi*f*x+p)*exp(-x/t)*c;','options',OPS,'problem','q','coefficients',{'A','k','B','f','p','t','c'});
        
        [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',q);
        
        display(f0.k);
        %         display(f0.A);
        %         display(f0.B);
        
    else
        f0=fit(fixed_short(:,1),fixed_short(:,2),mod_str);
    end
    
    %display(f0); %Show fit parameters
    
    if plotty
        
        %plot fit to see how well it matches
        figure()
        plot(fixed_short(strt:end,1),f0(fixed_short(strt:end,1)),'b',fixed_short(strt:end,1),fixed_short(strt:end,2),'r');
        %     hold on
        %     plot(fixed_short(strt:end,1),sgolayfilt(fixed_short(strt:end,2),1,25),'g');
        %     plot(out(:,1),out(:,2),'r')
        %     out=fixed_short;
    end
    
    flat=[fixed_short(strt:end,1) fixed_short(strt:end,2)-f0(fixed_short(strt:end,1))];
    %     short=[fixed_short(strt:end,1) fixed_short(strt:end,2)];
    
    if pulse
        %Use make_pulse to zero signal and mirror it about t=0 to create a
        %pulse shaped thing. This allows for the use of a filtering window in
        %the periodogram. Also cut off tails of the distribution that just look
        %like noise
        flat=make_pulse(flat(1:end,:));
        
        % Find the stuff we'll need to take the spectral profile
        num=length(flat(:,1));
        fs=num/(flat(end,1)-flat(1,1));
        p=14; %magnitude of zero padding (Was 15 from Jeff)
        pdsize=2^p-num-2; %more padding = smoother transform
        
        %Pad equally on both ends.
        pad_val=mean(flat(end-50:end,2));
        pad=zeros(pdsize,1);
        pad(1:end)=pad_val;
        tstep=flat(end,1)-flat(end-1,1);
        tpad=flat(end,1):tstep:flat(end,1)+(pdsize-1)*tstep;
        
        flat_pad=[-flipud(tpad') pad;flat(:,1) flat(:,2);tpad' pad];
        
        %         flat_pad=padarray(flat,pdsize,mean(flat(end-50:end,2)),'post');
        %         flat_pad=padarray(flat_pad,pdsize,mean(flat_pad(1:50,2)),'pre');
    else
        if derivative
            tstep=flat(end,1)-flat(end-1,1);
            d_flat=diff(flat(:,2))/tstep;
            flat=[flat(1:length(d_flat),1) d_flat];
        end
        
        % Find the stuff we'll need to take the spectral profile
        num=length(flat(:,1));
        fs=num/(flat(end,1)-flat(1,1));
        p=18; %magnitude of zero padding (Was 15 from Jeff)
        pdsize=2^p-num-2; %more padding = smoother transform
        
        %Only pad on the positive end
        pad_val=mean(flat(end-50:end,2));
        pad=zeros(pdsize,1);
        pad(1:end)=pad_val;
        tstep=flat(end,1)-flat(end-1,1);
        tpad=flat(end,1):tstep:flat(end,1)+(pdsize-1)*tstep;
        
        flat_pad=[flat(:,1) flat(:,2);tpad' pad];
        %         short_pad=[short(:,1) short(:,2);tpad' pad];
        
        %         flat_pad=padarray(flat,pdsize,mean(flat(end-50:end,2)),'post'); %pad with DC level of signal near end
    end
    
    if plotty
        figure()
        plot(flat(:,1),flat(:,2))
        %     plot(flat_pad(:,1),flat_pad(:,2));
    end
    
    
    nfft=length(flat_pad(:,2));
    
    %Find the Power Spectral density
    [psd freq]=periodogram(flat_pad(:,2),rectwin(nfft),nfft,fs); %periodogram method
    %     [psd freq]=periodogram(flat_pad(:,2),hamming(nfft,'periodic'),nfft,fs);
    
    %Use a hamming window and a Welchs method. Hamming does the best of the
    %ones I've tried and Welch does slightly better than the normal
    %periodogram.
    %     [psd freq]=pwelch(flat_pad(:,2),hamming(nfft,'periodic'),[],nfft,fs);
    %     [psd_dc freq_dc]=pwelch(short_pad(:,2),hamming(nfft,'periodic'),[],nfft,fs);
    
    %Don't save out DC spike in FFT
    amp=sqrt(psd(5:end));
    out=[freq(5:end) amp];
    
    if saveout
        dlmwrite('dat_spec.txt',out);
    end
    %out=f0;
    
    if plotfft
        figure()
        hold on
        plot(freq(5:end),amp,'r');
        xlim([0 1.7e9]);
        %xlim([0 1e11]);
        title(num2str(p));
    end
    
    %     figure()
    %     inv=ifft(sqrt(psd_dc-psd));
    %     norm_tim=1:length(inv);
    %     back_trans=[];
    %     back_trans(:,1)=5e-11*norm_tim;
    %     back_trans(:,2)=inv;
    %     plot(back_trans(:,1),back_trans(:,2))
    %
    %     saw_off=fixed_short(:,2)-inv(1:length(fixed_short(:,2)));
    %     figure()
    %     plot(fixed_short(:,1),saw_off(:))
    
    %     [max_fq ind_fq]=max(out(:,2));
    %     sz=150;
    %     %ind_fq=ind_fq-35;
    %     freq_sm=out(ind_fq-sz:ind_fq+sz,:);
    %     [ltz_fit ltz_par ltz_resnorm ltz_resid]=lorentzfit(freq_sm(:,1),freq_sm(:,2));
    %     hold on
    %     plot(freq_sm(:,1),ltz_fit,'r-');
    %     display(ltz_par(2));
    %
end

end