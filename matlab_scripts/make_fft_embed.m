function [fft,diffusivity] = make_fft_embed(num_in,pos_file,neg_file,mod,strt,grat,peak,pass)
%   NOTE: File command to edit tab autocomplete file is as follows
%       edit(fullfile(matlabroot,'toolbox/local/TC.xml'))

nelson=0;

derivative=1;
plotty=0;
plotfft=0;
saveout=0;

if nargin==3
    mod='therm';
    strt=5;
end

if nargin==4
    strt=6;
end

if nargin<8
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
    [fft,diffusivity]=make_fft_embed(num_in-1,pos_file,neg_file,mod,strt,grat,peak,pass);
else
    x_tot=pass{1}(:,1);
    y_tot=pass{1}(:,2);
    
    for j=2:length(pass)
        x_tot=x_tot+pass{j}(:,1);
        y_tot=y_tot+pass{j}(:,2);
    end
    
    dat=[x_tot/length(pass) y_tot/length(pass)];
    
    if saveout
        dlmwrite('dat_temp.txt',dat);
    end
    
    if plotty
        figure()
        plot(dat(:,1),dat(:,2),'r')
    end
    
    [mx ind]=max(dat(:,2));
    
    fixed_short=dat(ind:end,:);
    fixed_short(:,2)=fixed_short(:,2)-mean(fixed_short(end-50:end,2));
    
    custom=0;
    therm=0;
    thermsaw=0;
    
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
    elseif strcmp(mod, 'thermsaw')
        thermsaw=1;
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
    elseif strcmp(mod, 'thermsaw2')
        if length(peak)<=1
            thermsaw=1;
        else
            thermsaw2=1;
        end
        fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
    else
        mod_str='power1'; %single term power fit (goes something like 1/x)
    end
    
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
        diffusivity=f0.k;
    elseif therm
        q=2*pi/(grat*10^(-6));
        LB=[0 0 0];
        UB=[10 5*10^-4 0.1];
        ST=[0.05 1*10^-4 0];
        OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
        TYPE=fittype('A.*erfc(q*sqrt(k*x))+c;','options',OPS,'problem','q','coefficients',{'A','k','c'});
        [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',q);
        diffusivity=f0.k;
    elseif thermsaw
        q=2*pi/(grat*10^(-6));
        LB=[0 0 0 0 0 0];
        UB=[10 10^-3 10 2*pi 10^-6 1];
        ST=[0.05 1*10^-4 0.01 0 10^-5 0];
        OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
        TYPE=fittype('A.*erfc(q*sqrt(k*x))+B.*sin(2*pi*f*x+p)*exp(-x/t)+C;','options',OPS,'problem',{'q','f'},'coefficients',{'A','k','B','p','t','C'});
        
        [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',{q,peak});
        diffusivity=f0.k;
    elseif thermsaw2
        q=2*pi/(grat*10^(-6));
        LB=[0 0 0 0 0 0 0 0 0];
        UB=[10 10^-3 10 2*pi 10^-6 10 2*pi 10^-6 1];
        ST=[0.05 1*10^-4 0.01 0 10^-5 0.01 0 10^-5 0];
        OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
        TYPE=fittype('A.*erfc(q*sqrt(k*x))+B.*sin(2*pi*f1*x+p1)*exp(-x/t1)+C.*sin(2*pi*f2*x+p2)*exp(-x/t2)+D;','options',OPS,'problem',{'q','f1','f2'},'coefficients',{'A','k','B','p1','t1','C','p2','t2','D'});
        
        [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',{q,peak(1),peak(2)});
        diffusivity=f0.k;
        %display(f0)
    else
        f0=fit(fixed_short(:,1),fixed_short(:,2),mod_str);
        diffusivity=0;
    end
    
    if plotty
        
        %plot fit to see how well it matches
        figure()
        plot(fixed_short(strt:end,1),f0(fixed_short(strt:end,1)),'b',fixed_short(strt:end,1),fixed_short(strt:end,2),'r');
        
    end
    
    flat=[fixed_short(strt:end,1) fixed_short(strt:end,2)-f0(fixed_short(strt:end,1))];
    
    if plotty
        figure()
        plot(flat(:,1),flat(:,2))
    end
    
    % Time step info necessary for differentiation and flat padding
    tstep=flat(end,1)-flat(end-1,1);
    
    % If option selected, take transform of derivative of recorded signal
    % to filter out DC junk
    if derivative
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
    tpad=flat(end,1):tstep:flat(end,1)+(pdsize-1)*tstep;
    
    flat_pad=[flat(:,1) flat(:,2);tpad' pad];
    
    nfft=length(flat_pad(:,2));
    
    %Find the Power Spectral density
    
    %Use a hamming window and a Welchs method. Hamming does the best of the
    %ones I've tried and Welch does slightly better than the normal
    %periodogram.
%     [psd freq]=pwelch(flat_pad(:,2),hamming(nfft,'periodic'),[],nfft,fs);
    [psd freq]=periodogram(flat_pad(:,2),rectwin(nfft),nfft,fs); %periodogram method
    
    %Don't save out DC spike in FFT
    amp=sqrt(psd(5:end));
    fft=[freq(5:end) amp];
    
    if saveout
        dlmwrite('dat_spec.txt',out);
    end
    
    if plotfft
        figure()
        hold on
        plot(freq(5:end),amp,'r');
        xlim([0 1.7e9]);
    end
    
end

end