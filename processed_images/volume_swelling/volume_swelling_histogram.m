%
run void_data.m

do_both=1;

thickness_05=0.170; %in um or 230 or 105nm? +/- 20nm
thickness_10=0.200; %+/- 5nm
thickness_30=0.145; %+/- 2nm
thickness_50=0.220; %+/- 7nm
thickness_90=0.200; %+/- 7 nm %reconstructed width was 192.5nm

bin_thickness=0.6;

image_width_05=5.01;
image_width_10=4.11;
image_width_30=2.65;
image_width_50=4.3;
image_width_90=6.0;

area_05=image_width_05*bin_thickness;
area_10=image_width_10*bin_thickness;
area_30=image_width_30*bin_thickness;
area_50=image_width_50*bin_thickness;
area_90=image_width_90*bin_thickness;

distance_raw=bin_thickness:bin_thickness:5.4;
distance=zeros(1,length(distance_raw)*2+2);

for jj=1:length(distance)
    if jj==1
        distance(jj)=0;
    elseif jj==length(distance)
        distance(jj)=6.0;
    else
        distance(jj)=distance_raw(floor(jj/2));
    end
end

void_swelling_05=zeros(1,20);
void_swelling_10=zeros(1,20);
void_swelling_30=zeros(1,20);
void_swelling_50=zeros(1,20);
void_swelling_90=zeros(1,20);

for jj=5:7
    name_str=strcat('area_05dpa_bin',num2str(jj));
    void_swelling_05(jj*2)=calculate_swelling(eval(name_str),area_05,thickness_05);
    void_swelling_05(jj*2-1)=calculate_swelling(eval(name_str),area_05,thickness_05);
end
for jj=6:7
    name_str=strcat('area_10dpa_bin',num2str(jj));
    void_swelling_10(jj*2)=calculate_swelling(eval(name_str),area_10,thickness_10);
    void_swelling_10(jj*2-1)=calculate_swelling(eval(name_str),area_10,thickness_10);
end
for jj=1:8
    name_str=strcat('area_30dpa_bin',num2str(jj));
    void_swelling_30(jj*2)=calculate_swelling(eval(name_str),area_30,thickness_30);
    void_swelling_30(jj*2-1)=calculate_swelling(eval(name_str),area_30,thickness_30);
end
for jj=2:8
    name_str=strcat('area_50dpa_bin',num2str(jj));
    void_swelling_50(jj*2)=calculate_swelling(eval(name_str),area_50,thickness_50);
    void_swelling_50(jj*2-1)=calculate_swelling(eval(name_str),area_50,thickness_50);
end
for jj=2:9
    name_str=strcat('area_90dpa_bin',num2str(jj));
    void_swelling_90(jj*2)=calculate_swelling(eval(name_str),area_90,thickness_90);
    void_swelling_90(jj*2-1)=calculate_swelling(eval(name_str),area_90,thickness_90);
end

ave_swelling_05=mean(void_swelling_05(1:16));
ave_swelling_10=mean(void_swelling_10(1:16));
ave_swelling_30=mean(void_swelling_30(1:16));
ave_swelling_50=mean(void_swelling_50(1:16));
ave_swelling_90=mean(void_swelling_90(1:16));

ave_swelling_vec=[ave_swelling_05 ave_swelling_10 ave_swelling_30 ave_swelling_50 ave_swelling_90];

figure()
plot([5 10 30 50 90],ave_swelling_vec*100)

plot_lim=30;

figure('Position',[100 100 500 600])
%make subplots from the bottom up
subplot('Position',[0.12 0.08 0.85 0.18])
plot(distance,void_swelling_90*100,'-','Color',[0 0.5 0],'DisplayName','90 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 plot_lim],'k--')
%legend('90 dpa','Location','NW')
ylim([0 plot_lim])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'ytick',[5 15 25],...
    'LineWidth',1.25)
xlabel({'Distance from surface ($\mu$m)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')
%%%
subplot('Position',[0.12 0.26 0.85 0.18])
plot(distance,void_swelling_50*100,'-','Color',[1 0.5 0],'DisplayName','50 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 plot_lim],'k--')
ylim([0 plot_lim])
%legend('50 dpa','Location','NW')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)

%%%
subplot('Position',[0.12 0.44 0.85 0.18])
plot(distance,void_swelling_30*100,'k-','DisplayName','30 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 plot_lim],'k--')
%legend('30 dpa','Location','NW')
ylim([0 plot_lim])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)
ylabel({'Areal porosity (\%)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')
%%%
subplot('Position',[0.12 0.62 0.85 0.18])
plot(distance,void_swelling_10*100,'b-','DisplayName','10 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 plot_lim],'k--')
%legend('10 dpa','Location','NW')
ylim([0 plot_lim])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)

%%%
subplot('Position',[0.12 0.80 0.85 0.18])
plot(distance,void_swelling_05*100,'r-','DisplayName','5 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 plot_lim],'k--')
%legend('5 dpa','Location','NW')
ylim([0 plot_lim])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)

if do_both

void_swelling_05=void_swelling_05*100;
void_swelling_10=void_swelling_10*100;
void_swelling_30=void_swelling_30*100;
void_swelling_50=void_swelling_50*100;
void_swelling_90=void_swelling_90*100;

%depth data sampled at 0.6um bins

bin_step=0.6;
distance_raw=bin_step:bin_step:5.4;
distance=zeros(1,length(distance_raw)*2+2);

for jj=1:length(distance)
    if jj==1
        distance(jj)=0;
    elseif jj==length(distance)
        distance(jj)=6.0;
    else
        distance(jj)=distance_raw(floor(jj/2));
    end
end

%5 dpa

area_fraction_5=[0 0 0 0 0 0 0 0 1.176 1.176 3.336 3.336 1.796 1.796 0 0 0 0 0 0];

slices_w_voids_5=[1.176 3.336 1.796];

image_width_5=5.01;

total_area_5=image_width_5*4.8;
void_area_5=sum((slices_w_voids_5/100)*bin_step*image_width_5);
void_fraction_5=void_area_5/total_area_5;

%10 dpa

area_fraction_10=[0 0 0 0 0 0 0 0 0 0 4.769 4.769 10.931 10.931 0 0 0 0 0 0];

slices_w_voids_10=[4.769 10.931];

image_width_10=4.11;

total_area_10=image_width_10*4.8;
void_area_10=sum((slices_w_voids_10/100)*bin_step*image_width_10);
void_fraction_10=void_area_10/total_area_10;

%30 dpa

area_fraction_30=[2.061 2.061 3.196 3.196 1.734 1.734 1.808 1.808 2.390 2.390 5.253 5.253 6.301 6.301 5.627 5.627 0 0 0 0];

slices_w_voids_30=[2.061 3.196 1.734 1.808 2.390 5.253 6.301 5.627];

image_width_30=2.65;

total_area_30=image_width_30*4.8;
void_area_30=sum((slices_w_voids_30/100)*bin_step*image_width_30);
void_fraction_30=void_area_30/total_area_30;

%h50 dpa

area_fraction_50=[0 0 0.865 0.865 13.518 13.518 8.882 8.882 7.670 7.670 14.380 14.380 16.965 16.965 6.466 6.466 0 0 0 0];

slices_w_voids_50=[0.865 13.518 8.882 7.670 14.380 16.965 6.466];

image_width_50=4.3;

total_area_50=image_width_50*4.8;
void_area_50=sum((slices_w_voids_50/100)*bin_step*image_width_50);
void_fraction_50=void_area_50/total_area_50;

%90 dpa

area_fraction_90=[0 0 16.594 16.594 10.153 10.153 11.069 11.069 12.933 12.933 13.471 13.471 15.749 15.749 22.480 22.480 8.962 8.962 0 0];

slices_w_voids_90=[16.594 10.153 11.069 12.933 13.471 15.749 22.480];

image_width_90=6.0;

total_area_90=image_width_90*4.8;
void_area_90=sum((slices_w_voids_90/100)*bin_step*image_width_90);
void_fraction_90=void_area_90/total_area_90;

void_fraction_vec=[void_fraction_5 void_fraction_10 void_fraction_30 void_fraction_50 void_fraction_90];

figure()
plot([5 10 30 50 90],void_fraction_vec)


figure('Position',[100 100 500 600])
%make subplots from the bottom up
subplot('Position',[0.12 0.08 0.85 0.18])
plot(distance,area_fraction_90,'r-','LineWidth',1.5)
hold on
plot(distance,void_swelling_90,'b--','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k:','LineWidth',1.5)
%legend('90 dpa','Location','NW')
ylim([0 30])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'ytick',[5 15 25],...
    'LineWidth',1.25)
xlabel({'Distance from surface ($\mu$m)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')
%%%
subplot('Position',[0.12 0.26 0.85 0.18])
plot(distance,area_fraction_50,'r-','LineWidth',1.5)
hold on
plot(distance,void_swelling_50,'b--','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k:','LineWidth',1.5)
ylim([0 30])
%legend('50 dpa','Location','NW')
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)

%%%
subplot('Position',[0.12 0.44 0.85 0.18])
plot(distance,area_fraction_30,'r-','LineWidth',1.5)
hold on
plot(distance,void_swelling_30,'b--','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k:','LineWidth',1.5)
%legend('30 dpa','Location','NW')
ylim([0 30])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)
ylabel({'Swelling (\%)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')
%%%
subplot('Position',[0.12 0.62 0.85 0.18])
plot(distance,area_fraction_10,'r-','LineWidth',1.5)
hold on
plot(distance,void_swelling_10,'b--','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k:','LineWidth',1.5)
%legend('10 dpa','Location','NW')
ylim([0 30])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)

%%%
subplot('Position',[0.12 0.80 0.85 0.18])
pa=plot(distance,area_fraction_5,'r-','LineWidth',1.5,'DisplayName','Areal Porosity');
hold on
pv=plot(distance,void_swelling_05,'b--','LineWidth',1.5,'DisplayName','Volumetric Swelling');
hold on
plot([4.8 4.8],[0 30],'k:','LineWidth',1.5)
legend([pa pv],'Location','NW')
ylim([0 30])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)
end