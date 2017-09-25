%histogram of all swelling data

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

figure()
plot(distance,area_fraction_5,'DisplayName','5 dpa')
hold on
plot(distance,area_fraction_10,'DisplayName','10 dpa')
hold on
plot(distance,area_fraction_30,'DisplayName','30 dpa')
hold on
plot(distance,area_fraction_50,'DisplayName','50 dpa')
hold on
plot(distance,area_fraction_90,'DisplayName','90 dpa')
legend('show')

void_fraction_vec=[void_fraction_5 void_fraction_10 void_fraction_30 void_fraction_50 void_fraction_90];

figure()
plot([5 10 30 50 90],void_fraction_vec)


figure('Position',[100 100 500 600])
%make subplots from the bottom up
subplot('Position',[0.12 0.08 0.85 0.18])
plot(distance,area_fraction_90,'-','Color',[0 0.5 0],'DisplayName','90 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k--')
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
plot(distance,area_fraction_50,'-','Color',[1 0.5 0],'DisplayName','50 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k--')
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
plot(distance,area_fraction_30,'k-','DisplayName','30 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k--')
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
ylabel({'Areal porosity (\%)'},...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',20,...
    'FontName','Times')
%%%
subplot('Position',[0.12 0.62 0.85 0.18])
plot(distance,area_fraction_10,'b-','DisplayName','10 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k--')
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
plot(distance,area_fraction_5,'r-','DisplayName','5 dpa','LineWidth',1.5)
hold on
plot([4.8 4.8],[0 30],'k--')
%legend('5 dpa','Location','NW')
ylim([0 30])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Times',...
    'Xticklabel',{},...
    'ytick',[5 15 25],...
    'LineWidth',1.25)