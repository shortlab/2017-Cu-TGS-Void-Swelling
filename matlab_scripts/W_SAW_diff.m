%Script to process the fine-meshed calculation of Im[G33] for {001}
%tungsten. The mesh in this case varies from 600x600 to 2000x2000 points and the slowness
%bounds are [0,0.4] us/mm in both dimensions

%Base directory will need to be changed. Calculations on different mesh sizes provided 
%in /analytical_calculations/W_mesh_size/
cd('/Users/Cody/Documents/MIT/Short_Lab/LSAW/A_Every Code/EveryFixed/W_mesh_size');

diff=zeros(1,15);
diff_index=1;

for kk=600:100:2000
    if kk<1000
        num_str=strcat(num2str(0),num2str(kk));
    else
        num_str=num2str(kk);
    end
    
    file_str=strcat('W_quarter_100_',num_str,'.txt');
    
    slow=dlmread(file_str);
    [sx,sy]=size(slow);
    
    slow_filtered=zeros(sx,sy);
    num_pts=0;
    
    for i=1:sy
        for j=1:sx
            if slow(i,j)>60 %Threshold value gets whole SAW line
                slow_filtered(i,j)=100;
            else
                slow_filtered(i,j)=0;
                num_pts=num_pts+1;
            end
        end
    end
    
    slow_mags=zeros(num_pts,3);
    k=1;
    
    for i=1:sy
        for j=1:sx
            if slow_filtered(i,j)==0
                x_units=(j-0.5)*(0.4/(sy));
                y_units=(i-0.5)*(0.4/(sx));
                slow_mags(k,1)=sqrt(x_units^2+y_units^2);
                slow_mags(k,2)=j;
                slow_mags(k,3)=i;
                k=k+1;
            end
        end
    end
    
    max_speed=(1/(min(slow_mags(:,1))));
    min_speed=(1/(max(slow_mags(:,1))));
    
    diff(diff_index)=max_speed/min_speed-1;
    diff_index=diff_index+1;
    
    if kk==2000
        delta_s=0.4/kk;
        min_s=0.5*delta_s;
        max_s=(kk-0.5)*delta_s;
        figure()
        contourf(min_s:delta_s:max_s,min_s:delta_s:max_s,slow_filtered,10,'EdgeColor','none');
        hold on
        colormap(gray)
    end
end

display(diff)