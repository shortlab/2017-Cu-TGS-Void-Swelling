function [swelling]=calculate_swelling(void_list,area,thickness)
%all calculations done in units of microns

diameter=2*sqrt(void_list/pi)+.007;

d_cubed=diameter.^3;

swole=(pi/6)*sum(d_cubed);

swelling=swole/(area*thickness-swole);

end