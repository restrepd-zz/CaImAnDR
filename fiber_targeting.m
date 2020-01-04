close all
clear all
delta_angle=1;
d_from_mirror=[0 10 20 30 40]; %cm
delta_beam=(d_from_mirror*10)*sind(delta_angle);
distance_from_second_to_coupler=4.7; %cm
distance_from_first_to_second=21; %cm
d_mirrors=[distance_from_second_to_coupler distance_from_first_to_second+distance_from_second_to_coupler];
delta_beam_mirrors=(d_mirrors*10)*sind(delta_angle);
plot(d_from_mirror,delta_beam,'-b')
hold on
plot(d_mirrors,delta_beam_mirrors,'or')
title('Movement of the laser beam across the coupling lens for 1 degree angle')
xlabel('Distance from mirror (cm)')
ylabel('Change across coupling lens (mm)')

%For three clockwise turns the mirror moves to the left 33 mm at 64 cm from the post
degrees_per_turn=asind(3.3/64)/3