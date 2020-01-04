%Generates an Airy pattern offset by a couple of mm
close all
clear all
x = -55:0.2:55; %Size of the image and step-size, arbitrary units
I0 = 10;  %Peak amplitude of the disk
%Obtain intensities using Bessel function of the first kind
I = I0*(2*besselj(1,x)./(x)).^2; 
mult_factor=12;  %Makes the pattern more tight
mm_offset=2;
distance_mm=(x/mult_factor)-mm_offset;
figure(1)
plot(distance_mm,I)
hold on
%50 nW line
plot([distance_mm(1) distance_mm(end)],[0.05 0.05],'-r')
ylabel('mW')
xlabel('mm')

