function [T1, T2] = transforms(theta)

%function which gives a T1 and T2 given a theta

m = cosd(theta);
n = sind(theta);

%Transformation Matrices 

T1 = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, (m^2-n^2)];
T2 = [m^2, n^2, m*n; n^2, m^2, -m*n; -2*m*n, 2*m*n, (m^2-n^2)];