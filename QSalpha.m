function [Qb, Sb, alphak] = QSalpha(E1, E2, v12, G12, theta, aCTE, tCTE)

%% MAIN

%Creating a Reduced Compliance Matrix

S(1,1) = 1/E1;
S(2,2) = 1/E2;
S(3,3) = 1/G12;
S(1,2) = -v12/E1;
S(2,1) = -v12/E1;

%Creating the Reduced Stiffness Matrix

Q = S^-1;

m = cosd(theta);
n = sind(theta);

%Transformation Matrices 

T1 = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, (m^2-n^2)];
T2 = [m^2, n^2, m*n; n^2, m^2, -m*n; -2*m*n, 2*m*n, (m^2-n^2)];

%Creating Global Q and S Matrices

Qb = T1\Q*T2;
Sb = Qb^-1;

%Creating the Alpha Vector in xy system

alphak = T2\[aCTE; tCTE; 0];
end