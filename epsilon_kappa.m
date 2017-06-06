function [e0, kappa] = epsilon_kappa(E1, E2, v12, G12, a, t, thetas, lt, deltaT, Nx, Ny, Mx, My)
%% MAIN

A = zeros(3);
B = zeros(3);
D = zeros(3);
Nt = 0; Mt = 0;
N = length(thetas);
H = lt * N/2;

for k = 1:N
    thetak = thetas(k);
    [Qk, Sb, alphak] = QSalpha(E1, E2, v12, G12, thetak, a, t);
    zk = -H + (k)*lt; %z_k
    zk1 = -H + (k-1)*lt; %z_k-1
    
    Ak = Qk*(zk-zk1);
    A = A + Ak;
    
    Ntk = deltaT*(Qk * alphak)*lt; %Calculate Nthermal and Mthermal
    Mtk = .5*deltaT*((Qk*alphak)*(zk^2-zk1^2));
    
    Nt = Nt + Ntk;
    Mt = Mt + Mtk;
    
    %This matrix does not match up with the Test II solution. Everything else works
    Bk = Qk*((zk^2)-(zk1^2))/2;
    B = B + Bk;
    
    Dk = Qk*((zk^3)-(zk1^3))/3;
    D = D + Dk;
end
  
ABD = [A, B; B, D];

final = [Nt(1) + Nx; Nt(2) + Ny; Nt(3); Mt(1) + Mx; Mt(2)+ My; Mt(3)];
ek = ABD\final;
e0 = ek(1:3);
kappa = ek(4:6);




