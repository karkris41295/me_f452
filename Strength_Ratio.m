function [sr_top, sr_bot] = Strength_Ratio(E1, E2, v12, G12, a, t, thetas, lt, deltaT, Nx, Ny, Mx, My, F1, F2, F11, F22, F66)
%Function which calculates TSAI WU values for Project_Main
%NOTE: IN THIS FUNCTION ALL SR STANDS FOR TSAI WU

N = length(thetas);
H = lt * N/2;
sr_top = zeros(N,1);
sr_bot = zeros(N,1);

[e0, kappa] = epsilon_kappa(E1, E2, v12, G12, a, t, thetas, lt, deltaT, Nx, Ny, Mx, My);

for k = 1:N
    
    thetak = thetas(k);
    [Qk, ~, alphak] = QSalpha(E1, E2, v12, G12, thetak, a, t);
    
    %Getting transfrom ready
    [T1, ~] = transforms(thetak);
    %In an earlier iteration of the code, I had mixed up zk and zk1 for top and bottom
    %check at top
    zk1 =  -H + (k-1)*lt; %z_k-1
    eps_top = e0 + zk1*kappa;%top layer epsilon
    top_stress_g = Qk * (eps_top-alphak*deltaT); %top stress global
    top_stress_l = T1 * top_stress_g; %top stress local

    % check at bottom
    zk = -H + (k)*lt; %z_k
    eps_bot = e0 + zk*kappa;
    bottom_stress_g = Qk * (eps_bot-alphak*deltaT);% Sigma x, y, z
    bottom_stress_l = T1 * bottom_stress_g; %Sigma 1, 2, 3
       
    s1t = top_stress_l(1);
    s2t = top_stress_l(2);
    s12t = top_stress_l(3);
    s1b = bottom_stress_l(1);
    s2b = bottom_stress_l(2);
    s12b = bottom_stress_l(3);
    
    Eqn1 = F1*s1t + F2*s2t + F11*(s1t)^2 + F22*(s2t)^2 + F66*(s12t)^2;
    Eqn2 = F1*s1b + F2*s2b + F11*(s1b)^2 + F22*(s2b)^2 + F66*(s12b)^2;
    
    sr_top(k) = Eqn1; %THESE ARE TSAI WU VALUES BEING STACKED IN ARRAYS
    sr_bot(k) = Eqn2;

end