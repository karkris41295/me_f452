clc; clear;

% Given Input Load
N = [5;5;0] * 1e3; %N/m
M = [-3000;0;0];  %N

Q = [100, 25, 0; 25, 100, 0; 0, 0, 37.5]*1e9; %Pa
F1 = 0; F2 = 0; F11 = 4e-18; F22 = F11; F66 = 4e-16; %1/Pa^2

t = 2e-3; %Initial Thickness Guess in metres

for i = 1:15
    [T1, T2] = transforms(45); % using transforms function for transform
    Qb = T1\Q*T2;
    
    %ABD Calculations
    A = Qb*t;
    D = Qb*(t^3/6);
    eps0 = A\N; %A inverse times N
    kappa = D\M; %D inverse times M
    
    z_bot = t/2; %checking at bottom
    eps = eps0 + z_bot * kappa;
    
    % Stress in global coordinates
    sigma_x = Qb*eps;
    sigma_1 = T1*sigma_x;
    s1 = sigma_1(1);
    s2 = sigma_1(2);
    s12 = sigma_1(3);
    
    TsaiWu = F11*s1^2 + F22*s2^2 + F66*s12^2 ;
    SR = TsaiWu;
    
    if SR>1.01
        fprintf('Iteration %d Thickness %.4f SR %.3f\n',i, t*1000, SR)
        t= t+t/5; %Increse thickness if TsaiWu is large
        
    elseif SR < .99
        fprintf('Iteration %d Thickness %.4f SR %.3f\n',i, t*1000, SR)
        t= t-t/5; %Reduce thickness if TsaiWu is small
    else
        fprintf('Iteration %d Thickness %.4f SR %.3f\n',i, t*1000, SR)
        break;
    end
end

if (.99<SR && SR<1.01)
    fprintf('\nOptimized Thickness is %.3f mm \n',t*1000)
else
    disp('more iterations needed')    
end    