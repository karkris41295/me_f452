clear all; clc;
%% Constants
E1 = 38.6 * 10^9; %Pa (Axial Modulus)
E2 = 8.27 * 10^9 ; %Pa (Transverse Modulus)
v12 = .26; %Dimensionless (Poisson's Ratio)
G12 = 4.14 * 10^9; %Pa (Shear Modulus)
a = 8.6 * 10^-6; %m/C (Axial CTE)
t = 22.1 * 10^-6; %m/C (Transverse CTE)

thetas = [45, -45, 90, 0, 0, 90, -45, 45] %Deg (Angles of stacking)
N = length(thetas); %no. of layers
lt = 1e-3; %m (Layer Thickness) INITIAL GUESS

H = lt * N/2;
deltaT = -200;
Nx = 5000; Ny = 5000; %Nm^-1
Mx = -3000; %N
My = 0;

xt = 1062e6;
xc = -610e6;
yt = 31e6;
yc = -118e6;
s = 72e6;

F1 = (1/xt + 1/xc);
F2 = (1/yt + 1/yc);
F11 = -1/(xt*xc);
F22 = -1/(yt*yc);
F66 = 1/s^2;

for i = 1:100
    [tsai_top, tsai_bot] = Strength_Ratio(E1, E2, v12, G12, a, t, thetas, lt, deltaT, Nx, Ny, Mx, My, F1, F2, F11, F22, F66);
    
    tsai = max(max([tsai_top, tsai_bot]));
    
    if tsai > 1.01
        fprintf('Iteration %d, Thickness %.4f mm, Max Tsai Wu %.3f\n',i, lt*1000, tsai)
        lt= lt+lt/10; %Increse thickness if TsaiWu is large
    elseif tsai < .99
        fprintf('Iteration %d, Thickness %.4f mm, Max Tsai Wu %.3f\n',i, lt*1000, tsai)
        lt= lt-lt/10; %Reduce thickness if TsaiWu is small
    else
        fprintf('Iteration %d, Thickness %.4f mm, Max Tsai Wu %.3f\n',i, lt*1000, tsai)
        break;
    end
end

if (.99<tsai && tsai<1.01)
    fprintf('\nOptimized Thickness is %.3f mm \n',lt*1000)
else
    disp('more iterations needed')    
end    
