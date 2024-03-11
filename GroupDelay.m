function t = GroupDelay(f,phase) 
% Cálculo de Retardo de grupo 
%Hecho por Christopher J. Struck - 3 August 2000 
k = length(phase); 
% Unwrap phase data 
for n = 2:k 
while phase(n) - phase(n - 1) <= -180 
 phase(n) = phase(n) + 360; 
 end 
 while phase(n) - phase(n - 1) >= 180 
 phase(n)= phase(n) - 360; 
 end 
end 
% Calculate Group Delay 
for n = 2:k-1 
 t(n) = (-1/720) * (((phase(n) - phase(n - 1)) / (f(n) - f(n - 1)))... 
 + ((phase(n + 1) - phase(n)) / (f(n + 1) - f(n)))); 
end 
t(1) = (-1/360) * (((phase(2) - phase(1))/(f(2) - f(1)))); 
t(k) = (-1/360) * (((phase(k) - phase(k - 1))/(f(k) - f(k - 1)))); 
end