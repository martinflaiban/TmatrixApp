%% Funcion de Struve2
%Approximation of the Struve function H1 occurring in impedance calculation
%R.M AARTS & AUGUSTUS JANSSEN

function H = struve2(z)
H = (2/pi)-besselj(0,z)+((16/pi) -5)*(sin(z)./z)+(12-(36/pi))*(1-cos(z))./(z.^2);
end