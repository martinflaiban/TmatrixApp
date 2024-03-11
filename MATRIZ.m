function [H,f,ftg,tg,hlp,fmotor,fbox,fc,htf] = MATRIZ(BL,Re,Le,Mm,Cms,Rms,Ds,Lport,Dport,Rap,Vbox,type)
%% Función de cálculo %%
%% Agradecimientos al Ing. Facundo Ramón %%
%% Autor: Flaiban Martín %%
% esta función permite calcular la transferencia completa de un
% altoparlante en gabinete ventilado a partir de sus parámetros T-S
% BL Motor magnético [Tesla x Metro]
% Re Resistencia eléctrica de la bobina [ohms]
% Le Inductancia eléctrica de la bobina [Hy]
% Mm Masa mecánica del diafragma [kg]
% Cms Compliancia mecánica de la suspensión [N/metro]
% Rms resistencia mecánica del altavoz [Kg/s]
% Ds superficie efectiva del diafragma [m^2]
% Lport Largo del port [m]
% Dport Diametro del port [m]
% Rap = Resistencia acústica del port [Ns/m^5]
% Vbox = Volumen del gabinete [m^3]

%% Constantes acusticas
AC_c = 343; % Velocidad del sonido
AC_rho0 = 1.18; % Densidad del aire CNPT
AC_p0 = 101.325e3; % Presión estática 
AC_pref = 20e-6;
AC_Zs = AC_rho0*AC_c;
AC_gamma = 1.4; % Constante de los gases adiabáticos

a = Ds/2;      %radio m
motor_Sd = pi*(Ds/2)^2; %m^2
fmotor = 1/2/pi/sqrt(Cms*Mm) ; %resonancia del parlante al aire libre

%% Vectores auxiliares
k = @(w) w/AC_c;    %numero de onda como funcion de velocidad angular
ka = @(w) k(w)*a;   %numero de onda relacionado con radio de piston

%% Construccion de matrices
T_Re = [1 Re; 0 1];
T_Le = @(w) [1 (1j*w*Le); 0 1];
G_Bl = [0 BL; 1/BL 0];
T_Mm = @(w) [1 (1j*w*Mm); 0 1];
T_Cm = @(w) [1 (1/(1j*w*Cms)); 0 1];
T_Rm = [1 Rms; 0 1];
T_motor = @(w) T_Re*T_Le(w)*G_Bl*T_Mm(w)*T_Cm(w)*T_Rm; % Motor electromagnético 
T_Sd = [motor_Sd 0; 0 1/motor_Sd]; %transformador con factor Sd

%% Impedancia de radiación
% Impedancia de radiacion desde el punto de vista mecanico para un piston 
% plano en pantalla infinita.
% Z_rad: Valor cuantitativo de como el medio (aire) reacciona contra el
% movimiento de una superficie vibrante
% Se modela en el dominio MECANICO como una impedancia en SERIE.
% Ortega, Romero - Electroacustica, Altavoces y Microfonos. Eq. (3.32)

ZM_rad = @(w) pi*a^2*AC_Zs*...
    (1-(besselj(1,2*ka(w))/(ka(w))))+...
    1j*pi*a^2*AC_Zs*(struve2(2*ka(w))/((ka(w))));
T_Zm = @(w) [1 ZM_rad(w);0 1]; %%% Esta matriz es impedancia mecánica 
                               %%%    de radiación del aire
%% Impedancia Acustica Especifica para una onda esferica.
%  Relacion entre presion y velocidad volumetrica para una onda esferica 
%propagandose por el aire una distancia D. (Lampton 1978 Eq. (6.1))

d = 1;    %distancia de medicion en metros
Q = 2;    %factor de directividad (Q=2 -> semi-esfera)

ZA_rad = @(w) 1j*w*AC_rho0*Q/(4*pi*d)*exp(-1j*k(w)*d);

T_Za = @(w) [0 -ZA_rad(w); 1/ZA_rad(w) 0];

% Giro de fase debido a la propagacion en el aire.
Z_delay = @(w) 1j*exp(-1j*k(w)*d);

%% Tipo de gabinete
%% type ==2 vented box , type ==1 closed box , type ==0 infinite baffle %%
if type == 2
    S = pi*(Dport/2)^2 ; % Superficie del port en metros cuadrados
    ca = Vbox/AC_gamma/AC_p0; %compliancia acústica del gabinete
    ma = AC_rho0*Lport/S; %Masa acústica del port
    rap = Rap; %resistencia en serie
    fbox = 1/2/pi/sqrt(ca*ma) ; %resonancia de la caja
    Ca = @(w) 1/(1j*w*ca); % Reactancia del gabinete
    Ma = @(w) 1j*w*ma+rap; %Reactancia del port
    T_Ca = @(w) [1 0;1/Ca(w) 1]; %admitancia del gabinete
    T_Ma = @(w) [1 Ma(w); 0 1]; %impedancia del port
    T_minus = [-1 0;0 -1];
    T_box = @(w) T_Sd*T_minus*T_Ca(w)*T_Ma(w); %vented box matrix
    fc=0;
elseif type == 1
    ca = Vbox/AC_gamma/AC_p0; %compliancia
    Cas = Cms*a^2;
    Ca = @(w) a^2/(1j*w*ca); % Reactancia del gabinete
    alpha= Cas/ca ; 
    fc = fmotor*sqrt(1+alpha);
    T_Ca = @(w) [1 Ca(w);0 1]; %admitancia
    P = @(w) T_motor(w)*T_Ca(w)*T_Sd; %Closed box matrix
    fbox=0;
else
    P = @(w) T_motor(w)*T_Zm(w)*T_Sd; % Infinte bafle matrix
    fbox=0;
    fc=0;
end

    

%% Calculo de transferencia
fs= 16e3 ; % Frecuencia de muestreo
f_1 = 0.1; % Frecuencia inicial
n = 2^16 ; % Número de puntos para el cálculo, Resolución = fs/N
f= f_1:fs/n:fs/2 ; %% definir  la frec para fase lineal
W = 2*pi*f; % Velocidad angular (rad/s)
ZL = zeros(n,1); %Z Load
H = zeros(n,1); %Funcion de transferencia entre voltaje y presion
delay = zeros(n,1); %Delay a substraer
Vin = 2.83; %Tensión de entrada para cálculo de sensibilidad a 8ohms de impedancia nominal

for i = 1:length(f)
    if type ~= 2
        ZL(i) = ZA_rad(W(i));
        T = P(W(i));
        A = T(1,1);
        B = T(1,2);
        H(i) = ZL(i)/(A*ZL(i)+B);
    else
        t = T_Zm(W(i))*T_Sd; % Matriz del motor
        A = t(1,1);
        B = t(1,2);
        C = t(2,1);
        D = t(2,2);
        t_ = T_box(W(i)); % Matriz del gabinete
        A_ = t_(1,1);
        B_ = t_(1,2);
        C_ = t_(2,1);
        D_ = t_(2,2);
        T_sipo = [A+A_+(B*C_+B_*C-B*C-B_*C_)/(D+D_) (B_*D+B*D_)/(D+D_);
            (C*D_+C_*D)/(D+D_) D*D_/(D+D_)]; % Transformación serie-paralelo

        t = T_motor(W(i))*T_sipo*T_Za(W(i));
        A = t(1,1);
        H(i) = 1/A; % Transferencia del sistema en el dominio de la frecuencia con ZL = inf
    end 
    delay(i) = Z_delay(W(i));

end
H(1)=0+0i;
fir = Vin*(H)/AC_pref; % FIR , conjugar y hacer ifft(fir,'simmetry')
fir(1)=0+0i; % 
fir2=flip(conj(fir));
fir=[fir(1:length(f));fir2(length(f)+1:end)];
htf=ifft(fir,'symmetric');
E=sqrt(sum(htf.^2)); 
htf=htf./E;
[M,idx]=max(htf); idx= idx-10 ;
firn = 2047;
fade=round(firn/8);
win=hann(2*fade);
htf= htf(idx:firn+idx);
htf(end-fade:end)=htf(end-fade:end).*win(fade:end);
hlp = firtest(htf,2^11);
hlp(1:fade)=hlp(1:fade).*win(1:fade);
hlp(end-fade:end)=hlp(end-fade:end).*win(fade:end);

%%%%%%%
H = Vin*(H./delay)/AC_pref; %Calculo de SPL para visualización
H2=flip(conj(H));
H=[H(1:length(f));H2(length(f)+1:end)];
%%%%%%%
ftg= 0:fs/n:fs/2 ; %% definir  la frec
phase=rad2deg(angle(H(1:length(ftg))));
tg = GroupDelay(ftg,phase)*1000 ;
end
    