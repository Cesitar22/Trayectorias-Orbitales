%==========================================================================
%
% Programa desarrollado por el grupo 22, 4º curso GIA VA:
%
%   José Luis Flores Alameda      - joseluis.flores.alameda@alumnos.upm.es
%   Larissa Ruiz Sánchez          - larissa.ruiz.sanchez@alumnos.upm.es
%   César Eduardo Jiménez Gámez   - cesar.jgamez@alumnos.upm.es
%   Pedro Mederos Pulido          - p.mederos@alumnos.upm.es
%   Cristina Bernal Retamar       - cristina.bernalr@alumnos.upm.es
%
%==========================================================================

clc
clear all
close all

%% Datos del problema

%Región de interés que debe cubrir la misión
lon_min_deg = -9;                               %longitud mínima [deg]
lon_max_deg = 3;                                %longitud máxima [deg]
lat_min_deg = 36;                               %latitud mínima [deg]
lat_max_deg = 44;                               %latitud máxima [deg]

%Definición de las inclinaciones del problema
INC_a = lat_max_deg - 0.8*lon_max_deg;          %inclinación del apartado a [deg]
INC_b = 85;                                     %inclinación del apartado b [deg]
INC_varias = [INC_a, INC_b];                    %vector de inclinaciones del ejercicio [deg]


%==========================================================================
%=============================     INPUT       ============================
%==========================================================================
i1 = INC_b;                                     %selección del ángulo de inclinación ([deg] en el intervalo [0, 180]) que se desea utilizar para pintar la traza; debe introducirse INC_a o INC_b
alt_min = 100e3;                                %mínima altitud en [m] considerada como válida para las órbitas RGT que cumplen el criterio de revisita cada 5 días o menos

%Selección del punto de lanzamiento para representación de la traza
lambda_NA_0 = 0;                                %valor en [deg], en el intervalo [-180, 180]; longitud del nodo ascendente de la primera traza en el sistema de referencia ECEF
lat_0 = 0;                                      %valor en [deg], en el intervalo [-90, 90]; latitud inicial de la primera traza en el sistema de referencia ECEF
%==========================================================================
%==========================================================================


%Datos de la Tierra
T_T = 23*60*60+56*60+4;                         %periodo orbital siderio de la tierra [s]
r_T = 6378e3;                                   %radio de la Tierra en [m]
w_T = 7.27*10^-5;                               %velocidad angular de rotación de la Tierra [rad/s]
mu = 398601e9;                                  %constante gravitacional de la Tierra[m^3/s^2]
J2 = 1.0827*10^(-3);                            %perturbación debida al achatamiento de la Tierra

%Datos del satélite
eta = 10.5;                                     %semiángulo del cono de visibilidad [deg]

%Parámetros de cálculo y representación
error_representar = 10.5;                       %error utilizado para dibujar la cobertura entorno a la región de interés [deg]
points_traza = 100;                             %número de divisiones utilizadas para definir el vector tiempo de cada traza


%% APARTADO 3, ÓRBITAS RGT KEPLERIANAS

%Inicialización de variables
k_day2r = (1:5);                                %días hasta la repetición (RGT)
k_orb2r = (1:1000);                             %órbitas hasta la repetición (RGT)

semiejes_buenos = 0;                            %inicialización del vector de los semiejes válidos
k_orbPday_buenos = 0;                           %inicialización del vector de los k_orbPday correspondientes a los semiejes válidos
k_day2r_buenos = 0;                             %inicialización del vector de los k_day2r correspondientes a los semiejes válidos

for i = 1:length(k_day2r)
    
    k_day2r_i = k_day2r(i);         
    
    semiejes_ok = 0;
    k_obrPday_ok = 0;
    k_day2r_ok = 0;
    
    for j = 1:length(k_orb2r)
        
        k_orb2r_j = k_orb2r(j);
        k_orbPday_j = k_orb2r_j/k_day2r_i;
        
        a = (mu*(T_T/(2*pi*k_orbPday_j))^2)^(1/3);      %calculo del semieje mayor [m] de las órbitas Keplerianas
        
        if a>r_T+alt_min && a<r_T+1000e3                %coge los valores de 'a' que son validos para las orbitas de alt_min a 1000e3m de altitud (se suma el radio terrestre pues es el semieje)
            semiejes_ok(j) = a;                         %se almacena en un vector los valores de 'a' válidos
            k_obrPday_ok(j) = k_orbPday_j;              %se almacena en un vector los valores de 'k_orbPday_j' corresponientes a los 'a' válidos
            k_day2r_ok(j) = k_day2r_i;                  %se almacena en un vector los valores de 'k_day2r_i' corresponientes a los 'a' válidos
        end
        
    end
    
    semiejes_buenos = [semiejes_buenos semiejes_ok];    %vector que almacena los semiejes válidos (a en [m])  para todos los posibles valores de 'k_day2r' 
    k_orbPday_buenos = [k_orbPday_buenos k_obrPday_ok]; %vector que almacena los k_orbPday correspondientes a los semiejes válidos
    k_day2r_buenos = [k_day2r_buenos k_day2r_ok];       %vector que almacena los k_day2r correspondientes a los semiejes válidos
    
end

matriz_buena = [semiejes_buenos; k_orbPday_buenos; k_day2r_buenos];

k_orbPday = k_orbPday_buenos(find(k_orbPday_buenos ~= 0)); %se quitan los ceros al vector k_orbPday_buenos
a_kepler = (semiejes_buenos(find(semiejes_buenos ~= 0)));  %se quitan los ceros al vector semiejes_buenos, [m]

display('El número de órbitas RGT que cumplen el requisito de revisita son:')
num_orb = length(a_kepler)                                 %número de órbitas válidas

altitudes_buenas = (a_kepler-r_T)*10^-3;                   %altitudes en [km] de las órbitas Keplerianas válidas

%% APARTADO 3, ÓRBITAS RGT KEPLERIANAS CORREGIDAS CON J2

a_kepler_J2 = zeros(length(INC_varias),length(a_kepler));                       %inicialización de la matriz de órbitas RGT keplerianas con la corrección de J2; una fila por cada inclinación

for j = 1:length(INC_varias)
    
    inc_i = deg2rad(INC_varias(j));

    for i = 1:length(semiejes_buenos)                                           %variamos para todas las 'a' de Kepler obtenidas anteriormente, utilizando sus correspondientes k_orbPday

        a_kepler_i = semiejes_buenos(i);
        k_orbPday_i = k_orbPday_buenos(i);

        a_real = a_kepler_i;                                                    %inicialización de 'a', el semieje mayor
        K=1;
        for K=1:1000                                                            %para cada 'a' de Kepler, se calcula la 'a' correspondiente corregida con J2 
            e = 0;                                                              %excentricidad [adim]
            p =  a_real*(1-e^2);                                                %parámetro orbital [m]
            n = sqrt(mu/(a_real^3));                                            %velocidad angular media [rad/s]
            omega_mayus = -3*J2*n*r_T^2*cos(inc_i)/(2*p^2);                     %regresión de los nodos instantánea [rad/s]
            omega_minus = 3*J2*n*r_T^2*(4-5*(sin(inc_i))^2)/(4*p^2);            %avance del perigeo instantáneo [rad/s]
            M1_sec = 3*J2*n*r_T^2*sqrt(1-e^2)*(2-3*(sin(inc_i))^2)/(4*p^2);     %variación de la anomalía media instantánea [rad/s]
            T_omega = (T_T/k_orbPday_i)*(1/(1 - omega_mayus/w_T));
            term_1 = 2*pi/T_omega - M1_sec - omega_minus;                       %término 1 de la ecuación de T_omega
            a_iter = (mu/term_1^2)^(1/3);                                       %se calcula la nueva 'a'
            
            if abs(a_real-a_iter) < 10^-4
              break
            else
              a_real = a_iter;
              K = K+1;
            end
        end
        
        a_kepler_J2(j,i) = a_real;
       
    end
end

a_inc_a = a_kepler_J2(1,:);                 %semiejes [m] corregidos de las órbitas RGT para la inclinación a
a_inc_b = a_kepler_J2(2,:);                 %semiejes [m] corregidos de las órbitas RGT para la inclinación b
a_inc_a_ = (a_inc_a(find(a_inc_a > 0)));    %semiejes [m] corregidos de las órbitas RGT para la inclinación a, eliminando los NaN
a_inc_b_ = (a_inc_b(find(a_inc_b > 0)));    %semiejes [m] corregidos de las órbitas RGT para la inclinación b, eliminando los NaN

altitudes_inc_a = (a_inc_a - r_T)*10^-3;     %altitudes [km] corregidas de las órbitas RGT para la inclinación a
altitudes_inc_b = (a_inc_b - r_T)*10^-3;     %altitudes [km] corregidas de las órbitas RGT para la inclinación b
altitudes_inc_a_ = (a_inc_a_ - r_T)*10^-3;   %altitudes [km] corregidas de las órbitas RGT para la inclinación a, eliminando los NaN
altitudes_inc_b_ = (a_inc_b_ - r_T)*10^-3;   %altitudes [km] corregidas de las órbitas RGT para la inclinación b, eliminando los NaN

display('Las altitudes en [km] de las órbitas RGT corregidas con J2, para la inclinación A, son:')
altitudes_inc_a_
display('Las altitudes en [km] de las órbitas RGT corregidas con J2, para la inclinación B, son:')
altitudes_inc_b_

%% APARTADO 4, MÍNIMA MASA DE PROPULSANTE PARA LAS MANIOBRAS DE MANTENIMIENTO
%Se elige la órbita de mayor altitud (menor resistencia aerodinámica) para masa minima de propulsante de mantenimiento de órbita 
%JUSTIFICACIÓN: a_D = -rho*V^2/(2*B) A mayor altitud --> menor densidad y menor velocidad, por tanto la resistencia aerodinámica es menor cuanto mayor sea la altitud

[orb_min_prop_a, index_a] = max(altitudes_inc_a);        %se elige el valor máximo de altitud [km] y su posición en el vector
display('La altitud en [km] de la órbita RGT para mínima masa de propulsante, para la inclinación A, es:')
orb_min_prop_a
k_orbPday_a = k_orbPday_buenos(index_a);                 %k_orbPday_min_prop correspondiente a dicho máximo de altitud
k_day2r_a = k_day2r_buenos(index_a);                     %días hasta la repetición, k_day2r_min_prop, para la inclinación A
k_orb2r_a = k_orbPday_a*k_day2r_a;                       %órbitas hasta la repetición, k_orb2r_min_prop, para la inclinación A

[orb_min_prop_b, index_b] = max(altitudes_inc_b);        %se elige el valor máximo de altitud [km] y su posición en el vector
display('La altitud en [km] de la órbita RGT para mínima masa de propulsante, para la inclinación B, es:')
orb_min_prop_b
k_orbPday_b = k_orbPday_buenos(index_b);                 %k_orbPday_min_prop correspondiente a dicho máximo de altitud
k_day2r_b = k_day2r_buenos(index_b);                     %días hasta la repetición, k_day2r_min_prop, para la inclinación A
k_orb2r_b = k_orbPday_b*k_day2r_b;                       %órbitas hasta la repetición, k_orb2r_min_prop, para la inclinación A

orb_min_prop_tot = max(orb_min_prop_a,orb_min_prop_b);   %órbita RGT de mayor altitud entre ambas inclinaciones posibles

%% APARTADO 5, COBERTURA
lat_ss_deg = (lat_min_deg + lat_max_deg)/2;    %latitud del punto subsatélite intermedio de la región de interés [deg]
lon_ss_deg = (lon_min_deg + lon_max_deg)/2;    %longitud del punto subsatélite intermedio de la región de interés [deg]
delta_lon = abs(lon_min_deg - lon_max_deg)/2;  %semidistancia entre las longitudes mínima y máxima [deg]
delta_lat = abs(lat_min_deg - lat_max_deg)/2;  %semidistancia entre las longitudes mínima y máxima [deg]


%Representación de la órbita para la inclinación seleccionada
if i1 == INC_a
    h = orb_min_prop_a*10^3;                       %altitud del satélite [m] para el caso a
    N = k_orb2r_a;                                 %número de trazas a representar para el caso a
elseif i1 == INC_b
    h = orb_min_prop_b*10^3;                       %altitud del satélite [m] para el caso b
    N = k_orb2r_b;                                 %número de trazas a representar para el caso b
end

a = h + r_T;                                   %distancia del satélite al centro de la Tierra [m]


% Dibujar la traza
[lat_N, long_N, t] = GRUPO22_function_TRAZA(i1,a,lambda_NA_0,lat_0,N,points_traza);
title('Traza');
longitud = rad2deg(long_N);
latitud = rad2deg(lat_N);

hold on
plot([-delta_lon delta_lon] + lon_ss_deg, (delta_lat+lat_ss_deg)*[1 1], 'k', 'linewidth', 1.5), hold on
plot([-delta_lon delta_lon] + lon_ss_deg, (-delta_lat+lat_ss_deg)*[1 1], 'k', 'linewidth', 1.5), hold on
plot((delta_lon+lon_ss_deg)*[1 1], [-delta_lat delta_lat] + lat_ss_deg, 'k', 'linewidth', 1.5), hold on
plot((-delta_lon+lon_ss_deg)*[1 1], [-delta_lat delta_lat] + lat_ss_deg, 'k', 'linewidth', 1.5), hold on


% Dibujar el swath
figure
for i = 1:length(longitud)
    lat_1 = latitud(i);
    long_1 = longitud(i);
    
    if long_1>=lon_min_deg-error_representar && long_1<=lon_max_deg+error_representar && lat_1>=lat_min_deg-error_representar && lat_1<=lat_max_deg+error_representar;
        GRUPO22_function_SWATH(long_1, lat_1, h, eta, 1); hold on;
        warning('off','all');
    end
    
end
title('Cobertura')

hold on
plot([-delta_lon delta_lon] + lon_ss_deg, (delta_lat+lat_ss_deg)*[1 1], 'k', 'linewidth', 1.5), hold on
plot([-delta_lon delta_lon] + lon_ss_deg, (-delta_lat+lat_ss_deg)*[1 1], 'k', 'linewidth', 1.5), hold on
plot((delta_lon+lon_ss_deg)*[1 1], [-delta_lat delta_lat] + lat_ss_deg, 'k', 'linewidth', 1.5), hold on
plot((-delta_lon+lon_ss_deg)*[1 1], [-delta_lat delta_lat] + lat_ss_deg, 'k', 'linewidth', 1.5), hold on

