% i1=40;
% a = 10000e3;
% lambda_NA_0=0;
% lat_0=33;
% N=2;
% points_traza=200;
% 
% [lat_N, long_N, t] = GRUPO22_function_TRAZA(i1,a,lambda_NA_0,lat_0,N,points_traza);
% title('Traza');

longitud_ps = 175; %grados
latitud_ps = -60; %grados
h = 880.7e3;  %metros
eta = 50;   %grados
representar = 1;

r_T = 6371e3;
a = h + r_T;
rho = asin(r_T/a);


GRUPO22_function_SWATH(longitud_ps, latitud_ps, h, eta, representar)

eta = deg2rad(eta);
if eta<=rho
    eps = acos(sin(eta)/sin(rho));
    lambda1_ = pi/2 - eta - eps;     %radio angular [rad]
    lambda1 = rad2deg(lambda1_);
    display('Lambda es menor que lambda_0:')
    lambda1
elseif eta>rho
    lambda_0 = pi/2 - rho;
    lambda2_ = lambda_0;
    lambda2 = rad2deg(lambda2_);
    display('Lambda es mayor que lambda_0:')
    lambda2
end

