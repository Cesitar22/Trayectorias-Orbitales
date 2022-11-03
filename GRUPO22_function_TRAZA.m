function [lat_N, long_N, t] = GRUPO22_function_TRAZA(i1,a,long_0, lat_0, N, points_traza)

%==========================================================================
%
% Función desarrollada por el grupo 22, 4º curso GIA VA:
%
%   José Luis Flores Alameda      - joseluis.flores.alameda@alumnos.upm.es
%   Larissa Ruiz Sánchez          - larissa.ruiz.sanchez@alumnos.upm.es
%   César Eduardo Jiménez Gámez   - cesar.jgamez@alumnos.upm.es
%   Pedro Mederos Pulido          - p.mederos@alumnos.upm.es
%   Cristina Bernal Retamar       - cristina.bernalr@alumnos.upm.es
%
% ------
% INPUT:
% ------
%   i1           - (1x1) inclinación de la órbita [deg]
%   a            - (1x1) distancia del satélite al centro de la Tierra [m]
%   long_0       - (1x1) longitud del nodo ascendente de la primera traza en el sistema de referencia ECEF [deg]
%   lat_0        - (1x1) latitud inicial de la primera traza en el sistema de referencia ECEF [deg]
%   N            - (1x1) número de trazas a representar
%   points_traza - (1x1) número de divisiones utilizadas para definir el vector tiempo de cada traza
%
%==========================================================================

%% Datos

%Datos de la Tierra
w_T = 7.27*10^-5;               %velocidad angular de rotación de la Tierra [rad/s]
r_T = 6371e3;                   %radio Tierra [m]
mu = 398601e9;                  %constante gravitacional de la Tierra[m3/s2]

%Datos de la órbita
h = a-r_T;                      %semieje mayor de la órbita [m]
T = 2*pi*sqrt(a^3/mu);          %periodo de la órbita [s]
lambda = w_T*2*pi*sqrt(a^3/mu); %longitud entre nodos ascendentes [rad]
long_0 = deg2rad(long_0);       %longitud inicial [rad]
lat_0 = deg2rad(lat_0);         %latitud inicial [rad]
i1 = deg2rad(i1);               %inclinación de la órbita [rad]

%Datos del satélite
eta = deg2rad(40);            %semiángulo del cono [rad]
r_b = h*tan(eta);               %radio de visión barrido por el satélite [m]
R = h*tan(eta);                 %radio de visibilidad del satélite [m]

%Datos de trigonometría esférica
rho = asin(r_T/a);              %radio angular terrestre [rad]
eps = acos(sin(eta)/sin(rho));  %elevación [rad]
lambda_ = pi/2 - eta - eps;     %radio angular [rad]


%% Calculo de la órbita

t = linspace(0,T,points_traza);                                          %vector de tiempo [s]

%Definición de la longitud y la latitud de los puntos de la traza
long_N = 0;
lat_N = 0;
for j = 1:N
    for i=1:length(t)
        u(i) = sqrt(mu/(a^3))*t(i);
        lat(i) = asin(sin(i1)*sin(u(i)));

        if (i>=0*length(t)) && (i<=1/4*length(t))                   %Corrección del 1º cuadrante
        long(i) = long_0 - w_T*t(i) + atan(cos(i1)*tan(u(i)));

        elseif (i>1/4*length(t)) && (i<=1/2*length(t))              %Corrección del 2º cuadrante
        long(i) = long_0 - w_T*t(i) + atan(cos(i1)*tan(u(i)))+pi;

        elseif (i>1/2*length(t)) && (i<=3/4*length(t))              %Corrección del 3º cuadrante
        long(i) = long_0 - w_T*t(i) + atan(cos(i1)*tan(u(i)))-pi;

        elseif (i>3/4*length(t)) && (i<=1*length(t))                %Corrección del 4º cuadrante
        long(i) = long_0 - w_T*t(i) + atan(cos(i1)*tan(u(i)));
        end

        %Corrección de la traza al intervalo [-180º, 180º]
        if long(i)<-pi
            long(i) = long(i) + 2*pi;
        end

        if long(i)>pi
            long(i) = long(i) - 2*pi;
        end

    end
    
    long_0 = long(end);        %la longitud inicial de la siguiente traza es la longitud final de la anterior
    long_N = [long_N long];    %se almacena en un vector los valores de longitud de todas las trazas
    lat_N = [lat_N lat];       %se almacena en un vector los valores de latitud de todas las trazas
    
end

    long_N = long_N(2:end);
    lat_N = lat_N(2:end);

    for pos_long_0 = 1:length(lat_N)
        if lat_0<lat_N(pos_long_0)
            break
        end
            pos_long_0;
    end

    lat_N_auxiliar = lat_N(1:pos_long_0);
    long_N_auxiliar = long_N(1:pos_long_0);
    
    lat_N = lat_N(pos_long_0:end);
    long_N = long_N(pos_long_0:end);
    long_N = long_N - long_N(1) + deg2rad(long_0);
    
    lat_N = [lat_N lat_N_auxiliar];
    long_N = [long_N long_N_auxiliar+long_N(end)];
    
    %Corrección de la traza al intervalo [-180º, 180º]
    for i = 1:length(long_N)
        if long_N(i)<-pi
            long_N(i) = long_N(i) + 2*pi;
        end

        if long_N(i)>pi
            long_N(i) = long_N(i) - 2*pi;
        end
    end

%% Representación de la órbita
figure
plot_ground_track(rad2deg(lat_N),rad2deg(long_N),'Earth')
title ('Traza del punto subsatélite')


function plot_ground_track(lat,lon,opts,planet)
    
    % ------------------------------------
    % Sets (or defaults) plotting options.
    % ------------------------------------
    
    % sets line color (defaults to default MATLAB color)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'color')
        color = [0,0.4470,0.7410];
    else
        color = opts.color;
    end
    
    % sets line style (defaults to solid line)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'line_style')
        line_style = '-';
    else
        line_style = opts.line_style;
    end
    
    % sets line width (defaults to 1.5)
    if (nargin < 3) || isempty(opts) || ~isfield(opts,'line_width')
        line_width = 1.5;
    else
        line_width = opts.line_width;
    end
    
    % sets background (defaults to Earth coastlines)
    if (nargin < 4) || isempty(planet)
        planet = 'Earth Coastlines';
    end
    
    % --------------------------------------
    % Draws background of ground track plot.
    % --------------------------------------
    
    % sets background of ground track plot
    if strcmpi(planet,'Earth Coastlines')
        
        % loads Earth topographic data
        load('topo.mat','topo');

        % rearranges Earth topopgrahic data so it goes from -180 to 180 
        % deg longitude
        topoplot = [topo(:,181:360),topo(:,1:180)];

        % plots Earth map by making a contour plot of topographic data at
        % elevation of 0
        contour(-180:179,-90:89,topoplot,[0,0],'black');
        
    elseif ~strcmpi(planet,'Blank')
        
        % reads in image
        if strcmpi(planet,'Earth Cloudy')
            cdata = imread('images/earth.png')+imread('images/clouds.png');
        elseif strcmpi(planet,'Earth Night')
            cdata = imread('images/earthnight.png');
        elseif strcmpi(planet,'Earth Night Cloudy')
            cdata = imread('images/earthnight.png')+0.1*...
                imread('images/clouds.png');
        else
            cdata = imread(strcat('images/',lower(planet),'.png'));
        end
        
        % plots background
        image('CData',flipud(cdata),'XData',[-180,180],'YData',[-90,90]);
        
        % determines grid style based on background
        if strcmpi(planet,'Sun')
            grid_line_width = 1;
            grid_color = 'k';
        elseif strcmpi(planet,'Moon')
            grid_line_width = 0.75;
            grid_color = 'y';
        elseif strcmpi(planet,'Mercury')
            grid_line_width = 0.75;
            grid_color = 'y';
        elseif strcmpi(planet,'Venus')
            grid_line_width = 0.5;
            grid_color = 'k';
        elseif strcmpi(planet,'Earth')
            grid_line_width = 0.5;
            grid_color = 'y';
        elseif strcmpi(planet,'Earth Cloudy')
            grid_line_width = 1.5;
            grid_color = 'y';
        elseif strcmpi(planet,'Earth Night')
            grid_line_width = 0.5;
            grid_color = [0.5,0.5,0.5];
        elseif strcmpi(planet,'Earth Night Cloudy')
            grid_line_width = 0.5;
            grid_color = [0.5,0.5,0.5];
        elseif strcmpi(planet,'Mars')
            grid_line_width = 0.75;
            grid_color = 'y';
        elseif strcmpi(planet,'Jupiter')
            grid_line_width = 1;
            grid_color = 'k';
        elseif strcmpi(planet,'Saturn')
            grid_line_width = 0.65;
            grid_color = 'k';
        elseif strcmpi(planet,'Uranus')
            grid_line_width = 0.5;
            grid_color = 'k';
        elseif strcmpi(planet,'Neptune')
            grid_line_width = 0.5;
            grid_color = 'w';
        elseif strcmpi(planet,'Pluto')
            grid_line_width = 0.65;
            grid_color = 'g';
        end
        
        % manually adds grid
        hold on;
        for i = 1:11
            plot([-180+i*30,-180+i*30],[-90,90],'color',grid_color,...
                'linewidth',grid_line_width,'linestyle',':');
        end
        for i = 1:5
            plot([-180,180],[-90+i*30,-90+i*30],'color',grid_color,...
                'linewidth',grid_line_width,'linestyle',':');
        end
    
    end
    
    % ----------------------------------------
    % Plotting ground track / axis formatting.
    % ----------------------------------------
    
    % determines indices where ground track crosses figure border (i.e.
    % from 180 to -180 or -180 to 180) to avoid "jumps" in the plot
    j = [];
    for i = 2:length(lon)
        if ((lon(i) > 170) && (lon(i-1) < -170)) || ((lon(i) < -170) &&...
                (lon(i-1) > 170))
            j = [j,i];
        end
    end
    
    % adds last index to "j" in order to plot ground track between last
    % figure border crossing and the last input longitude
    j = [j,length(lon)];
    
    % plots groundtrack (starts new plot every time ground track crosses
    % left border or right border)
    hold on;
    plot(lon(1:(j(1)-1)),lat(1:(j(1)-1)),'color',color,'linestyle',...
        line_style,'linewidth',line_width);
    for i = 1:(length(j)-1)
        plot(lon(j(i):(j(i+1)-1)),lat(j(i):(j(i+1)-1)),'color',color,...
           'linestyle',line_style,'linewidth',line_width,...
           'handlevisibility','off');
    end

    % axis formatting
    axis equal
    grid on
    if strcmpi(planet,'Earth Coastlines') || strcmpi(planet,...
            'Blank')
        ax = gca;
        ax.GridColor = [0.35,0.35,0.35];
        ax.GridAlpha = 1;
        ax.GridLineStyle = ':';
    end
    xlim([-180,180]);
    xticks(-180:30:180);
    ylim([-90,90]);
    yticks(-90:30:90);
    xlabel('Longitude, $\lambda\;[^{\circ}]$','interpreter','latex',...
        'fontsize',18);
    ylabel('Latitude, $\phi\;[^{\circ}]$','interpreter','latex',...
        'fontsize',18);
    
end

end