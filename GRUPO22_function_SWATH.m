function [xx_, y] = GRUPO22_function_SWATH(longitud_ps, latitud_ps, h, eta, representar)

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
%   latitud_ps       - (1x1) latitud del punto subsatélite [deg]
%   longitud_ps      - (1x1) longitud del punto subsatélite [deg]
%   h                - (1x1) altitud de la órbita [m]
%   eta              - (1x1) semiángulo del cono de visibilidad [deg]
%   representar      - (1) representa el swath, (0) no representa el swath
%
%==========================================================================

%Datos de la Tierra
w_T = 7.27*10^-5;               %velocidad angular de la Tierra [rad/s]
r_T = 6371e3;                   %radio Tierra [m]
mu = 398601e9;                  %constante gravitacional de la Tierra[m3/s2]

%Datos de la órbita
a = h + r_T;                    %distancia al centro de la Tierra [m]
T = 2*pi*sqrt(a^3/mu);          %periodo de la órbita [s]
lambda = w_T*2*pi*sqrt(a^3/mu); %longitud entre nodos ascendentes [rad]
lat_ps = deg2rad(latitud_ps);   %latitud del punto subsatélite [rad]
long_ps = deg2rad(longitud_ps); %longitud del punto subsatélite [rad]

%Datos del satélite
eta = deg2rad(eta);             %semiángulo del cono [rad]
r_b = h*tan(eta);               %radio de visión barrido por el satélite [m]
R = h*tan(eta);                 %radio de visibilidad del satélite [m]

%Datos de trigonometría esférica
rho = asin(r_T/a);              %radio angular terrestre [rad]
eps = acos(sin(eta)/sin(rho));  %elevación [rad]
lambda_ = pi/2 - eta - eps;     %radio angular [rad]


%% Representación del swath
% figure
plot_ground_track(rad2deg(lat_ps),rad2deg(long_ps),'Earth')
hold on

p = 0:pi/120:2*pi; %azimut [rad]

%Inicialización de coordenadas del swath
x = zeros(1, length(p));
y = zeros(1, length(p));

%Cálculo y representación del swath
for j = 1:length(p)
    y(j) = asin(cos(lambda_)*sin(lat_ps)+sin(lambda_)*cos(lat_ps)*cos(p(j)));                         %latitud [rad]
    x(j) = acos((cos(lambda_) - sin(y(j))*sin(lat_ps))/(cos(y(j))*cos(lat_ps))) + long_ps;            %longitud [rad]
    x_sym(j) = -acos((cos(lambda_) - sin(y(j))*sin(lat_ps))/(cos(y(j))*cos(lat_ps))) + long_ps;       %longitud simétrica[rad]
    x(j) = rad2deg(x(j));
    y(j) = rad2deg(y(j));
    x_sym(j) = rad2deg(x_sym(j));
end

xx = [x(1:(length(x)-3)/2) x_sym((length(x_sym)-1)/2:length(x_sym))];
xx_ = real(xx);

xx_pos1 = zeros(1,length(xx));
xx_pos2 = zeros(1,length(xx));
xx_neg1 = zeros(1,length(xx));
xx_neg2 = zeros(1,length(xx));

for j = 1:length(p)
%Corrección en caso de que longitud <-180
    if xx(j)>=-180
        xx_pos1(j) = xx(j);
    elseif xx(j)<-180
        xx_pos1(j) = -180;
    end
    if xx(j)>=-180
        xx_pos2(j) = 180;
    elseif xx(j)<-180
        xx_pos2(j) = xx(j)+360;
    end
%Corrección en caso de que longitud >180
     if xx(j)<=180
         xx_neg1(j) = xx(j);
     elseif xx(j)>180
         xx_neg1(j) = 180;
     end
     if xx(j)<=180
         xx_neg2(j) = -180;
     elseif xx(j)>180
         xx_neg2(j) = xx(j)-360;
     end
end

%% Representación del swath

representar = representar;
if representar == 1
    if min(xx_)<-180
        fill(xx_pos1, y,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
        fill(xx_pos2, y,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
        title ('Swath del satélite');
    elseif max(xx_)>180
        fill(xx_neg1, y,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
        fill(xx_neg2, y,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
        title ('Swath del satélite');
     else
         fill(xx, y,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
         title ('Swath del satélite');
    end

    %Cuadrado corrector de los polos
    if lat_ps>0
        if xx_(1) - xx_(length(p)) > 300 | latitud_ps == 90
            x_ = [-180 180 180 -180];
            y_ = [max(y) max(y) 90 90];
            fill(x_,y_,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
        end
    elseif lat_ps<0
        if xx_((length(xx_)-1)/2) - xx_(((length(xx_)+1)/2)) > 300 | latitud_ps == -90
            x_ = [-180 180 180 -180];
            y_ = [-90 -90 min(y) min(y)];
            fill(x_,y_,'b', 'EdgeAlpha', 0,'FaceAlpha', 0.3);
        end
    end
end



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
