%% =========================================================================
%
%   Modelo simple de FDTD de un recinto en 2D con excitación Ricker y PML
%
%   Se incluyen 3 puntos de recepción y un ejemplo básico para introducir
%   objetos rigidos que interactuan con el campo acústico.
%
%   Jose Manuel Requena Plens (09/09/2019)
%
% =========================================================================
clc, clear, close all;

%% GENERAL
% Caracteristicas del medio
rho = 1.21;     % Densidad del medio [kg/m3]
c = 341;        % Velocidad de propagación [m/s]
K = (c^2)*rho;  % Modulo de compresibilidad
% Dimensiones del mapa en metros
lx = 5;         % Ancho [m]
ly = 3;         % Alto [m]
%
Origen = [1,2.5]; % Origen de la fuente
% Posicion de los receptores
rec1 = [4,0.25]; % Posicion del receptor 1 [m,m]
rec2 = [1,0.25];% Posicion del receptor 2 [m,m]
rec3 = [2,2]; % Posicion del receptor 3 [m,m]
% Tamaño/Ancho de la PML en metros
lPML = 0.25;    % m
% Parámetros de calculo
dh = .01;       % Definición espacial [m]
dt = dh/341/2;  % Definición temporal [s]
ts = 15;        % Tiempo de simulación [ms]

%% MATRICES
% Dimensiones de la matriz del mapa
nPML = round(lPML/dh);
nx   = round(lx/dh)+nPML*2; % Por 2 ya que la PML esta arriba y abajo
ny   = round(ly/dh)+nPML*2; % Por 2 ya que la PML esta a izquierda y derecha
% Crear matriz del mapa
p  = zeros(nx,ny);   % Matrices de presion
px = zeros(nx,ny);
py = zeros(nx,ny);
ux = zeros(nx+1,ny); % Matriz de velocidad de particula en x
uy = zeros(nx,ny+1); % Matriz de velocidad de particula en y

%% POSICION DE LA FUENTE
posNy = round((Origen(2))/dh) + nPML;
posNx = round(Origen(1)/dh) + nPML;


%% EXCITACIÓN
lenT = ts*10^-3/dt; % Longitud del vector de tiempo
a    = 2000/(sqrt(pi)/2)*4;
t    = ((1:lenT)/(1/dt)-4/a); % Vector de tiempos
w    = -(exp(-a^2*(t.^2)/2).*(a^2*(t.^2)-1)); % Ricker

%% PML
% Gradiente de impedancias desde el valor del medio hasta un porcentaje de
% este
gammamax = 0.5; % Maxima reducción de la impedancia / porcentaje
% PML Izquierda y derecha
gammaux = zeros(nx+1,ny);
gammaux(1:nPML,:) = repmat(gammamax*((nPML:-1:1)'/nPML).^2,1,ny); % Izquierda
gammaux(1+end-nPML:end,:) = repmat(gammamax*((1:1:nPML)'/nPML).^2,1,ny); % Derecha
gammax = (gammaux(1:end-1,:)+gammaux(2:end,:))/2; % Conjunto
% PML Superior e inferior
gammauy = zeros(nx,ny+1);
gammauy(:,1:nPML) = repmat(gammamax*((nPML:-1:1)/nPML).^2,nx,1); % Inferior
gammauy(:,1+end-nPML:end) = repmat(gammamax*((1:1:nPML)/nPML).^2,nx,1); % Superior
gammay = (gammauy(:,1:end-1)+gammauy(:,2:end))/2; % Conjunto

%% Inicializacion de los vectores de los puntos de escucha
h = zeros(1,lenT);
h2 = zeros(1,lenT);
h3 = zeros(1,lenT);

%% Condiciones de contorno
ux(1,:)   = -p(1,:)/rho/c;
ux(end,:) = p(end,:)/rho/c;
uy(:,1)   = -p(:,1)/rho/c;
uy(:,end) = p(:,end)/rho/c;

%% Calculo
fig1 = figure('Color',[1,1,1]);

colormap(jet(256))
for tt=1:lenT
    % Presion
    px = px.*(1-gammax)-K*dt/dh*diff(ux);
    py = py.*(1-gammay)-K*dt/dh*diff(uy')';
    
    % Excitación
        px(posNx,posNy) = w(tt);
        py(posNx,posNy) = w(tt);
    p = px+py;
    
    % Velocidad
    ux(2:nx,:) = ux(2:nx,:).*(1-gammaux(2:nx,:))-dt/rho/dh*diff(p);
    uy(:,2:ny) = uy(:,2:ny).*(1-gammauy(:,2:ny))-dt/rho/dh*diff(p')';
    
    % Objetos (Se puede eliminar o modificar sin introducir valores fuera
    % de las dimensiones del mapa)
    % Barrera de 30cm de ancho
    orX = 1; % Origen en X
    orY = 0; % Origen en Y
    wid = 0.3; % Ancho de la barrera
    hei = 1.5; % Altura de la barrera
    vecX = (round(orX/dh):round((orX+wid)/dh))+nPML;
    vecY = (round(orY/dh)+1:round(hei/dh))+nPML;
    ux(vecX,vecY) = 0;
    uy(vecX,vecY) = 0;
    
    % Respuesta al impulso
    h(tt) = p(round(rec1(1)/dh)+nPML,round(rec1(2)/dh)+nPML);
    h2(tt) = p(round(rec2(1)/dh)+nPML,round(rec2(2)/dh)+nPML);
    h3(tt) = p(round(rec3(1)/dh)+nPML,round(rec3(2)/dh)+nPML);
    
    % Representación gráfica
    if tt/10==round(tt/10)
        splmap = 10*log10(abs(p).^2/(2e-5)^2);
        %splmap = splmap - max(splmap(:));
        pcolor(((1:nx)-nPML)*dh,((1:ny)-nPML)*dh,splmap');
        hold on;
        plot([0,0,nx-2*nPML,nx-2*nPML,0]*dh,[0,ny-2*nPML,ny-2*nPML,0,0]*dh,'r--','linewidth',2)
        plot([rec1(1),rec2(1),rec3(1)],[rec1(2),rec2(2),rec3(2)], 'ro', 'MarkerSize', 2,'MarkerFaceColor','r');
        text([rec1(1),rec2(1),rec3(1)]+0.05,[rec1(2),rec2(2),rec3(2)],{'Rec1','Rec2','Rec3'},'Color','white')
        hold off
        axis equal;axis([-nPML,nx-nPML,-nPML,ny-nPML]*dh)
        set(gca,'Clim',[50 110]);shading flat,cc = colorbar;
        title(['Tiempo = ' num2str(round((tt)*1000*dt)) ' ms']);
        xlabel('X [m]'),ylabel('Y [m]'),cc.Label.String = 'SPL (dB)';
        drawnow
    end
end

% Respuesta al impulso / Punto de escucha
time = seconds(t);
maxx = max([abs(h),abs(h2),abs(h3)]);
h = h./maxx; h2 = h2./maxx; h3 = h3./maxx;
T = array2timetable([h',h2',h3'],...% Datos iniciales
'RowTimes',time); % Columna de tiempo
fig= figure;
s = stackedplot(T,...
    'Title','Señal recibida',...
    'DisplayLabels',{'Receptor 1', 'Receptor 2', 'Receptor 3'},...
    'XLabel','Señal recibida');
s.GridVisible = 1;
s.AxesProperties(1).YLimits = [-1.2,1.2];
s.AxesProperties(2).YLimits = [-1.2,1.2];
s.AxesProperties(3).YLimits = [-1.2,1.2];
fig.Children.FontName = 'Arial';
fig.Children.FontSize = 12;