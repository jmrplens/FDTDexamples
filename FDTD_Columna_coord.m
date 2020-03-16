%% =========================================================================
%
%   Modelo FDTD de una columna de altavoces direccionable.
%
%   Se incluyen 3 puntos de recepción.
%
%   Jose Manuel Requena Plens (16/03/2020)
%
% =========================================================================
clc, clear, close all;
% Version de la columna donde se maximiza la eficiencia en el punto de
% escucha optimizando inclinacion y curvatura

%% GENERAL
% Caracteristicas del medio
rho = 1.21;     % Densidad del medio (kg/m3)
c = 341;        % Velocidad de propagación (m/s)
k = (c^2)*rho;  % Modulo de compresibilidad
% Frecuencia de emision
freq = 4000;    % (Hz)
% Tipo de excitación ('pulse' o 'CW')
excit = 'CW';
% Dimensiones del mapa en metros
lx = 5;         % m
ly = 3;         % m
% Columna
N = 32;         % Numero de elementos de la columna
d = 0.04;       % Distancia entre elementos (de centro a centro) (m)
Origen = [1,2.5]; % Posicion de origen de la columna (m,m)
% Punto de escucha
Foco = [4,0.5];
% Posicion de los receptores
rec1 = [4,0.25]; % Posicion del receptor 1 (m,m)
rec2 = [1,0.25];% Posicion del receptor 2 (m,m)
rec3 = [lx-1,Origen(2)-(N-1)*d/2]; % Receptor a 0 grados
% Tamaño/Ancho de la PML en metros
lPML = 0.25;    % m
% Parámetros de calculo
dh = .01;       % Definición espacial (m)
dt = dh/341/2;  % Definición temporal (s)
ts = 15;        % Tiempo de simulación (ms)

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

%% POSICIONES DE LOS ELEMENTOS
% Posiciones de los elementos en la matriz del mapa
posNy = round((Origen(2)-(0:N-1) * d)/dh) + nPML;
posNx = round(Origen(1)/dh) + nPML;

%% CÁLCULO DE LOS RETARDOS TEMPORALES (dly)
% Calculo del angulo para emitir hacia el punto de escucha
theta = -rad2deg(atan((Origen(2)-d*N/2-Foco(2)+lPML)/(Foco(1)-Origen(1))));
% Inclinación
x = sin(-theta*pi/180)*d; % Calculo de distancia simulada del segundo elemento
D = x * (0:N-1)'; % Vector de distancias para cada elemento (de arriba a abajo)
dlytilt = D/c; % Vector de tiempos/delays para cada elemento por la inclinacion

% Calculo de la distancia desde el origen de la columna (centro de la
% columna inclinada) hasta el punto de escucha para utilizarlo de radio
% para la focalización
rarc = -sqrt((Foco(1)-(Origen(1)-D(round(N/2))))^2+(Foco(2)-(Origen(2)-d*N/2))^2);
% Arco
thetaarc = linspace(-pi/2-deg2rad(-theta), pi/2-deg2rad(-theta), 1000); % Vector de angulos para una semicircunferencia
Ltot = (N-1)*d; % Longitud total de la columna
% Semicircunferencia
x = (rarc * cos(thetaarc) + Origen(1)-rarc);
y = (rarc * sin(thetaarc) + Origen(2)-Ltot/2);
% Distancias
Darc = zeros(N,1);
for m = 1:N
    [~,idx]=min(abs(y-(Origen(2)-(m-1)*d)));
    Darc(m) = Origen(1)-abs(x(idx));
end
dlyarc = Darc/c; % Vector de tiempos/delays producidos por arquear

% Vector de delays
dly = dlyarc+dlytilt;

% Representacion de las fuentes virtuales frente a la real
% posVirtual = zeros(N,2);
% for m = 1:N
%     posVirtual(m,:) = [Origen(1)-D(m)-Darc(m),Origen(2)-(m-1)*d+lPML];
% end
% figb = figure('Color',[1,1,1]);
% plot(posVirtual(:,1),posVirtual(:,2),'ro', 'MarkerSize', 2,'MarkerFaceColor','r')
% hold on
% plot((ones(size(posNy))*posNx-nPML)*dh,posNy*dh,'ko', 'MarkerSize', 2,'MarkerFaceColor','k')
% axis equal,grid on
% axis([0,lx,0,ly])
% title('Posición de las fuentes')
% xlabel('X [m]'),ylabel('Y [m]')
% legend('Virtual','Real')

%% EXCITACIÓN (la misma señal con diferente retardo para cada elemento)
% El retardo se aplica en el vector 't'
lenT = ts*10^-3/dt; % Longitud del vector de tiempo
aa = freq/(sqrt(pi)/2)*4;
t=((1:lenT)/(1/dt)-4/aa); % Vector de tiempos
w = zeros(lenT,N);
for n = 1:N
    % Funcion que define la excitacion, una para cada n, se incluye el
    % retardo para cada elemento
    switch excit
        case 'pulse'
            w(:,n) = -(exp(-aa^2*((t-dly(n)).^2)/2).*(aa^2*((t-dly(n)).^2)-1)); % Impulso / Ricker con retardo
        case 'CW'
            w(:,n) = cos(2*pi*freq*t)+1i*sin(2*pi*freq*t); % Seno complejo
            % Retardo aplicado
            w(:,n) = w(:,n).*exp(dly(n).*-1j*2*pi*freq);
    end
end

%% PML
% Gradiente de impedancias desde el valor del medio hasta un porcentaje de
% este
gammamax = 0.5; % Maxima reducción de la impedancia / porcentaje
% PML Superior e Inferior
gammaux = zeros(nx+1,ny);
gammaux(1:nPML,:) = repmat(gammamax*((nPML:-1:1)'/nPML).^2,1,ny); % Izquierda
gammaux(1+end-nPML:end,:) = repmat(gammamax*((1:1:nPML)'/nPML).^2,1,ny); % Derecha
gammax = (gammaux(1:end-1,:)+gammaux(2:end,:))/2; % Conjunto
% PML Izquierda y Derecha
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
    px = px.*(1-gammax)-k*dt/dh*diff(ux);
    py = py.*(1-gammay)-k*dt/dh*diff(uy')';
    
    % Excitación
    for n = 1:N
        px(posNx,posNy(n)) = w(tt,n);
        py(posNx,posNy(n)) = w(tt,n);
    end
    p = px+py;
    
    % Velocidad
    ux(2:nx,:) = ux(2:nx,:).*(1-gammaux(2:nx,:))-dt/rho/dh*diff(p);
    uy(:,2:ny) = uy(:,2:ny).*(1-gammauy(:,2:ny))-dt/rho/dh*diff(p')';
    
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
        plot(repmat(posNx-nPML,1,N)*dh,(posNy-nPML)*dh, 'wo', 'MarkerSize', 2,'MarkerFaceColor','w');
        plot([rec1(1),rec2(1),rec3(1)],[rec1(2),rec2(2),rec3(2)], 'bo', 'MarkerSize', 4,'MarkerFaceColor','b');
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
fig= figure('Color',[1,1,1]);
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