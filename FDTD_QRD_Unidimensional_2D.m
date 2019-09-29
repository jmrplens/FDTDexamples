%% =========================================================================
%
%   Modelo FDTD de difusor QRD en recinto 2D con excitación Ricker o
%   continua
%
%
%   Jose Manuel Requena Plens (29/09/2019)
%
% =========================================================================
clc, clear, close all;

%% GENERAL
% Caracteristicas del medio
rho = 1.21;     % Densidad del medio [kg/m3]
c = 341;        % Velocidad de propagación [m/s]
k = (c^2)*rho;  % Modulo de compresibilidad
% Frecuencia de emision
freq = 3000;    % (Hz)
% Tipo de excitación ('pulse' o 'CW')
excit = 'pulse';
% Diseño del difusor QRD
N = 7;
fd = 1000;
fw = 5*fd;
% ¿Difusor o Panel plano?
tipo = 'difusor'; % 'difusor' o 'plano'
% Dimensiones del mapa en metros
lx = 1;         % [m]
ly = 1;         % [m]
% Posicion de la fuente
Origen = [0.5,0.99]; % [m,m]
% Tamaño/Ancho de la PML en metros
lPML = 0.5;    % m
% Parámetros de calculo
dh = .01;       % Definición espacial [m]
dt = dh/341/2;  % Definición temporal [m]
ts = 7;        % Tiempo de simulación [m]

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
a   = freq/(sqrt(pi)/2)*4;
t    = ((1:lenT)/(1/dt)-4/a); % Vector de tiempos
switch excit
    case 'pulse'
        w = -(exp(-a^2*(t.^2)/2).*(a^2*(t.^2)-1)); % Ricker
    case 'CW'
        w = cos(2*pi*freq*t);
end

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

%% Condiciones de contorno
ux(1,:)   = -p(1,:)/rho/c;
ux(end,:) = p(end,:)/rho/c;
uy(:,1)   = -p(:,1)/rho/c;
uy(:,end) = p(:,end)/rho/c;

%% PARÁMETROS DIFUSOR QRD
% Dimensiones
lambdad = c/fd;         % Longitud de onda (profundidad) [m]
lambdaw = c/fw;         % Longitud de onda (ancho) [m]
wide    = lambdaw/2;    % Ancho de cavidad [m]
widedh  = round(wide/dh);
% Vector de profundidades
n0      = ceil(N/2);    % Offset para centrar las secuencias
n       = ceil((1:N))-(n0);
sn      = mod(n.^2,N);
dn      = (sn.*lambdad)/(N*2);      % Profundidades [m]
% Otros
offsetdh = round(wide*(N/2)/dh);    % Para centrar la posición del difusor
maxdh   = round(max(dn)/dh);        % Maxima profundidad
centdh  = round(lx/2/dh);           % Centro del mapa

%% MASCARA DEL DIFUSOR
mascdifusor = zeros(nx,ny); % Inicialización de la mascara
for n = 1:N+1
    % Lineas verticales
    mascdifusor( (widedh*(n-1))+1+centdh-offsetdh+nPML,(1:maxdh)+nPML ) = 1;
    % Líneas horizontales
    if n<=N
        switch tipo
            case 'plano'
                dnh = 0; % Profundidad de la cavidad n
            case 'difusor'
                dnh = round(dn(n)/dh); % Profundidad de la cavidad n
        end
        xhor = (((widedh*(n-1)):(widedh*(n)) )+1+centdh-offsetdh)+nPML;
        yhor = maxdh-dnh+nPML;
        mascdifusor( xhor,yhor ) = 1;
    end
end
% Línea horizontal en la parte inferior del difusor
mascdifusor( ((0:(widedh*N) )+1+centdh-offsetdh)+nPML,nPML ) = 1;
% Diferenciación de la mascara
mascarax=1-(diff(mascdifusor).'~=0).';
mascaray=1-(diff(mascdifusor.')~=0).';
% Inversión de la mascara
mascdifusor = 1-mascdifusor;

%% REPRESENTACION MASCARAS Y PML
f = figure('Color',[1,1,1]);
subplot(1,3,1),pcolor(((1:nx)-nPML)*dh,((1:ny)-nPML)*dh,mascdifusor'), hold on
plot([0,0,nx-2*nPML,nx-2*nPML,0]*dh,[0,ny-2*nPML,ny-2*nPML,0,0]*dh,'r--','linewidth',2)
shading flat, axis equal,axis([-nPML,nx-nPML,-nPML,ny-nPML]*dh),title('Máscara')
subplot(1,3,2),pcolor(((1:nx-1)-nPML)*dh,((1:ny)-nPML)*dh,mascarax'), hold on
plot([0,0,nx-2*nPML,nx-2*nPML,0]*dh,[0,ny-2*nPML,ny-2*nPML,0,0]*dh,'r--','linewidth',2)
shading flat, axis equal,axis([-nPML,nx-nPML,-nPML,ny-nPML]*dh),title('Máscara en x')
subplot(1,3,3),pcolor(((1:nx)-nPML)*dh,((1:ny-1)-nPML)*dh,mascaray'), hold on
plot([0,0,nx-2*nPML,nx-2*nPML,0]*dh,[0,ny-2*nPML,ny-2*nPML,0,0]*dh,'r--','linewidth',2)
shading flat, axis equal,axis([-nPML,nx-nPML,-nPML,ny-nPML]*dh),title('Máscara en y')
fprintf('pulsa cualquier tecla para continuar');
pause
close(f)

%% Calculo
fig1 = figure('Color',[1,1,1]);
axis equal;axis([0,nx-nPML*2,0,ny-nPML*2]*dh)
xlabel('X [m]'),ylabel('Y [m]')
box on
colormap(flipud(bone(256)))
for tt=1:lenT
    % Presion
    px = px.*(1-gammax)-k*dt/dh*diff(ux).*mascdifusor;
    py = py.*(1-gammay)-k*dt/dh*diff(uy')'.*mascdifusor;
    
    % Excitación
    px(posNx,posNy) = w(tt);
    py(posNx,posNy) = w(tt);
    p = px+py;
    
    % Velocidad
    ux(2:nx,:) = ux(2:nx,:).*(1-gammaux(2:nx,:))-dt/rho/dh*diff(p).*mascarax;
    uy(:,2:ny) = uy(:,2:ny).*(1-gammauy(:,2:ny))-dt/rho/dh*diff(p')'.*mascaray;
    
    % Representación gráfica
    if tt/10==round(tt/10)
        map = abs(p);
        surface(((1:nx)-nPML)*dh,((1:ny)-nPML)*dh,abs(mascdifusor-1)'*5)
        hold on
        pcolor(((1:nx)-nPML)*dh,((1:ny)-nPML)*dh,map'), hold off
        set(gca,'Clim',[0 0.5]);shading interp
        title(['Tiempo = ' num2str(round((tt)*1000*dt)) ' ms']);
        drawnow
    end
end