
clear
clc
close all  

% VENTANAS NÚMERICAS
N = 2^9 ;                     % Número de paso
w0 = .5E-3;                    % Cintura del haz
lambda = 633E-9 ; % longitud de onda
k0= 2*pi/lambda; % Numero de onda
L = 4 * w0;                  % La ventana o el tamaño de la malla donde se vera el haz esta en función de la cintura
dR = (pi*w0^2)/lambda  ;      % Distancia Reyleigh
Lz = 1*dR;                    % Distancia en Z
dx = 2 * L/N  ;               % Nyquist 
kmax = pi/dx;                 % Paso a mi espacio de frecuencias 
dz = Lz/N ;                   % Tamaño de los pasos
dKx= 2*kmax/N;

% CREACIÓN DE COORDENADAS 

% EN ESPACIO REAL
x = (-N/2:1:N/2-1) * dx   ;                 % Cuando multiplicas NV * dx, nmormalizas el vector NV y le das el tamaño de tu ventana de trabajo
y = (-N/2:1:N/2-1) * dx  ;                  % Esto lo hacemos para que nuestros vectores de x y y 
z = 0:dz:Lz;
[X, Y] = meshgrid(x,y);      % Matrices con coordenadas por aparte 
rho2 = X.^2 + Y.^2;

%% Paso 1: Campo GB Normalizado
Um0 = exp(-rho2/w0^2);
%Cte_Norm_Um0 = sqrt(1/(trapz(y,trapz(x,abs(Um0).^2,2),1))); % Ctw de Normalizacion
Cte_Norm_Um0 = sqrt(1/(trapz(y,trapz(x,Um0.*conj(Um0),2),1))); % Ctw de Normalizacion

Um0 = Um0.*Cte_Norm_Um0 ; % Campo Um0 Normalizado
IUm0 = trapz(y,(trapz(x,abs(Um0).^2))); % Checar normalizacion campo Um0

%% Paso 2: Propagar GB y ponerle fase

f = .3*dR;  % Distancia Real
ttz = f/dz; % Pasos
Fase = exp(1i*atan2(Y,X)); % Fase que le meto
Nz = N;

Prpg = PropagadorFourier(Nz,Lz,dx,Um0,k0,Fase,ttz,1);
Um1 =  Prpg(:,:,ceil(ttz+f/dz)); % Campo despues de propagarse un foco de distancia 
% IUm1 = trapz(y,(trapz(x,abs(Um1).^2))); % Checar normalizacion campo Um1


figure(2)
uu = Prpg(end/2,:,:); % Grafica del propagador con un corte transverzal en y. Vemos X contra Z 
Ures = reshape(uu(:),size(uu,2),size(uu,3));
imagesc(z,y,abs(Ures));colorbar; colormap(hot); 
line([z(ceil(ttz)),z(ceil(ttz))], [y(1),y(end)], 'Color', 'r')
line([z(ceil(ttz))+f,z(ceil(ttz))+f], [y(1),y(end)], 'Color', 'b')
title('Propagación en Z')

figure(3)
imagesc(x,y,abs(Um0).^2), axis equal; colormap(hot); colorbar
title('Vista transversal del Campo INICIAL')
figure(4)
imagesc(x,y,angle(Um0)), axis equal; 
title('Fase del Campo Inicial')

figure(5)
imagesc(x,y,abs(Prpg(:,:,ceil(ttz+1))).^2), axis equal; colormap(hot); colorbar
title('Vista transversal del Campo tras haber interactuado con la fase INSTANTANEAMENTE')
figure(6)
imagesc(x,y,angle(Prpg(:,:,ceil(ttz+1)))), axis equal; colormap(hot); colorbar
title('Fase del Campo tras haber interactuado con la fase INSTANTANEAMENTE')

figure(7)
imagesc(x,y,abs(Um1).^2), axis equal; colormap(hot); colorbar
title('Vista transversal del Campo tras haber interactuado con la fase UN FOCO DESPUÉS')
figure(8)
imagesc(x,y,angle(Um1(1:2:end,1:2:end))), axis equal; colorbar
title('Fase del Campo tras haber interactuado con la fase  UN FOCO DESPUÉS')


%% Paso 3: Comparacion con LG10
% 
% HOL = conj(LaguerreGauss(X,Y,1,0,w0));
% Lente = ThinLense(X,Y,f,k0);
% 
% Hol_Um1 = (Um1.*HOL).*sqrt(1/(trapz(y,trapz(x,abs(Um1.*HOL).^2,2),1)));
% 
% Prpg2 = PropagadorFourier(Nz,Lz,dx,(Hol_Um1),k0,Lente,ttz,0);
% 
% Um2 =  Prpg2(:,:,ceil(ttz)); % Campo despues d einteractuar con el HOL
% Wght = abs(Um2(N/2+1,N/2+1)).^2;
% 
% figure(4)
% uu = Prpg2(end/2,:,:); % Grafica del propagador con un corte transverzal en y. Vemos X contra Z 
% Ures = reshape(uu(:),size(uu,2),size(uu,3));
% imagesc(z,y,abs(Ures).^0.5);colorbar; colormap(hot); 
% line([z(ceil(ttz)),z(ceil(ttz))], [y(1),y(end)], 'Color', 'r')
% %line([z(ceil(ttz))+f,z(ceil(ttz))+f], [y(1),y(end)], 'Color', 'b')
% title('Propagación en Z')
% 
% figure(5)
% imagesc(x,y,abs(Um2).^2), axis equal; colormap(hot); colorbar
% title('Vista transversal del Campo tras un foco de haber interactuado con la fase')

%% Comparacion con LG con indices variantes

Cc = 6;
m = -Cc:1:Cc;
p = 0:1:size(m,2)-1;

Wght = zeros(size(p,2),size(m,2));

%CASO PRUEBA 
%Um1 = (LaguerreGauss(X,Y,2,1,w0)+0.5*LaguerreGauss(X,Y,3,1,w0));
%Um1 = (Um1).*sqrt(1/(trapz(y,trapz(x,abs(Um1).^2,2),1)));
Norm = trapz(y,trapz(x,abs(Um1).^2,2),1);

for Mm = 1:size(m,2)
    for Pp = 1:size(p,2)

        HL = conj(LaguerreGauss(X,Y,m(Mm),p(Pp),w0));
        Um2 =  Um1.*HL; % Campo interactuando con el Holograma
        Hol_Um1 = (Um2).*sqrt(1/(trapz(y,trapz(x,abs(Um2).^2,2),1)));
        HolPrpg = trapz(y,trapz(x,Hol_Um1,2),1); % Propagacion pasando por el Holograma
        Wght(Pp,Mm) = (HolPrpg);
    end
end


Wght =((Wght.*(1./sqrt(sum(sum(abs(Wght).^2)))))); % Normalizacion de los pesos

% Reconstruccción del Haz
Haz = 0;
for Mm = 1:size(m,2)
    for Pp = 1:size(p,2)
        Peso = Wght(Pp,Mm).*(LaguerreGauss(X,Y,m(Mm),p(Pp),w0));
        Haz = Peso + Haz;
    end
end 

figure(2)
imagesc(x,y,abs(Haz).^2); colormap(hot);
title('Reconstrucción')
figure(4)
imagesc(x,y,angle(Haz)) ;colormap(hot);
title('Reconstrucción')
figure(3)
imagesc(x,y,abs(Um1).^2) ;colormap(hot);
title('Original')
figure(5)
imagesc(x,y,angle(Um1)) ;colormap(hot);
title('Original')
figure (6)
stem3(m,p,abs(Wght).^2)
title('Descomposición modal de LGmp con m = [-3,3] y p = 0')

%%
surf(m,p,Wght)
Xx = m;
Yy = p;
Zz = zeros(size(Xx,2));
U = zeros(size(Xx,2));
V = zeros(size(Xx,2));
W =Wght;


quiver3(Xx,Yy,Zz,U,V,W); axis equal
%% Parte 5, Reconstruccion del Haz Inicial
