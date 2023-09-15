clc;
clear; 
tic
%% Constantes: 

L_target = 1; %cm
dz = 0.1; % cm
N_laminas = L_target/dz; 
realE = 230; % MeV
sigmaE = 0.1; % MeV
Rad_length_Pb = 0.5612; % cm 
mc2_p = 938.27208816; % MeV/c²
tabla = importdata('StoppingPowerTable');   
Densidad_Pb = 11.34; %g/cm3
% Calculo la perdida de energia con tabla (NIST): 
V_energia = tabla(:,1); %MeV
V_perdida = tabla(:,2); %MeV*cm2/g
V_perdida = V_perdida * Densidad_Pb; %MeV/cm

%% Elipse: 
radio_corto = 1.6; %Optimizar
%Parametrizo la elipse: 
L_target2 = 1.16; 
y = (L_target2-radio_corto:0.01:L_target2); % <- Cambio aquí
x1 = real(sqrt(3^2*(1-(y-L_target2).^2/radio_corto^2))); 
x2 = real(-sqrt(3^2*(1-(y-L_target2).^2/radio_corto^2))); 
xt = [flip(x2),x1]; % <- Cambio aquí
xt = round(xt,4);
yt = [flip(y),y]; 
%plot(xt,yt);

%% PHSP de entrada experimental: 
%inPHSP1 = readtable('MapaDosisExpMasParticulas.txt'); %No deja subirlo a
%github por su tamaño
%inPHSP1 = readtable('MapaDosisExp.txt');
inPHSP1 = readtable('MapaDosisExp.txt');
inPHSP1 = table2array(inPHSP1); 
inPHSP1(:,[3,7,8,9,10]) = []; 
Npart = size(inPHSP1,1);
%% PHSP de entrada Haz Puntual: 
% Npart = 200000;
% sigma0_cm = 0.001;
% inPHSP1 = zeros(Npart,5);
% inPHSP1(:,1:2) = normrnd(0,sigma0_cm,Npart, 2); % ¿Porque de 0 a 1?
% inPHSP1(:,3:4) = 0; % Haz perpendicular al target
% inPHSP1(:,5) = normrnd(realE,sigmaE,Npart,1);
%% PHSP de entrada Haz Gaussiano: 
% sigma0_cm = 0.5;
% inPHSP1 = zeros(Npart,5);
% inPHSP1(:,1:2) = normrnd(0,sigma0_cm,Npart, 2); % ¿Porque de 0 a 1?
% inPHSP1(:,3:4) = 0; % Haz perpendicular al target
% inPHSP1(:,5) = normrnd(realE,sigmaE,Npart,1);
%% PHSP de salida del box (primera parte del colimador): 

inPHSP = inPHSP1; 
outPHSP = zeros(size(inPHSP));

for i = 1:N_laminas

[vectorEnergia] = CalculateNewEnergy_V2 (dz,V_energia, V_perdida, inPHSP); 
outPHSP(:,5) = vectorEnergia; 

% 2nd calculate new sigma: 
meanE = mean(inPHSP(:,5));
tao = meanE/mc2_p; % adimensional
pv = (tao+2)/(tao+1)*meanE; % MeV
sigmaTheta = (14.1/pv)*sqrt(dz/Rad_length_Pb)*(1+(1/9)*log10(dz/Rad_length_Pb)); % rad

% 2nd calculate new angles
outPHSP(:,3:4) = inPHSP(:,3:4) + normrnd(0,sigmaTheta,Npart,2);

% 3rd calculate new positions using new angles  %% 
outPHSP(:,1) = inPHSP(:,1) + dz.*tan(outPHSP(:,3));  
outPHSP(:,2) = inPHSP(:,2) + dz.*tan(outPHSP(:,4));

inPHSP = outPHSP;

end

%% PHSP salida tras curvatura (segunda parte del colimador): 

dz2 = 0.1; %cm
N_laminas2 = (max(yt)-L_target)/dz2; 
onesIdx = (1:length(inPHSP(:,1))); %Los que están dentro del poligono
ncIdx = []; %Los que están fuera del poligono
%Creamos planos rectangulares mediante vertices: 
xv = zeros(1,4); 
zv = zeros(1,4); 
%Para proyectar en aire:
finalPos = zeros(Npart,2);

for i = 1:N_laminas2

[vectorEnergia] = CalculateNewEnergy_V3 (dz2,V_energia, V_perdida, inPHSP,onesIdx); 
outPHSP(onesIdx,5) = vectorEnergia; 

% 2nd calculate new sigma: 
meanE = mean(inPHSP(onesIdx,5));
tao = meanE/mc2_p; % adimensional
pv = (tao+2)/(tao+1)*meanE; % MeV
sigmaTheta = (14.1/pv)*sqrt(dz2/Rad_length_Pb)*(1+(1/9)*log10(dz2/Rad_length_Pb)); % rad

% 2nd calculate new angles
Npart = length(onesIdx);
outPHSP(onesIdx,3:4) = inPHSP(onesIdx,3:4) + normrnd(0,sigmaTheta,Npart,2);

% 3rd calculate new positions using new angles  
outPHSP(onesIdx,1) = inPHSP(onesIdx,1) + dz2.*tan(outPHSP(onesIdx,3));  
outPHSP(onesIdx,2) = inPHSP(onesIdx,2) + dz2.*tan(outPHSP(onesIdx,4));

inPHSP = outPHSP;

%Cojo los indices de las partículas que siguen dentro del Pb: 
xq = inPHSP(:,1);  
yq = inPHSP(:,2);
%Para que no haya problemas con los ya proyectados: 
xq(ncIdx) = -10; 
yq(ncIdx) = -10;
%Calculo el plano que hay que coger para generar nuestro poligono
%(rectangulo):  
Dist_plano = L_target+i*dz2; 
%Busco los indices de los más cercanos:  
%1º Me quedo solo con la primera mitad del vector yt:
yap = yt(1:length(yt)/2);
closest = interp1(yap,yap,Dist_plano,'nearest');
idx = zeros(1,2);
idx(1) = find(yap==closest);
%Busco el simétrico:
idx(2) = length(yap) + (length(yap)-idx(1)) + 1;
%%
%Creo los vertices del poligono: 
xv(1) = xt(idx(1)); 
xv(2) = xt(idx(1));
xv(3) = xt(idx(2)); 
xv(4) = xt(idx(2));
zv(1) = -1.5;
zv(2) = 1.5; 
zv(3) = 1.5; 
zv(4) = -1.5; 

%Veo que cordenadas están dentro del polígono: 
in = inpolygon(xq,yq,xv,zv);
%Busco los indices de los ceros (fuera) y unos (dentro): 
%En este caso, es al reves lo que se encuentren dentro ya no se dispersan
%más: 
onesIdx = find(in==0);
ncIdx = find(in==1); 

%Los que den 0 proyectamos ya al aire: 
dzAire_cm = 4 - Dist_plano; % A 4 cm
finalPos(ncIdx,1) = outPHSP(ncIdx,1) + dzAire_cm.*tan(outPHSP(ncIdx,3));  
finalPos(ncIdx,2) = outPHSP(ncIdx,2) + dzAire_cm.*tan(outPHSP(ncIdx,4));

end
dzAire_cm = 4 - L_target+radio_corto; 
finalPos(onesIdx,1) = outPHSP(onesIdx,1) + dzAire_cm.*tan(outPHSP(onesIdx,3));
finalPos(onesIdx,2) = outPHSP(onesIdx,2) + dzAire_cm.*tan(outPHSP(onesIdx,4));
%% Plots: 
%Mapas de Dosis
paso = 0.01693; %cm
%paso = 0.03;
Xedges = (-1.19357:paso:1.19357); %cm
Yedges = (-0.73646:paso:0.73646); %cm
figure(1)
subplot(3,1,1)
h1 = histogram2(inPHSP1(:,1),inPHSP1(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
clear title
title('Mapa PHSP entrada')
subplot(3,1,2)
h2 = histogram2(outPHSP(:,1),outPHSP(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
clear title
title('Mapa PHSP salida'); 
subplot(3,1,3)
h3 = histogram2(finalPos(:,1),finalPos(:,2),Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','off');
clear title
title('Mapa PHSP salida + 1 cm aire'); 

%Saco las matrices
Matrix_h1 = h1.Values;
Matrix_h2 = h2.Values;
Matrix_h3 = h3.Values;

%Aplico una rotación y un flip, para dejarlas igual que las experimentales:
Matrix_h1 = rot90(Matrix_h1); 
Matrix_h1 = flip(Matrix_h1,1); 
Matrix_h2 = rot90(Matrix_h2); 
Matrix_h2 = flip(Matrix_h2,1); 
Matrix_h3 = rot90(Matrix_h3); 
Matrix_h3 = flip(Matrix_h3,1); 

%% Analisis: 
%Cargo experimentales: 
load('processedDoseMaps.mat');
IMap1 = doseMaps{1};
%cGy a Gy:
IMap1 = IMap1./100; 

%Mapas Matlab:
%Paso de fluencia a Dosis: Dosis = 0.1602*fluencia*Mass Stopping Power [Gy] 
%Fluencia: 
Area_bin = paso^2; %cm2
Matrix_h1 = Matrix_h1/Area_bin;
Matrix_h3 = Matrix_h3/Area_bin; % En fluencia
%Stopping Power:
E_mean_h1 = mean(inPHSP1(:,5));
E_mean_h3 = mean(outPHSP(:,5)); 
SP_h1 = interp1(V_energia,V_perdida,E_mean_h1);
SP_h3 = interp1(V_energia,V_perdida,E_mean_h3);
%Matriz en Dosis: 
Matrix_h1 = Matrix_h1*SP_h1*0.1602; %Gy
Matrix_h3 = Matrix_h3*SP_h3*0.1602; %Gy

%Factor multiplicación: No se necesita el paso a Dosis de arriba entonces? Da lo mismo lo incluya o no 
M_IMap1 = mean2(IMap1);
M_h1 = mean2(Matrix_h1);
Factor_h = M_IMap1/M_h1;

Matrix_h1 = Matrix_h1*Factor_h; 
Matrix_h3 = Matrix_h3*Factor_h; 

M1 = Matrix_h1; 
M3 = Matrix_h3; 

%Ahora ya tengo las matrices con la dosis correcta: %
%% Datos area completa: 
thresholdBajo = 10; 
threshold_h1 = (thresholdBajo / 100) * max(max(M1)); 
M1(M1<threshold_h1) = 0; 
threshold_h = (thresholdBajo / 100) * max(max(M3)); 
M3(M3<threshold_h) = 0; 

mascara1 = (M1 ~= 0);
mascara2 = (M3 ~= 0);

Media1 = mean(M1(mascara1));
Media2 = mean(M3(mascara2));

StdDesv1 = std(M1(mascara1)); 
StdDesv2 = std(M3(mascara2)); 

Area1 = sum(M1(mascara1))*paso*paso; %cm
Area2 = sum(M3(mascara2))*paso*paso; %cm

%% Datos Iso80:

AA_h1 = prctile(M1(mascara1),98,'all'); 
AA_h3 = prctile(M3(mascara2),98,'all');
 
M1A = M1/AA_h1; 
M3A = M3/AA_h3; 

M1(M1A<0.8) = 0; 
M3(M3A<0.8) = 0; 

mascara1Iso = (M1 ~= 0);
mascara2Iso = (M3 ~= 0);


Media80M1 = mean(M1(mascara1Iso));
Media80M3 = mean(M3(mascara2Iso)); 

Std80M1 = std(M1(mascara1Iso));
Std80M3 = std(M3(mascara2Iso));

%Area: 
A1 = sum(M1(mascara1Iso))*paso*paso; %cm
A2 = sum(M3(mascara2Iso))*paso*paso; %cm
 
%Guardo en tabla:
Resul = zeros(2,6); 
Resul(1,1) = StdDesv1; 
Resul(2,1) = StdDesv2; 
Resul(1,2) = Area1; 
Resul(2,2) = Area2; 
Resul(1,3) = Media1; 
Resul(2,3) = Media2; 
Resul(1,4) = Std80M1; 
Resul(2,4) = Std80M3; 
Resul(1,5) = A1; 
Resul(2,5) = A2; 
Resul(1,6) = Media80M1; 
Resul(2,6) = Media80M3; 
Resul = array2table(Resul); 
ResultadosIso80 = renamevars(Resul, ["Resul1","Resul2","Resul3","Resul4","Resul5","Resul6"],["StdTot [Gy]","AreaTot [cm2]","DosisMediaTot [Gy]","StdIso80 [Gy]","AreaIso80 [cm2]","DosisMediaIso80 [Gy]"]);
% Nombres de filas
rowNames = {'Previo Colimador', 'Plano de Pelets'};
% Asignar nombres a las filas
ResultadosIso80.Properties.RowNames = rowNames;

%% Datos Pelet: 
% Creo ROI: 
% %Las células por el método de cultivo (se les da vueltas) quedan
% %depositadas en forma de triangulo: 
% xi = zeros(1,3); 
% yi = zeros(1,3);
% Xedges(end) = []; %cm
% Yedges(end) = []; 
% Px1 = [-0.1 0.12 -0.05]; %cm
% Py1 = [-0.1 -0.14 0.15]; 
% closest = interp1(Xedges,Xedges,Px1(1),'nearest');
% xi(1) = find(Xedges==closest);
% closest = interp1(Xedges,Xedges,Px1(2),'nearest');
% xi(2) = find(Xedges==closest);
% closest = interp1(Xedges,Xedges,Px1(3),'nearest');
% xi(3) = find(Xedges==closest);
% closest = interp1(Yedges,Yedges,Py1(1),'nearest');
% yi(1) = find(Yedges==closest);
% yi(1) = length(Yedges)-yi(1);
% closest = interp1(Yedges,Yedges,Py1(2),'nearest');
% yi(2) = find(Yedges==closest);
% yi(2) = length(Yedges)-yi(2);
% closest = interp1(Yedges,Yedges,Py1(3),'nearest');
% yi(3) = find(Yedges==closest);
% yi(3) = length(Yedges)-yi(3);
% BW = poly2mask(xi,yi,length(Yedges),length(Xedges)); 
% %imshow(BW)
% %Eppendorf lado izquierdo: 
% xi2 = zeros(1,3); 
% yi2 = zeros(1,3); 
% Px2 = [-0.2 -0.42 -0.25]; 
% Py2 = [-0.1 -0.14 0.15]; 
% closest = interp1(Xedges,Xedges,Px2(1),'nearest');
% xi2(1) = find(Xedges==closest);
% closest = interp1(Xedges,Xedges,Px2(2),'nearest');
% xi2(2) = find(Xedges==closest);
% closest = interp1(Xedges,Xedges,Px2(3),'nearest');
% xi2(3) = find(Xedges==closest);
% closest = interp1(Yedges,Yedges,Py2(1),'nearest');
% yi2(1) = find(Yedges==closest);
% yi2(1) = length(Yedges)-yi2(1);
% closest = interp1(Yedges,Yedges,Py2(2),'nearest');
% yi2(2) = find(Yedges==closest);
% yi2(2) = length(Yedges)-yi2(2);
% closest = interp1(Yedges,Yedges,Py2(3),'nearest');
% yi2(3) = find(Yedges==closest);
% yi2(3) = length(Yedges)-yi2(3);
% BW2 = poly2mask(xi2,yi2,length(Yedges),length(Xedges)); 
% %imshow(BW2)

% Dimensiones de la imagen
Xedges(end) = []; %cm
Yedges(end) = []; 
height = length(Yedges);
width = length(Xedges);
% Tamaño del píxel en cm
pixelSize = 0.01693;
% Parámetros de la media luna convertidos a cm
radio_cm = 0.573 / 2 ; % cm
EspacioEntrePocillos_cm = 0.1; % cm
% Ahora mapeamos estas dimensiones en cm a píxeles
radio_pix = radio_cm / pixelSize;
EspacioEntrePocillos_pix = EspacioEntrePocillos_cm / pixelSize;
center_x = width/2 + radio_pix + EspacioEntrePocillos_pix/2;
center_y = height/2;
% Crea una matriz de coordenadas
[X, Y] = meshgrid(1:width, 1:height);
% Ecuación del círculo (solo media luna)
BW = (X - center_x).^2 + (Y - center_y).^2 <= radio_pix^2 & Y >= center_y;
imshow(BW);
% Pocillo izquierdo:
center_x2 = width/2 - radio_pix - EspacioEntrePocillos_pix/2;
BW2 = (X - center_x2).^2 + (Y - center_y).^2 <= radio_pix^2 & Y >= center_y;
imshow(BW2);
imshow(BW + BW2);

%% Histograma de dosis dentro de la zona: 
%Pocillo derecho: 
DosisPocillo1 = Matrix_h3(BW);
figure(2);
histogram(DosisPocillo1,80);
xlabel('Dosis Gy');
title('Histograma Pocillo derecho');
%Pocillo izquierdo: 
DosisPocillo2 = Matrix_h3(BW2);
figure(3);
histogram(DosisPocillo2,80);
xlabel('Dosis Gy');
title('Histograma Pocillo izquierdo');

% Dosis en pelets:
DosisPocilloDerecho = mean(Matrix_h3(BW));
DosisPocilloIzquierdo = mean(Matrix_h3(BW2));

% Tasa de Dosis en Pelets:
tIrradiation_s = 0.01; %s
TasaDosisPocilloDerecho = DosisPocilloDerecho/tIrradiation_s;
TasaDosisPocilloIzquierdo = DosisPocilloIzquierdo/tIrradiation_s;

% Desviacion estandard en Pelets:
StdPocilloDerecho = std(Matrix_h3(BW)); 
StdPocilloIzquierdo = std(Matrix_h3(BW2)); 

%Ponemos resultados en tabla: 
ResulPelet = zeros(2,3); 
ResulPelet(1,1) = DosisPocilloDerecho; 
ResulPelet(2,1) = DosisPocilloIzquierdo; 
ResulPelet(1,2) = StdPocilloDerecho; 
ResulPelet(2,2) = StdPocilloIzquierdo; 
ResulPelet(1,3) = TasaDosisPocilloDerecho;
ResulPelet(2,3) = TasaDosisPocilloIzquierdo;
ResulPelet = array2table(ResulPelet); 
ResultadosPelets = renamevars(ResulPelet, ["ResulPelet1","ResulPelet2", "ResulPelet3"],["Dosis [Gy]","Std [Gy]", "TasaDosis [Gy/s]"]);
% Nombres de filas
rowNames = {'Pocillo Derecho', 'Pocillo Izquierdo'};
% Asignar nombres a las filas
ResultadosPelets.Properties.RowNames = rowNames;


%% Contour80: 
thresholdBajo = 10; 
threshold_h1 = (thresholdBajo / 100) * max(max(Matrix_h1)); 
Matrix_h1(Matrix_h1<threshold_h1) = 0; 
threshold_h = (thresholdBajo / 100) * max(max(Matrix_h3)); 
Matrix_h3(Matrix_h3<threshold_h) = 0; 

Matrix_h1(Matrix_h1==0) = nan; 
Matrix_h3(Matrix_h3==0) = nan; 

AA_h1 = prctile(Matrix_h1,98,'all'); 
AA_h3 = prctile(Matrix_h3,98,'all');
 
Matrix_h1 = Matrix_h1/AA_h1; 
Matrix_h3 = Matrix_h3/AA_h3; 

% figure(3)
% [C,h] = contour(Xedges, Yedges, Matrix_h1, [0.5 0.8], '-r');
% hold on;
% %contour(xpos, ypos, D2, [0.9 0.9], '-g');
% %contour(xpos, ypos, D3, [0.9 0.9], '-k');
% clabel(C,h)
% [C,h] = contour(Xedges, Yedges, Matrix_h3, [0.5 0.8], '-b');
% clabel(C,h)
% %contour(xpos, ypos, D5, [0.9 0.9], '-m');s
% areaOfInterest = [-1 1 -0.5 0.5];
% xlabel('Ypos [mm]')
% ylabel('Xpos [mm]')
% grid on
% axis(areaOfInterest)
% hold on;
% boundaries1 = bwboundaries(BW);
% boundaries2 = bwboundaries(BW2);
% % Obtén las coordenadas x e y del primer contorno
% Px1 = boundaries1{1}(:, 2);
% Py1 = boundaries1{1}(:, 1);
% % Obtén las coordenadas x e y del segundo contorno
% Px2 = boundaries2{1}(:, 2);
% Py2 = boundaries2{1}(:, 1);
% patch(Px1,Py1,'green','FaceAlpha',0.5);
% hold on; 
% patch(Px2,Py2,'green','FaceAlpha',0.5);
% legend('Nude beam', 'Using compensator')
% title('Isodose lines from Matlab')

% % Crear una figura
% figure(3)
% hold on;
% % Primero, grafica las medias lunas
% boundaries1 = bwboundaries(BW);
% boundaries2 = bwboundaries(BW2);
% % Obtén las coordenadas x e y del primer contorno
% Px1 = boundaries1{1}(:, 2);
% Py1 = boundaries1{1}(:, 1);
% patch(Px1, Py1, 'green', 'FaceAlpha', 0.5);
% % Obtén las coordenadas x e y del segundo contorno
% Px2 = boundaries2{1}(:, 2);
% Py2 = boundaries2{1}(:, 1);
% patch(Px2, Py2, 'green', 'FaceAlpha', 0.5);
% % Ahora, superponer las líneas de contorno
% [C, h] = contour(Xedges, Yedges, Matrix_h1, [0.5 0.8], '-r');
% clabel(C, h)
% [C, h] = contour(Xedges, Yedges, Matrix_h3, [0.5 0.8], '-b');
% clabel(C, h)
% xlabel('Ypos [mm]')
% ylabel('Xpos [mm]')
% grid on
% areaOfInterest = [-1 1 -0.5 0.5];
% axis(areaOfInterest)
% legend('Pocillo', 'Nude beam', 'Using compensator')
% title('Isodose lines from Matlab')

figure(3)
hold on;
% Muestra las medias lunas usando imshow con transparencia
BW_combined = BW + BW2;
imshow(BW_combined, 'InitialMagnification', 'fit', 'XData', [1 width], 'YData', [1 height]);
set(gca, 'YDir', 'normal');  % Asegúrate de que la dirección del eje Y sea normal
h_img = imshow(BW_combined); % Devuelve un objeto de imagen
set(h_img, 'AlphaData', BW_combined*0.5); % 0.5 es la transparencia
% Ahora, superponer las líneas de contorno
hold on;
[C, h] = contour(Xedges, Yedges, Matrix_h1, [0.5 0.8], '-r');
clabel(C, h);
hold on;
[C, h] = contour(Xedges, Yedges, Matrix_h3, [0.5 0.8], '-b');
clabel(C, h);
xlabel('Ypos [mm]');
ylabel('Xpos [mm]');
grid on;
areaOfInterest = [-1 1 -0.5 0.5];
axis(areaOfInterest);
legend('Pocillo', 'Nude beam', 'Using compensator');
title('Isodose lines from Matlab');

%% Perfiles: 
Xedges = Xedges*10; %cm -> mm
Yedges = Yedges*10; %mm

indX0_h = find(Xedges>0,1);
indY0_h = find(Yedges>0,1);

figure(4);
subplot(2,2,1);
plot(Yedges,Matrix_h1(:,indX0_h));
title('Perfil vertical entrada');
legend('Matlab');
xlabel('Ypos (mm)')
ylabel('Dosis (Gy)')
subplot(2,2,2);
plot(Xedges,Matrix_h1(indY0_h,:));
title('Perfil horizontal entrada');
legend('Matlab');
xlabel('Xpos (mm)')
ylabel('Dosis (Gy)')
subplot(2,2,3);
plot(Yedges,Matrix_h3(:,indX0_h));
title('Perfil vertical salida');
legend( 'Matlab');
xlabel('Ypos (mm)')
ylabel('Dosis (Gy)')
subplot(2,2,4);
plot(Xedges,Matrix_h3(indY0_h,:));
title('Perfil horizontal salida');
legend( 'Matlab');
xlabel('Xpos (mm)')
ylabel('Dosis (Gy)')

toc
