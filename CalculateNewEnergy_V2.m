function [vectorNuevasEnergias] = CalculateNewEnergy_V2 (dz,V_energia, V_perdida, inPHSP)

stepSize = dz.*(1 + (sin(inPHSP(:,3))).^2 + (sin(inPHSP(:,4))).^2).^(1/2); %preguntar
stoppingPower_MeVcm = interp1(V_energia,V_perdida,inPHSP(:,5));

vector = inPHSP(:,5) - stoppingPower_MeVcm.*stepSize;
vectorNuevasEnergias = vector; 