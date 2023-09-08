function [vectorNuevasEnergias] = CalculateNewEnergy_V3 (dz,V_energia, V_perdida, inPHSP,onesIdx)

stepSize = dz.*(1 + (sin(inPHSP(onesIdx,3))).^2 + (sin(inPHSP(onesIdx,4))).^2).^(1/2); %preguntar
stoppingPower_MeVcm = interp1(V_energia,V_perdida,inPHSP(onesIdx,5));

vector = inPHSP(onesIdx,5) - stoppingPower_MeVcm.*stepSize;
vectorNuevasEnergias = vector; 