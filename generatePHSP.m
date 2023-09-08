clear variables

%% 0. Define variable attributes
VSAD = -190; % mm
VSADY = -190; 
E = 230;
dE = 0.5;
outputFile = 'MapaDosisExpV2';
DesiredNParticles = 7.0390e+09;
thresoldPercent = 10; % Removes points with less than this % of the dose

%% 1. Load image and turn into something understandable
load('processedDoseMaps.mat');

Xpos = Xpositions{1};
Ypos = Ypositions{1};
IMap = doseMaps{1};

% 1. Normalize
IMap = IMap / sum(IMap, 'all');
% 2. Delete points under threshold and renormalize
threshold = (thresoldPercent / 100) * max(max(IMap)); %1
IMap(IMap<threshold) = 0;
IMap = IMap / sum(IMap, 'all');

%% 2. Create particles
validMask = IMap>=threshold;
XposMat = ones(size(Ypos)) * Xpos;
YposMat = Ypos * ones(size(Xpos));
validXpos = XposMat(validMask);
validYpos = YposMat(validMask);
validNParts = round(DesiredNParticles * IMap(validMask));
NParticles = sum(validNParts);

%% 3 Loop to generate final PHSP file

% Common attributes
finalZ = 0;
dCosX = sign(VSAD) * validXpos ./ sqrt(validXpos.^2 + VSAD^2);
dCosY = sign(VSAD) * validYpos ./ sqrt(validYpos.^2 + VSADY^2);
% Para meter dispersión angular: 
% TetaX = acosd(dCosX); 
% TetaY = acosd(dCosY); 
% TetaDX = normrnd(TetaX, ¿???? ); 
% TetaDY = normrnd(TetaY, ¿???? );
% dCosTX = cos(TetaDX); 
% dCosTX = cos(TetaDX);

weight = 1;   
PID = 2212; %PROTONES
flag9 = 1;
flag10 = 1;

% Aux variables
maxE = 0;
minE = Inf;
outputFilePHSP = [outputFile '.phsp'];
outputFileHeader = [outputFile '.header'];

fID = fopen(outputFilePHSP, 'w');
for i=1:numel(validNParts)
    JLim = validNParts(i);
   
    for j=1:JLim
        thisX = 0.1*validXpos(i) + normrnd(0,0.1*pxSize_mm);
        thisY = 0.1*validYpos(i) + normrnd(0,0.1*pxSize_mm);
        thisE = normrnd(E,dE);
        if thisE > maxE
            maxE = thisE;
        end
        if thisE < minE
            minE = thisE;
        end
        lineFormat = '%f\t%f\t%f\t%f\t%f\t%f\t%i\t%i\t%i\t%i\n';
        fprintf(fID,lineFormat,thisX,thisY,finalZ,dCosX(i),dCosY(i),thisE,weight,PID,flag9,flag10);
    end
end
fclose(fID);

%% 4. Write header, necessary to be used in TOPAS
fIDh = fopen(outputFileHeader, 'w');
fprintf(fIDh,'TOPAS ASCII Phase Space\n\n');

fprintf(fIDh,'Number of Original Histories: %i\n',NParticles);
fprintf(fIDh,'Number of Original Histories that Reached Phase Space: %i\n',NParticles);
fprintf(fIDh,'Number of Scored Particles: %i\n\n',NParticles);

fprintf(fIDh,'Columns of data are as follows:\n');
fprintf(fIDh,' 1: Position X [cm]\n');
fprintf(fIDh,' 2: Position Y [cm]\n');
fprintf(fIDh,' 3: Position Z [cm]\n');
fprintf(fIDh,' 4: Direction Cosine X\n');
fprintf(fIDh,' 5: Direction Cosine Y\n');
fprintf(fIDh,' 6: Energy [MeV]\n');
fprintf(fIDh,' 7: Weight\n');
fprintf(fIDh,' 8: Particle Type (in PDG Format)\n');
fprintf(fIDh,' 9: Flag to tell if Third Direction Cosine is Negative (1 means true)\n');
fprintf(fIDh,'10: Flag to tell if this is the First Scored Particle from this History (1 means true)\n\n');

fprintf(fIDh,'Number of proton: %i\n\n', NParticles);

fprintf(fIDh,'Minimum Kinetic Energy of proton: %f MeV\n\n', minE);

fprintf(fIDh,'Maximum Kinetic Energy of proton: %f MeV\n', maxE);
fclose(fIDh);
