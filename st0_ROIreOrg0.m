%% Setup parameters
%%% load unreorganized data
dataIN='./Input/long_task.mat'; D0=load(dataIN,'Dcsfr','MotionInf'); D0=D0.Dcsfr; 
%%% define output filename
dataOUT='./Input/long_task_csfr.mat'; 
%%% load parcel label and network spreedsheet
AtlasTable = readtable(['./resources/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv']);
Label='ROI_Label'; % this should be variable name of the numerical label of ROIs
System='Network'; % this is the variable name of the numerical label of functional networks
%%% add function files
addpath('./QPPfv0922/')
%% code
AtlasTable.Label=eval(['AtlasTable.' Label]); AtlasTable.System=eval(['AtlasTable.' System]);
[NetLB, ~, ROI2Net]=unique(AtlasTable.System, 'stable'); 
Nroi=size(AtlasTable,1); [nsbj, nscn]=size(D0);

if ~isfield(D0, 'MotionInf'),MotionInf=cell(size(D0)); end
ct=0;
for i=1:nsbj
    for j=1:nscn
        ct=ct+1;
        if ~isfield(D0, 'MotionInf')
            MotionInf{i,j}=[1:size(D0{i,j},2)];
        end
    end
end
save(dataOUT,'D0','MotionInf','ROI2Net','NetLB');