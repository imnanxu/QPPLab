%% Setup parameters
%%% load unreorganized data
dataIN='./Input/B_GSR_HCPR3.mat'; D0=load(dataIN,'B','MotionInf'); D0=D0.B; 
%%% define output filename
dataOUT='./Input/HCPR3gsr.mat'; 
%%% load parcel label and network spreedsheet
AtlasTable = readtable(['./resources/Glasser2016_360Parcels_7Networks.xlsx']);
Label='Parcel_ID'; % this should be variable name of the numerical label of ROIs
System='Network'; % this is the variable name of the numerical label of functional networks
%%% add function files
addpath('./QPPfv0922/')
%% code
AtlasTable.Label=eval(['AtlasTable.' Label]); AtlasTable.System=eval(['AtlasTable.' System]);
AtlasTable = sortrows(AtlasTable,'Label','ascend');
AtlasTable_newLabel = sortrows(AtlasTable,'System','ascend');
[~, ia, ~]=unique(AtlasTable_newLabel.System); 
Nroi=size(AtlasTable,1); ibY=[0, ia(2:end)', Nroi]; nY=length(ibY)-1;
for nid=1:nY, iG2Y{nid}=[ibY(nid)+1:ibY(nid+1)]; G2Y(ibY(nid)+1:ibY(nid+1))=nid; end
[YLB, net_index, G2Y0]=unique(AtlasTable.System);  % AtlasTable.System(ia(i))=net_id{i}; net_id{net_id_index(i)}=AtlasTable.System(i);
[nsbj, nscn]=size(D0);

figure; ct=0;
for i=1:nsbj
    for j=1:nscn
        ct=ct+1;
        data0=D0{i,j};
        data=[G2Y0, data0];
        data_sort=sortrows(data, 1);
        data_sort=data_sort(:,2:end);

        D0{i,j}=data_sort;   
        subplot(nsbj,nscn,ct)
        title(['sbj' num2str(i) ' scn' num2str(j)])
        imagesc(corr(data_sort')); 
        plotNets(YLB,ibY,30,1)
        if ~isfield(D0, 'MotionInf')
            MotionInf{i,j}=[1:size(data_sort,2)];
        end
    end
end
save(dataOUT,'D0','MotionInf','ibY','iG2Y','G2Y','nY','YLB');