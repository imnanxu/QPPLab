clear; clc; close all
%% Setup data input, output directory & running method
dataext='HumanVisual_test'; % extended filename=[data '_' ext];
% dataext='HCPR3gsr_test'; % extended filename=[data '_' ext];
runM=1; % QPP running method
% runM: 1 -concatenate all D{i,j} as a whole group and detect group QPP
%       2 -concatenate all D{i,:} and detect QPP from all scans of each subject
%       3 -concatenate all D{:,j} and detect QPP from all subjects of each
%       scan
rbstScrn=1; % control for robust QPP detection
%rbstScrn:  1 - scan all possible initial segments for robust QPP detection
%           0 - scan randomly slected initial segments for fast QPP detection
%% Automatically load data & other hidden parameters
fprintf('QPP & FC Visualizations\n'); 
p2param=['Params_' dataext '.mat']; load(p2param); addpath(p2qppf);
% load(p2data, 'D0','MotionInf','nY','G2Y','ibY','iG2Y','YLB'); [nsbj,nscn]=size(D0); 
% 
% ssg=ones(nP,1); ssg(2:end)=PL(2:end); tres=0.7; 
% paramQPPf2=param_QPPf2(nP, PL, cthph, sdph, s);
% paramQPPf4=param_QPPf4(nP, PL, ibY, iG2Y, fz); 
% ITPstp=20; % step to show progress of algorithm when running ITP times

d2O='./Output/';  % directory to outputs files
if runM==1
    Ng=1; a0=[d2O 'GrpQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save GrpQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Grp'; 
elseif runM==2
    Ng=nsbj; a0=[d2O 'SbjQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save SbjQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Sbj';
elseif runM==3   
    Ng=nscn; a0=[d2O 'ScanQPP/']; if ~exist(a0,'dir'), mkdir(a0); end % dir2 save ScanQPPs  
    a=[a0 'interm/']; if ~exist(a,'dir'), mkdir(a); end % save intermediate files
    indn='Scn';
else, error('Unidenfied value for runM.\n')
end
p2S=cell(Ng,nP); p2S0=cell(Ng,1); % pth2 save SbjQPPs
for ig=1:Ng % gp id
    for ip=1:nP % QPP id
        p2S{ig,ip}=[a dataext '_' indn num2str(ig) '_rbst' num2str(rbstScrn) '_qpp' num2str(ip)];
    end
    p2S0{ig}=[a0 dataext '_' indn num2str(ig) '_rbst' num2str(rbstScrn) '_QPPs'];
     
end
%% Visualization
ct=0; ksalpha=0.05; bin=0.08;
for ig=1:Ng
    load(p2S0{ig},'QPPs','TMXs','Cs','METs', 'QPPas', 'TMXas', 'METas', 'Cas', ...
        'Ds', 'Drs', 'Crs', 'FCrs','nY', 'G2Y', 'ibY','iG2Y','YLB', ...
        'ntlist','ITP','PL','ssg','tres',...        
        'cth13','cth45','paramQPPf2','paramQPPf4');      
    for ip=1:nP  
        ct=ct+1;
        f1=figure(1); %qpps
        subplot(Ng,nP,ct); imagesc(QPPs{ip,1}); plotNets(YLB,ibY,PL(ip),0);
        title(['QPP #' num2str(ip)]);

        f2=figure(2); %reverse qpps
        subplot(Ng,nP,ct); imagesc(QPPs{ip,2}); plotNets(YLB,ibY,PL(ip),0);
        title(['rph QPP #' num2str(ip)]);
        if ip==1, ylabel('reverse phase (rph)'); end
        
        
        f3=figure(3); %qpp sliding correlations        
        subplot(Ng,nP,ct); C1=Cs(ip,:); 
        TMX1=TMXs{ip,1}; Met1=METs{ip,1}; TMX2=TMXs{ip,2}; Met2=METs{ip,2};
        plot(C1,'b'); hold on; plot(TMX1,C1(TMX1),'bv'); plot(TMX2,C1(TMX2),'bo'); 
        axis([0 length(C1) -1 1]); 
        set(gca,'XTick',ntlist(1):ntlist(2):sum(ntlist),'YTick',-1:0.2:1); grid on
        title({['QPP #' num2str(ip)],['median max: ' num2str(0.01*round(100*Met1(1))) ...
        ', median \Deltat_{max}:' num2str(0.1*round(10*Met1(2))) ' tps' ...
        ', #max:' num2str(Met1(3))], ['median min: ' num2str(0.01*round(100*Met2(1))) ...
        ', median \Deltat_{min}:' num2str(0.1*round(10*Met2(2))) ' tps' ...
        ', #min:' num2str(Met2(3))]},'FontSize', 8,'fontweight','normal');
        if ip==1, ylabel(['sliding correlation timecourse']); end
           
        f5=figure(5); %qpp sliding correlations before and after qpp regression     
        subplot(Ng,nP,ct); Cr2=Crs{ip}(ip,:);
        plot(C1,'b'); hold on; plot(Cr2,'r'); axis([0 length(C1) -1 1]); 
        set(gca,'XTick',ntlist(1):ntlist(2):sum(ntlist),'YTick',-1:0.2:1); grid on
        title(['QPP #' num2str(ip)],'fontweight','normal');        
        legend('before qpp regression','after ~');
        if ip==1,  ylabel(['sliding correlation timecourse']); end
        
        f6=figure(6); %qpp sliding correlations before and after qpp regression     
        subplot(Ng,nP,ct); 
        [kl,js,ks]=hist2distns(C1, Cr2,'before qpp reg', 'after ~', 'b','r', bin, ksalpha);
        title({['QPP #' num2str(ip)],['KL div=' num2str(kl) ...
        ', JS dist=' num2str(js) ', KS(' num2str(ksalpha) ')=' num2str(ks)]},'fontweight','normal','FontSize', 8);
        legend('before qpp regression','after ~'); 
        if ip==1, ylabel('timecourse histogram'); end
    
        f7=figure(7); %qpp sliding correlations before and after qpp regression     
        subplot(Ng,nP,ct); FC=corr(Ds{ip}'); FC1=triu(FC,0); FC2=tril(FCrs{ip},-1);
        imagesc(FC1+FC2); plotNets(YLB,ibY,PL(ip),1); title(['QPP #' num2str(ip)]);
        if ip==1, ylabel({['FC matrix before (upper)'], 'and after (lower) qpp reg'}); end
    
        
        
    end
    
end