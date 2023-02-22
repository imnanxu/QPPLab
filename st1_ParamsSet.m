clear; clc;  close all
%% For saving & reading path-names
data='HumanVisual'; % the filename of the data file
data='HCPR3gsr';
p2data=['./Input/' data '.mat']; % path of data file which must includes following parameters:
% D0:   a nsbj X nscn cell matrix. Each cell has a nroi X ntimepoints 
%       matrix of EPI timeseries
% MotionInf: a nsbj X nscn cell matrix. Each cell >=1 segments of
%       timepoints without significant motions.
% nY:   # of total networks
% G2Y:  a (nroi X 1) vector including the network# of each ROI
% ibY:  a (nY+1 X 1) vector including the last ROI label of each 
%       functional network;  always set ibY(1)=0.
% iG2Y: a (nY X 1) cell vector including the list of ROI labels for each network.
% YLB: a (nY X 1) cell vector including the shorthand label for each network

ext=['test']; % filename extension for the parameter file
p2param=['Params_' data '_' ext]; % parameter filename
p2qppf='./QPPfv0922/'; %directory of QPP functions
%% Set up global parameters 
%%%%%% for QPP detection %%%%%%
% nP: total # of QPPs to detect (E.g., only detect the primary QPP (QPP1) if nP=1; detect both QPP1 & QPP2 if nP=2, etc.)
% PL: a (nP X 1) vector of QPP widow length, ~20s in humans. E.g.,
% PL(ip)=20s/TR.
nP=5; PL=[30; 30; 30; 29; 29]; 

%% For QPP detection: QPPf1detectRbs.m
% cth13 & cth45: correlation threshold for QPP1-QPP3 & for QPP4-QPP5;
% fd: 1 or 0
%    1 -fast QPP detection; selected limited number of starting points
%    0 -robust detection; select all possible starting points
cth13=[0.1, 0.2]; cth45=[0.1, 0.2]; fd=1; 

%% For QPP phase adjustment: QPPf2phadj.m
% cthph: similarity threshold when phase-adjusting (phadj) a QPP
% s: 1 or 0
%    1 -strict phase adjustment
%    0 -relaxed phase adjustment
cthph=0.88; s=0; 

% Reference parcel(s) ID for QPP phase adjustment: 
% One can >=1 parcel(s). These parcels will have QPP waveform starting from 
% rising postive values (e.g. sine wave).
sdph=cell(nP,1); 
sdph{1}={163;18}; sdph{2}={163; 87}; sdph{3}={18;87}; sdph(4:5)={{18}}; 
%% For network analysis (e.g., QPPf4regscnRbst.m)
% fz: 1 or 0. 
%     1 -output matrix FCr will be the pearson correlation matrix
%     0 -output matrix FCr will be the Fisher Z-Transformation of the pearson correlaion
fz=1; 
save(p2param);