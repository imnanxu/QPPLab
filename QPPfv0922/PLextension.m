function [PLh, PLc, PLe]=PLextension(PL)
% developed by Nan Xu on June 23, 2022.
PLh=round(PL/2)+[0 -rem(PL,2)]; % pad length for temporally extending a QPP
PLc=round(PL/2)+(1:PL); % range of interest in an extended QPP, matches PLh
PLe=PL+sum(PLh); % length of an extended QPP (derivable but saved/read)
