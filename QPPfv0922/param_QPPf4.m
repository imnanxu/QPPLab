function paramQPPf4=param_QPPf4(nP, PL, ibY, iG2Y, nz)    
    PLc=cell(nP,1); 
    for ip=1:nP
        [~, PLc{ip}, ~]=PLextension(PL(ip));
    end
    paramQPPf4.PLc=PLc;
    paramQPPf4.ibY=ibY;
    paramQPPf4.iG2Y=iG2Y;
    paramQPPf4.fz=nz;
