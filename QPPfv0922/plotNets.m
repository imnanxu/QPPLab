function plotNets(YLB,ibY,PL,net)
nX2=ibY(end);
PLhs=round(PL/2)+1; PLhe=round(PL/2)+PL;
tick_pos_msk=diff(ibY)/2+ibY(1:end-1);
hold on
for i=2:length(YLB)        
    plot([0 nX2*2],[ibY(i) ibY(i)],'k'); 
end             
set(gca,'XTick',[PLhs (PLhs+PLhe)/2 PLhe],'XTickLabel',[1 round(PL/2) PL],'YTick',tick_pos_msk,'YTickLabel',YLB)

if net==1    
    for i=2:length(YLB)        
        plot([ibY(i) ibY(i)],[0 nX2],'k'); 
    end        
    axis equal; axis([0 nX2 0 nX2]); 
    set(gca,'XTick',tick_pos_msk,'XTickLabel',YLB)
    xtickangle(45)
end
hold off;