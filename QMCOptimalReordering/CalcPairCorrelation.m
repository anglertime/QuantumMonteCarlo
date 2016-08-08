% function pairData = CalcPairCorrelation (coords,runningHistogram)
% Calculating the pair correlations.
% global Box
% for ptcl=1:length(coords)
%     for ptcl2=1:length(coords)
%         diff=coords(ptcl2)-coords(ptcl1);
%         diff=put_in_box(diff);
%         dist=sqrt(dot(diff,diff));
%         index=floor(dist/Box(1)*length(runningHistogram)); 
%         if (index>0 && index < length(runningHistogram))
%             runningHistogram(index)=runningHistogram(index)+1; 
%         end 
%      end 
% end
% pairData = runningHistogram; 