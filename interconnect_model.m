function [R,C0]=interconnect_model(devicemodel_allocation,area_model,X,optimizing_stage,optimizing_stage_weight,input_total_area,Ri,Ci)
% total_w1=0;
% total_n1=0;
% for i_stage=1:length(optimizing_stage)
%     total_w1=total_w1+optimizing_stage_weight(i_stage)*sum(X(3:2:optimizing_stage(i_stage)*2+1,i_stage));
%     total_n1=total_n1+optimizing_stage_weight(i_stage)*optimizing_stage(i_stage);
% end
% total_w2=0;
% total_n2=0;
% for i_stage=1:length(optimizing_stage)
%     total_w2=total_w2+optimizing_stage_weight(i_stage)*sum(X(4:2:optimizing_stage(i_stage)*2+2,i_stage));
%     total_n2=total_n2+optimizing_stage_weight(i_stage)*optimizing_stage(i_stage);
% end
% current_ave_w1=total_w1/total_n1;
% current_ave_w2=total_w2/total_n2;
[total_area, G_A ] = chip_area_model(devicemodel_allocation,area_model, X,optimizing_stage,optimizing_stage_weight );
ratio=sqrt(   total_area / input_total_area );
C0=Ci*ratio;
R=Ri*ratio;