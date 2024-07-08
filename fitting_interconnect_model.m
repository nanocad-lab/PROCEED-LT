function [wire_model]=fitting_interconnect_model(file1,file2)
Data=dlmread(file1);
wire_model(1:2)=polyfit(Data(:,1),Data(:,2),1);
% plot(Data(:,1),Data(:,2),'*');
% hold on;
% plot(Data(:,1),Data(:,1)*wire_model(1)+wire_model(2));
Data=[];
Data=dlmread(file2);
wire_model(3:4)=polyfit(Data(:,1),Data(:,2),1);
wire_model(9) = 0; % flag for HGI
% plot(Data(:,1),Data(:,2),'*');
% hold on;
% plot(Data(:,1),Data(:,1)*wire_model(1)+wire_model(2));
