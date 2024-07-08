clear;
close all;
mode=1; %1:DP 2:DA 3:AP
switch mode
    case 1 
        c1=1;
        c2=2;
    case 2
        c1=1;
        c2=4;
    case 3
        c1=4;
        c2=2;
end
flag_plot=1;
file_list=['file_list'];
fid=fopen(file_list);
tline=fgetl(fid);
delay_correction=1;
stage=[];
[t r]=strtok(tline);
while(t>0)
    stage(length(stage)+1)=str2num(t);
    [t r]=strtok(r);
end
n_bins=length(stage);
n_file=0;
write_name=fgetl(fid);
all_data=[];
f_DVFS=0;
while(~feof(fid))
    tline=fgetl(fid);
    if(strcmp(tline,'DVFS'))
        f_DVFS=1;
        break;
    end
    n_file=n_file+1;
    file(n_file).result_file=tline;
    file(n_file).data=dlmread(tline);
    file(n_file).data(:,1)=log10(delay_correction)+file(n_file).data(:,1);
    file(n_file).data(:,1+3+2*n_bins+2*sum(stage))=n_file;
    if( n_file>1)
    if(length(all_data(1,:) ) > length(file(n_file).data(1,:)) )
        file(n_file).data(:,length(file(n_file).data(1,:))+1:length(all_data(1,:)) ) = 0;
    else
        all_data(:,length(all_data(1,:))+1 : length(file(n_file).data(1,:))) = 0;
    end
    end
    all_data=[all_data;file(n_file).data];
end
if(f_DVFS)
    n_DVFS=0;
    while(~feof(fid))
        tline=fgetl(fid);
        n_DVFS=n_DVFS+1;
        DVFS(n_DVFS).result_file=tline;
        DVFS(n_DVFS).data=dlmread(tline);
    end
    Th_put=[];
    ith=0;
    for i_DVFS=1:n_DVFS
        for i_c=1:length(DVFS(i_DVFS).data(1,:))/2
            if(isempty(find(Th_put==DVFS(i_DVFS).data(1,i_c*2-1))))
                Th_put=[Th_put,DVFS(i_DVFS).data(1,i_c*2-1)];
                ith=length(Th_put);
                DVFS(ith).final_data=[];
            else
                ith=find(Th_put==DVFS(i_DVFS).data(1,i_c*2-1));
            end
            DVFS(ith).final_data=[DVFS(ith).final_data;DVFS(i_DVFS).data(2:length(DVFS(i_DVFS).data(:,1)),i_c*2-1:i_c*2)];
        end
    end
    for ith=1:length(Th_put)
        temp=sortrows(DVFS(ith).final_data,-1);
        line=1;
        while(line<length(temp))
            if(temp(line+1,2)>=temp(line,2))
                temp(line+1,:)=[];
            else
                line=line+1;
            end
        end
        DVFS(ith).final_data=temp;
    end
end
        
all_data=sortrows(all_data,[c1,c2]);
line=1;
while(line<length(all_data(:,c1)))
    if(all_data(line+1,c2)>=all_data(line,c2))
        all_data(line+1,:)=[];
    else
        line=line+1;
    end
end

type_plot(1).color=['r'];
type_plot(2).color=['g'];
type_plot(3).color=['b'];
type_plot(4).color=['y'];
type_plot(1).marker=['x'];
type_plot(2).marker=['o'];
type_plot(3).marker=['.'];
type_plot(4).marker=['s'];

    h1=figure;
    ha1=axes;
    plot(ha1,all_data(:,c1),all_data(:,c2),'ok','Markersize',14);
    for i_file=1:n_file
        hold on;
        plot(ha1,file(i_file).data(:,c1),file(i_file).data(:,c2),[type_plot(i_file).color,type_plot(i_file).marker]);
    end
   

    
if(flag_plot)
    h2=figure;
    ha2=axes;
    plot(all_data(:,c1),-all_data(:,4+2),'*r','Markersize',12);
    hold on;
    plot(all_data(:,1),-all_data(:,4+2*n_bins-2+2*sum(stage)-2*max(stage)+2),'ob','Markersize',12);
    legend('1st vtshift','2nd vtshift');
%     title('Vt allocation','Fontweight','b'); 
    xlabel('Log(Delay)','fontsize',16);
    ylabel('Vt','fontsize',16);
    set(ha2,'Fontweight','b','fontsize',14);
end
delay_vdd_vt(:,1)=all_data(:,1);
if(flag_plot)
    h3=figure;
    ha3=axes;
    for i=1:n_bins
        hold on;
        if(i==1)
            delay_vdd_vt(:,2) = all_data(:,4+2*i-2+2*sum(stage(1:i-1))+1)+all_data(:,4+2*i-2+2*sum(stage(1:i-1))+2) -0.46;
        end
        if(i==3)
            delay_vdd_vt(:,3) = all_data(:,4+2*i-2+2*sum(stage(1:i-1))+1)+all_data(:,4+2*i-2+2*sum(stage(1:i-1))+2) -0.025;
        end
        plot(all_data(:,1),all_data(:,4+2*i-2+2*sum(stage(1:i-1))+1)+all_data(:,4+2*i-2+2*sum(stage(1:i-1))+2),[type_plot(i).color,'.'],'Markersize',9);
    end
    xlabel('Delay','fontsize',12);
    ylabel('Vdd-Vt,shift','fontsize',16);
    set(ha3,'Fontweight','b','fontsize',16);
end
dlmwrite('DelayVSVdd_Vt.txt',delay_vdd_vt);
if(flag_plot)
    h5=figure;
    ha5=axes;
    whole_sizes=zeros(length(all_data(:,1)),1);
    for i=1:n_bins
        for i_size = 1:i
            whole_sizes = whole_sizes + 2* all_data(:, 4 +2*i + 2* sum(stage(1:i-1)) + i_size*2-1) + all_data(:, 4 +2*i + 2* sum(stage(1:i-1)) + i_size*2);
        end 
        
    end
    plot(all_data(:,1),whole_sizes,[type_plot(i).color,'.'],'Markersize',9);
    xlabel('Delay','fontsize',12);
    ylabel('All sizes','fontsize',16);
    set(ha3,'Fontweight','b','fontsize',16);
    delay_size_data = 10.^all_data(:,1);
    delay_size_data (:,2) = whole_sizes;
    i_p=2;
    while( i_p <= length(delay_size_data(:,1)))
        if( delay_size_data(i_p,2) > delay_size_data(i_p-1,2))
            delay_size_data(i_p,:) = [];
            i_p = i_p-1;
        end
        i_p =i_p+1;
    end
    dlmwrite('Delay_sizes.txt',delay_size_data,'\t');
    
end
all_data(:,1:3)=10.^all_data(:,1:3);
dlmwrite(write_name,all_data,' ');
write_data(:,1) = all_data(:,c1);
c2=4;
write_data(:,2) = all_data(:,c2);
switch mode
    case 1
        DP_file=['Delay_Power_',write_name];
        dlmwrite(DP_file,write_data,'\t');
    case 2
        DP_file=['Delay_Area_',write_name];
        dlmwrite(DP_file,write_data,'\t');
    case 2
        DP_file=['Area_Power_',write_name];
        dlmwrite(DP_file,write_data,'\t');
end
if(f_DVFS)
h4=figure;
ha4=axes;
for ith=1:length(Th_put)
    DVFS_file=[num2str(Th_put(ith)),'_DVFS_',write_name];
    dlmwrite(DVFS_file,DVFS(ith).final_data,'\t');
    hold on;
    plot(log10(DVFS(ith).final_data(:,1)),log10(DVFS(ith).final_data(:,2)),type_plot(mod(ith,4)+1).color,'linewidth',2);
    DVFS_legend(ith,:)=[num2str(Th_put(ith),0.01),'Throughput'];
end
legend(DVFS_legend);
for i_d=1:n_DVFS
    for i_t=1:length(DVFS(i_d).data(1,:))/2
        hold on;
        plot(log10(DVFS(i_d).data(:,i_t*2-1)),log10(DVFS(i_d).data(:,i_t*2)),[type_plot(mod(i_d,4)+1).marker,type_plot(mod(i_t,4)+1).color]);
    end
end

title('DVFS');
fclose all;
end

    
    

