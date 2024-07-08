load('Current_results');
load('Design_information');
n_bins=length(optimizing_stage);
fid=fopen(txtname,'w');
    for j=1:1:length(J(:,1))
        fprintf(fid,'%d %d %d %d\n',log10(delay_weight)+J(j,1),J(j,2),J(j,4),J(j,5)*area_scale);
    end
    fclose(fid);
    stage_delay=['stage_delay_',txtname];
    dlmwrite(stage_delay,'J_DB');
    fclose all;    
 Xtxtname=['X_',txtname];
 fid=fopen(Xtxtname,'w');
 for j=1:1:length(J(:,1))
     for i_stage=1:n_bins
         for i=1:2*optimizing_stage(i_stage)+2
             fprintf(fid,'%d ',DB(j).X(i,i_stage));
         end
     end
     fprintf(fid,'\n');
 end
 fclose(fid);

 
