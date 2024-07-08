% Modified by CZC

fid=fopen(spname,'a');
fprintf(fid,'\n.data cv\nvd Vtn Vtp ');
for i=1:n
    fprintf(fid,'w%i ',i);
end

fprintf(fid,'\n %d %d %d ',Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)-flag_vtp*dv/2);
for i=4:n+3     % 要打印的width数量！
    wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);       % number of fingers
    if wid <2
        wid=2;
    end
    fprintf(fid,'%d ',wid);
end

% 一阶导
for j=1:n+3
    if(j==1&&flag_vdd)
        fprintf(fid,'\n %d %d %d ',Xc(1,i_stage)+flag_vdd*dv/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)-flag_vtp*dv/2);
        for i=4:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);       % number of fingers
             if wid <2
                wid=2;
             end
            fprintf(fid,'%d ',wid);
        end
    end
    if(j==2&&flag_vtn)
        fprintf(fid,'\n %d %d %d ',Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)+flag_vtn*dv/2, Xc(3,i_stage)-flag_vtp*dv/2);
        for i=4:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
    if(j==3&&flag_vtp)
        fprintf(fid,'\n %d %d %d ',Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)+flag_vtp*dv/2);
        for i=4:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
    if(j>3)
        fprintf(fid,'\n %d %d %d ',Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)-flag_vtp*dv/2);
        for i=4:j-1     
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
        % optimize the j th width
        wid = ceil((1+dwu*flag_w)*Xc(j,i_stage)/width_eff_mult);
        if wid <2
            wid=2;
        end
        fprintf(fid,'%d ',wid);
        for i=j+1:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
end

% 二阶导
for j=1:n+3
    if((j==1&&flag_vdd==0)||(j==2&&flag_vtn==0)||(j==3&&flag_vtp==0))
        continue;
    end
    if(j==1)
        fprintf(fid,'\n %d %d %d ', Xc(1,i_stage)+flag_vdd*dv*3/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)-flag_vtp*dv/2);
        for i=4:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
    if(j==2)
        fprintf(fid,'\n %d %d %d ', Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)+flag_vtn*dv*3/2, Xc(3,i_stage)-flag_vtp*dv/2);
        for i=4:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
    if(j==3)
        fprintf(fid,'\n %d %d %d ', Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)+flag_vtp*dv*3/2);
        for i=4:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
    if(j>3)
        fprintf(fid,'\n %d %d %d ', Xc(1,i_stage)-flag_vdd*dv/2, Xc(2,i_stage)-flag_vtn*dv/2, Xc(3,i_stage)-flag_vtp*dv/2);     % modified by czc: Xc(2,i_stage)-flag_vtp*dv/2 -> Xc(3,i_stage)-flag_vtp*dv/2
        for i=4:j-1
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
        wid = ceil((1+(dwu+dwl+dwu)*flag_w)*Xc(j,i_stage)/width_eff_mult);      % The j th width that needs to be optimized
        if wid <2
            wid=2;
        end
        fprintf(fid,'%d ',wid);
        for i=j+1:n+3
            wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
    end
    for kt=j+1:n+3      % 二阶导，(n+3+1)*(n+3)/2   个导数组合
        if((kt==1&&flag_vdd==0)||(kt==2&&flag_vtn==0)||(kt==3&&flag_vtp==0))        % really critical!!! flag==0就会跳过导数项的计算，所以在derivation_calculation.m中，就要补上缺失导数项的位置。属实傻逼。
            continue;
        end
        fprintf(fid,'\n ');
        for i=1:j-1
            if(i<=3)
                if(i==1)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vdd*dv/2);
                elseif (i==2)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vtn*dv/2);
                elseif (i==3)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vtp*dv/2);
                end
            else
                wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
                if wid <2
                    wid=2;
                end
                fprintf(fid,'%d ',wid);
            end
        end
        if(j<=3)    % j th is the first variable needs to be calculated the derivative
            fprintf(fid,'%d ',Xc(j,i_stage)+dv/2);  % 这里不用乘flag，因为flag为0就会继续下一轮循环（continue）
        else
            wid = ceil((1+dwu*flag_w)*Xc(j,i_stage)/width_eff_mult);    % j th需要被optimize，如果是j th是voltage就优化电压，如果j th是width就优化width
            if wid <2
                wid=2;
            end
        fprintf(fid,'%d ',wid);
        end
        for i=j+1:kt-1
            if(i<=3)
                if(i==1)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vdd*dv/2);
                elseif(i==2)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vtn*dv/2);
                elseif(i==3)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vtp*dv/2);
                end
            else
                wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
                if wid <2
                    wid=2;
                end
                fprintf(fid,'%d ',wid);
            end
        end
        if(kt<=3)       % kt is the second variable needs to be calculated the derivative
            fprintf(fid,'%d ',Xc(kt,i_stage)+dv/2);
        else
            % wid = ceil((1+dwu*flag_w)*Xc(j,i_stage)/width_eff_mult);        % modified by czc: j here should be kt
            wid = ceil((1+dwu*flag_w)*Xc(kt,i_stage)/width_eff_mult);   
            if wid <2
                wid=2;
            end
            fprintf(fid,'%d ',wid);
        end
        for i=kt+1:n+3
            if(i<=3)
                if(i==1)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vdd*dv/2);
                elseif(i==2)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vtn*dv/2);
                elseif(i==3)
                    fprintf(fid,'%d ',Xc(i,i_stage)-flag_vtp*dv/2);
                end
            else
                wid = ceil((1-dwl*flag_w)*Xc(i,i_stage)/width_eff_mult);
                if wid <2
                    wid=2;
                end
                fprintf(fid,'%d ',wid);
            end
        end
    end
end
fprintf(fid,'\n.enddata\n');
fclose(fid);
