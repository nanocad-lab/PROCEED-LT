flag_vt_eq_vd_bndry = 0;
flag_vd_eq_vt_bndry = 0;
decision_change_bndry = 0;
% Select the vd_bndry and vt_bndry of the initial point
if (strcmp(vdSelect, '2vd') == 1)
    differentIndex = find(Xc(1,1:n_bins) ~= Xc(1,1), 1, 'first') - 1;
    if isempty(differentIndex)
        old_vd_bndry = round(mean(optimizing_stage));
        flag_vd_eq_vt_bndry = 1;
    else
        old_vd_bndry = differentIndex;
    end
else
    old_vd_bndry = min(optimizing_stage);
end
if (strcmp(vtSelect, '2vt') == 1)
    differentIndex = find(Xc(2,1:n_bins) ~= Xc(2,1), 1, 'first') - 1;
    if isempty(differentIndex)
        % old_vt_bndry = round(mean(optimizing_stage));
        old_vt_bndry = round(mean(optimizing_stage));
        flag_vt_eq_vd_bndry = 1;
    else
        old_vt_bndry = differentIndex;
    end
else
    old_vt_bndry = min(optimizing_stage);
end
vdd1 = Xc(1,old_vd_bndry);
vdd2 = Xc(1,old_vd_bndry+1);
vtn1 = Xc(2,old_vt_bndry);
vtn2 = Xc(2,old_vt_bndry+1);
vtp1 = Xc(3,old_vt_bndry);
vtp2 = Xc(3,old_vt_bndry+1);


Xc_old = Xc;
curr_delay = max(J_DB(l,1:n_bins));
% Get deta of all kinds of combinations of bounadries using quick algorithms.
vd_bndry_left = 0;
vd_bndry_right = max(optimizing_stage);
vt_bndry_left = 0;
vt_bndry_right = max(optimizing_stage);
bndry_pd_set = [];
% for bndry_pd_set array arrangement
delay_ind = 1;
power_total_ind = 2;
power_leak_ind = 3;
vd_bndry_ind = 4;
vt_bndry_ind = 5;
delay_each_stage_ind = 6:6+n_bins-1;
power_dyn_each_stage_ind = 6+n_bins:6+2*n_bins-1;
power_leak_each_stage_ind = 6+2*n_bins:6+3*n_bins-1;
bndry_diff_ind = 6+3*n_bins;

fprintf("Reselecting the bins allocation vd_bndry!!!\n");
fprintf("The old vd_bndry is %i\n", old_vd_bndry);
fprintf("The old vt_bndry is %i\n", old_vt_bndry);

for vd_bndry_cal = [vd_bndry_left, vd_bndry_right]        % actually you only need to calculate two boundries (1 and 10-1)
    if (strcmp(vdSelect, '1vd') == 1 && strcmp(vtSelect, '1vt') == 1)
        break;
    end
    if (strcmp(vdSelect, '2vd') == 1)
        Xc(1,1:vd_bndry_cal) = vdd1;
        Xc(1,vd_bndry_cal+1:n_bins) = vdd2;
    end
    for vt_bndry_cal = [vt_bndry_left, vt_bndry_right]
        % bndry_pd_set_row = (vd_bndry_cal-1)*(max(optimizing_stage)-1) + vt_bndry_cal;
        bndry_pd_set_row = (vd_bndry_cal) * (max(optimizing_stage)+1) + vt_bndry_cal + 1;
        if (strcmp(vtSelect, '2vt') == 1)
            Xc(2,1:vt_bndry_cal) = vtn1;
            Xc(2,vt_bndry_cal+1:n_bins) = vtn2;
            Xc(3,1:vt_bndry_cal) = vtp1;
            Xc(3,vt_bndry_cal+1:n_bins) = vtp2;
        end
        % Calculate the delay and power of the new vd_bndry
        Delay_each_stage = [];
        Power_dyn_each_stage = [];
        Power_leak_each_stage = [];
        current_E = 0;
        current_L = 0;
        Power_total = 0;
        
        Xkc = Xc;
        detail_delay_and_power_simulation_ng_bndry_fast;
        Delay_each_stage = Delay2;
        Power_dyn_each_stage = optimizing_stage_weight .* Energy2;
        Power_leak_each_stage = optimizing_stage_weight .* (-Leakage2);
        current_E = sum(Power_dyn_each_stage);
        current_L = sum(Power_leak_each_stage);

        % for i_stage=1:n_bins
        %     fprintf ('vd boundry: %i, vt boundry: %i, i_stage:%i - Power_avg:%d, Power_leak:%d\n', vd_bndry_cal, vt_bndry_cal, i_stage, Power_dyn_each_stage(i_stage), Power_leak_each_stage(i_stage));
        % end
        max_D = max(Delay_each_stage);
        Power_total=current_E/(delay_weight*max_D)+current_L; % final total power
        current_E = Power_total;
        % bndry_diff_ind: column {4+3*n_bins+1}. It is the difference between the current vd_bndry and the previous vd_bndry
        if strcmp(vdSelect, '1vd') == 1
            bndry_diff = abs(vt_bndry_cal - old_vt_bndry);
        elseif strcmp(vtSelect, '1vt') == 1
            bndry_diff = abs(vd_bndry_cal - old_vd_bndry);
        else
            bndry_diff = max(abs(vd_bndry_cal - old_vd_bndry), abs(vt_bndry_cal - old_vt_bndry));
        end
        % Record the delay and power of the new vd_bndry
        % delay, total power, leakage power, vd_bndry, delay_each_stage, Dynamic power_each_stage, leakage power_each_stage
        bndry_pd_set(bndry_pd_set_row,:) = [log10(max_D), log10(Power_total), log10(current_L), vd_bndry_cal, vt_bndry_cal, Delay_each_stage, Power_dyn_each_stage, Power_leak_each_stage, bndry_diff];      
        % DB(l).X=Xc;
    end
end

% To get all the data
% index arrangement: [delay_ind, power_total_ind, power_leak_ind, vd_bndry_ind, vt_bndry_ind, delay_each_stage_ind, power_dyn_each_stage_ind, power_leak_each_stage_ind]
% Four different combinations of vd_bndry and vt_bndry:
% vd_bndry_left, vt_bndry_left
% vd_bndry_left, vt_bndry_right
% vd_bndry_right, vt_bndry_left
% vd_bndry_right, vt_bndry_right
bndry_pd_set_row_comb1 = (vd_bndry_left) * (max(optimizing_stage)+1) + vt_bndry_left + 1;
bndry_pd_set_row_comb2 = (vd_bndry_right) * (max(optimizing_stage)+1) + vt_bndry_left + 1;
bndry_pd_set_row_comb3 = (vd_bndry_left) * (max(optimizing_stage)+1) + vt_bndry_right + 1;
bndry_pd_set_row_comb4 = (vd_bndry_right) * (max(optimizing_stage)+1) + vt_bndry_right + 1;
for vd_bndry_cal = vd_bndry_left : vd_bndry_right
    for vt_bndry_cal = vt_bndry_left : vt_bndry_right
        % bndry_pd_set_row = (vd_bndry_cal-1)*(max(optimizing_stage)-1) + vt_bndry_cal;
        bndry_pd_set_row = (vd_bndry_cal) * (max(optimizing_stage)+1) + vt_bndry_cal + 1;
        if ((vd_bndry_cal == vd_bndry_left && vt_bndry_cal == vt_bndry_left) || ...
            (vd_bndry_cal == vd_bndry_left && vt_bndry_cal == vt_bndry_right) || ...
            (vd_bndry_cal == vd_bndry_right && vt_bndry_cal == vt_bndry_left) || ...
            (vd_bndry_cal == vd_bndry_right && vt_bndry_cal == vt_bndry_right))
            continue;
        end
        if (vd_bndry_cal <= vt_bndry_cal)   % pattern 4, pattern 3, pattern 1
            Delay_each_stage = [bndry_pd_set(bndry_pd_set_row_comb4,delay_each_stage_ind(1:vd_bndry_cal)), ...      % Numbers of bins in pattern 4
                                bndry_pd_set(bndry_pd_set_row_comb3,delay_each_stage_ind(vd_bndry_cal+1:vd_bndry_cal+(vt_bndry_cal-vd_bndry_cal))), ...     % Numbers of bins in pattern 3
                                bndry_pd_set(bndry_pd_set_row_comb1,delay_each_stage_ind(vd_bndry_cal+(vt_bndry_cal-vd_bndry_cal)+1:end))];     % Numbers of bins in pattern 1
            Power_dyn_each_stage = [bndry_pd_set(bndry_pd_set_row_comb4,power_dyn_each_stage_ind(1:vd_bndry_cal)), ...      % Numbers of bins in pattern 4
                                    bndry_pd_set(bndry_pd_set_row_comb3,power_dyn_each_stage_ind(vd_bndry_cal+1:vd_bndry_cal+(vt_bndry_cal-vd_bndry_cal))), ...     % Numbers of bins in pattern 3
                                    bndry_pd_set(bndry_pd_set_row_comb1,power_dyn_each_stage_ind(vd_bndry_cal+(vt_bndry_cal-vd_bndry_cal)+1:end))];     % Numbers of bins in pattern 1
            Power_leak_each_stage = [bndry_pd_set(bndry_pd_set_row_comb4,power_leak_each_stage_ind(1:vd_bndry_cal)), ...      % Numbers of bins in pattern 4
                                    bndry_pd_set(bndry_pd_set_row_comb3,power_leak_each_stage_ind(vd_bndry_cal+1:vd_bndry_cal+(vt_bndry_cal-vd_bndry_cal))), ...     % Numbers of bins in pattern 3
                                    bndry_pd_set(bndry_pd_set_row_comb1,power_leak_each_stage_ind(vd_bndry_cal+(vt_bndry_cal-vd_bndry_cal)+1:end))];     % Numbers of bins in pattern 1
        else    % pattern 4, pattern 2, pattern 1
            Delay_each_stage = [bndry_pd_set(bndry_pd_set_row_comb4,delay_each_stage_ind(1:vt_bndry_cal)), ...      % Numbers of bins in pattern 4
                                bndry_pd_set(bndry_pd_set_row_comb2,delay_each_stage_ind(vt_bndry_cal+1:vt_bndry_cal+(vd_bndry_cal-vt_bndry_cal))), ...     % Numbers of bins in pattern 2
                                bndry_pd_set(bndry_pd_set_row_comb1,delay_each_stage_ind(vt_bndry_cal+(vd_bndry_cal-vt_bndry_cal)+1:end))];     % Numbers of bins in pattern 1
            Power_dyn_each_stage = [bndry_pd_set(bndry_pd_set_row_comb4,power_dyn_each_stage_ind(1:vt_bndry_cal)), ...      % Numbers of bins in pattern 4
                                    bndry_pd_set(bndry_pd_set_row_comb2,power_dyn_each_stage_ind(vt_bndry_cal+1:vt_bndry_cal+(vd_bndry_cal-vt_bndry_cal))), ...     % Numbers of bins in pattern 2
                                    bndry_pd_set(bndry_pd_set_row_comb1,power_dyn_each_stage_ind(vt_bndry_cal+(vd_bndry_cal-vt_bndry_cal)+1:end))];     % Numbers of bins in pattern 1
            Power_leak_each_stage = [bndry_pd_set(bndry_pd_set_row_comb4,power_leak_each_stage_ind(1:vt_bndry_cal)), ...      % Numbers of bins in pattern 4
                                    bndry_pd_set(bndry_pd_set_row_comb2,power_leak_each_stage_ind(vt_bndry_cal+1:vt_bndry_cal+(vd_bndry_cal-vt_bndry_cal))), ...     % Numbers of bins in pattern 2
                                    bndry_pd_set(bndry_pd_set_row_comb1,power_leak_each_stage_ind(vt_bndry_cal+(vd_bndry_cal-vt_bndry_cal)+1:end))];     % Numbers of bins in pattern 1
        end
        max_D = max(Delay_each_stage);
        current_L = sum(Power_leak_each_stage);
        current_E = sum(Power_dyn_each_stage);
        Power_total = current_E / (delay_weight * max_D) + current_L;
        current_E = Power_total;
        if strcmp(vdSelect, '1vd') == 1
            bndry_diff = abs(vt_bndry_cal - old_vt_bndry);
        elseif strcmp(vtSelect, '1vt') == 1
            bndry_diff = abs(vd_bndry_cal - old_vd_bndry);
        else
            bndry_diff = max(abs(vd_bndry_cal - old_vd_bndry), abs(vt_bndry_cal - old_vt_bndry));
        end
        bndry_pd_set(bndry_pd_set_row,:) = [log10(max_D), log10(Power_total), log10(current_L), vd_bndry_cal, vt_bndry_cal, Delay_each_stage, Power_dyn_each_stage, Power_leak_each_stage, bndry_diff];
    end
end



old_bndry_pd_set_row = (old_vd_bndry) * (max(optimizing_stage)+1) + old_vt_bndry + 1;
old_bndry_pd_set = bndry_pd_set(old_bndry_pd_set_row,:);
old_delay = 10^bndry_pd_set(old_bndry_pd_set_row, delay_ind);
old_power_total = 10^bndry_pd_set(old_bndry_pd_set_row, power_total_ind);



% Find the best vd_bndry
if (strcmp(vdSelect, '1vd') == 1 && strcmp(vtSelect, '1vt') == 1)
    new_vd_bndry = old_vd_bndry;
end
if (strcmp(vtSelect, '1vt') == 1)
    new_vt_bndry = old_vt_bndry;
end
if (strcmp(vdSelect, '2vd') == 1 || strcmp(vtSelect, '2vt') == 1)
    % We need to redefine the vtn/vtp/vdd allocation
    % Sort the bndry_pd_set, found the best vd_bndry
    bndry_pd_set=sortrows(bndry_pd_set,[delay_ind power_total_ind]); % 1: delay, 2: total power, 4: vd_bndry_ind, 5: vt_bndry_ind
    del_row_ind = [];
    for i = 2:size(bndry_pd_set,1)    % for each row
        % if one point's delay and power all greater than the previous point, obviously it is not a Pareto point, delete it
        if ((bndry_pd_set(i,delay_ind) > bndry_pd_set(1,delay_ind) && bndry_pd_set(i,power_total_ind) >= bndry_pd_set(1,power_total_ind)) || ...
                (bndry_pd_set(i,delay_ind) >= bndry_pd_set(1,delay_ind) && bndry_pd_set(i,power_total_ind) > bndry_pd_set(1,power_total_ind)))
            del_row_ind = [del_row_ind, i];
        end
    end
    % fprintf("Going to delete the vd_bndry selection: %i\n", del_row_ind);
    bndry_pd_set(del_row_ind,:) = [];
    bndry_pd_set=sortrows(bndry_pd_set, [delay_ind power_total_ind bndry_diff_ind]);
    % fprintf('Selecting from the boundries: %i, %i\n', bndry_pd_set(:,vd_bndry_ind:vt_bndry_ind));
    new_vd_bndry = bndry_pd_set(1,vd_bndry_ind);
    new_vt_bndry = bndry_pd_set(1,vt_bndry_ind);

    if (flag_vt_eq_vd_bndry == 1)
        new_vt_bndry = new_vd_bndry;
    end
    if (flag_vd_eq_vt_bndry == 1)
        new_vd_bndry = new_vt_bndry;
    end
    
    fprintf('Selection ends! The old [vd_bndry,vt_bndry] is [%d,%d], The new [vd_bndry,vt_bndry] is [%d,%d]\n', old_vd_bndry, old_vt_bndry, new_vd_bndry, new_vt_bndry);
    new_delay = 10^bndry_pd_set(1,delay_ind);
    new_power_total = 10^bndry_pd_set(1,power_total_ind);
    
    delay_improvement = (old_delay-new_delay)/old_delay*100;
    power_improvement = (old_power_total-new_power_total)/old_power_total*100;
    fprintf('The old delay is %d, The new delay is %d, improvement: %.2f%%\n', old_delay, new_delay, delay_improvement);
    fprintf('The old power_total is %d, The new power_total is %d, improvement: %.2f%%\n', old_power_total, new_power_total, power_improvement);
    
    if (delay_improvement < 0 || power_improvement < 0)
        decision_change_bndry = 0;
        fprintf('The new boundry is not strictly better than the old boundry, old boundry reserved.\n');
        new_vd_bndry = old_vd_bndry;
        new_vt_bndry = old_vt_bndry;
    elseif ((new_vd_bndry ~= old_vd_bndry) || (new_vt_bndry ~= old_vt_bndry))
        % Check if the point with new vd/vt allocation boundary is over the current Pareto front
        decision_change_bndry = is_point_over_Pareto_front(J, bndry_pd_set(1,delay_ind), bndry_pd_set(1,power_total_ind), mode);
        if (decision_change_bndry == 0)
            fprintf('The new boundry is over current Pareto Front, old boundry reserved.\n');
            new_vd_bndry = old_vd_bndry;
            new_vt_bndry = old_vt_bndry;
        else
            fprintf('The new boundry is under current Pareto Front, new boundry reserved.\n');
        end
    end

    % Update the vdd_allocation
    allocation_len = max(new_vd_bndry, max(optimizing_stage)-new_vd_bndry);
    if(strcmp(vdSelect, '2vd') == 1)
        if (new_vd_bndry == 0)
            new_vd_bndry = new_vd_bndry + 1;
        end
        if (new_vd_bndry == max(optimizing_stage))
            new_vd_bndry = new_vd_bndry - 1;
        end
        vdd_allocation = [[1:new_vd_bndry, zeros([1, allocation_len-new_vd_bndry])]; ...
                            [new_vd_bndry+1:max(optimizing_stage), zeros([1, allocation_len-(max(optimizing_stage)-new_vd_bndry)])]];
    end
    % Update the Xc
    Xc(1,1:new_vd_bndry) = vdd1;
    Xc(1,new_vd_bndry+1:n_bins) = vdd2;

    % Update the vtn_allocation
    allocation_len = max(new_vt_bndry, max(optimizing_stage)-new_vt_bndry);
    if(strcmp(vtSelect, '2vt') == 1)
        if (new_vt_bndry == 0)
            new_vt_bndry = new_vt_bndry + 1;
        end
        if (new_vt_bndry == max(optimizing_stage))
            new_vt_bndry = new_vt_bndry - 1;
        end
        vtn_allocation = [[1:new_vt_bndry, zeros([1, allocation_len-new_vt_bndry])]; ...
                            [new_vt_bndry+1:max(optimizing_stage), zeros([1, allocation_len-(max(optimizing_stage)-new_vt_bndry)])]];
        vtp_allocation = vtn_allocation;
    end
    % Update the Xc
    Xc(2,1:new_vt_bndry) = vtn1;
    Xc(2,new_vt_bndry+1:n_bins) = vtn2;
    Xc(3,1:new_vt_bndry) = vtp1;
    Xc(3,new_vt_bndry+1:n_bins) = vtp2;
    

    %%%Check max metric violation
    [Current_area, temp_G_A] = chip_area_model(devicemodel_allocation, area_model, Xc, optimizing_stage, optimizing_stage_weight);
    f5=0;
    switch mode
        case 1
            if( Current_area /area_scale  > max_limit_metric )
                f5=1;
            end
        case 2
            if( current_E > max_limit_metric )
                f5=1;
            end
        case 3
            if( max(bndry_pd_set(1,4:4+n_bins-1)) > max_limit_metric )
                f5=1;
            end
    end
    if (f5 ==1)
        cout_max_metric_vio = cout_max_metric_vio + 1;
        switch mode
            case 1
                message=['violate max limit metric ', num2str(Current_area ),' > ', num2str(max_limit_metric * area_scale)];
            case 2
                message=['violate max limit metric ', num2str(current_E ),' > ', num2str(max_limit_metric )];
            case 3
                message=['violate max limit metric ', num2str( max(J_DB(nof(k-1)+ij,:)) * delay_weight ),' > ', num2str(max_limit_metric *delay_weight)];
        end
        Xc = Xc_old;
    else
        if (decision_change_bndry == 0)
            % The new boundry point is over current Pareto Front, old boundry reserved.
            % Reserve the old DB, J_DB, J
            X = Xc_old;
            DB(l).X=X;                                     % Pareto Point
            J_DB(l,:)=old_bndry_pd_set(1,delay_each_stage_ind);     % Stage delay
            % J(l,1)=old_bndry_pd_set(1,delay_ind);                 % delay
            % J(l,2)=old_bndry_pd_set(1,power_total_ind);                 % total power
            % J(l,4)=old_bndry_pd_set(1,power_leak_ind);        % leakage power
            % J(l,5)=Current_area /area_scale;            % area
        else
            % The new boundry point is under current Pareto Front, new boundry reserved.
            % Update the DB, J_DB, J
            X = Xc;
            DB(l).X=X;                                     % Pareto Point
            J_DB(l,:)=bndry_pd_set(1,delay_each_stage_ind);     % Stage delay     
            J(l,1)=bndry_pd_set(1,delay_ind);                 % delay
            J(l,2)=bndry_pd_set(1,power_total_ind);                 % total power
            J(l,4)=bndry_pd_set(1,power_leak_ind);        % leakage power
            J(l,5)=Current_area /area_scale;            % area
        end
    end
end