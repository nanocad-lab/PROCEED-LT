Delay1 = max(Delay_i);

Xkc = X;

detail_delay_and_power_simulation_ng_bndry_fast;

if flag_IR_drop == 1
    fprintf('Finish the initial simulation for IR drop power estimation\n');
    current_E=sum(optimizing_stage_weight.*(Energy2));      % Here is still dynamic energy
    current_L=sum(optimizing_stage_weight.*(-Leakage2));
    max_D=max(Delay2);
    current_E=current_E/(delay_weight*max_D)+current_L;     % (Then you don't need to do initial simulation!!! Save a lot of time.)
    Power_total_now = current_E;

    detail_delay_and_power_simulation_ng_bndry_fast_IRdrop;
end

Delay_each_stage = Delay2;

current_E = sum(optimizing_stage_weight .* Energy2);

current_L = sum(optimizing_stage_weight .* (-Leakage2));

max_D = max(Delay_each_stage);

Power_total = current_E / (delay_weight * max_D) + current_L; % final total power