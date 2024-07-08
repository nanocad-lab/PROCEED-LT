LDH=sortrows(LDH,1);
optimizing_stage=sortrows(optimizing_stage,-1);
optimizing_stage_divide(1)=length(find(LDH<=depth_critical*optimizing_stage(1)/max(optimizing_stage)));
for i_stage=2:n_bins-1
    optimizing_stage_divide(i_stage)=length(find(LDH<=depth_critical*optimizing_stage(i_stage)/max(optimizing_stage)))-sum(optimizing_stage_divide(1:i_stage-1));
end
optimizing_stage_divide(n_bins)=length(LDH)-sum(optimizing_stage_divide(1:n_bins-1));   % optimizing_stage_divide is the number of paths in each stage
for i_stage=1:n_bins
    % total depth in each stage / total logic depth in the whole circuit * total gates / i_stage
    optimizing_stage_weight(i_stage)=(total_gates/sum(LDH))*sum(LDH(sum(optimizing_stage_divide(1:i_stage-1))+1:sum(optimizing_stage_divide(1:i_stage))))/optimizing_stage(i_stage);
    optimizing_stage_gate_num(i_stage) = optimizing_stage_weight(i_stage) * optimizing_stage(i_stage);
    fprintf('Optimizing stage weight %i = %2.2f\n', i_stage, optimizing_stage_weight(i_stage));
end
delay_weight=depth_critical/max(optimizing_stage);
% delay_weight = depth_critical/max(optimizing_stage)/1.5;
fprintf('> Depth Crit: %i - Delay Weight Value: %2.2f\n\n', depth_critical,delay_weight);

