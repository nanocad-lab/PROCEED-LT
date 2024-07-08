function [decision_add_point] = is_point_over_Pareto_front(J, log_delay, log_power, mode)
    % This function checks if the new point is over the Pareto front
    switch mode
        case 1
            rankc1 = 1;     % power
            rankc2 = 2;     % delay
        case 2
            rankc1 = 1;
            rankc2 = 5;
        case 3
            rankc1 = 5;
            rankc2 = 2;
    end
    Jtt = J;

    ori_len_Jtt = size(Jtt, 1);
    
    Jtt(ori_len_Jtt+1, 1) = log_delay;
    Jtt(ori_len_Jtt+1, 2) = log_power;
    id_new_point = ori_len_Jtt+1;
    Jtt(ori_len_Jtt+1, 3) = id_new_point;
    Jtt(ori_len_Jtt+1, 4:end) = 0;

    new_len_Jtt = ori_len_Jtt + 1;

    Jtt=sortrows(real(Jtt),[rankc1 3]);       % sort by delay and then by id
    f5 = 0;
    for i=1:1:new_len_Jtt
        if(Jtt(i,rankc1)==38)
            f5=1;
            break;
        end
    end
    if(f5==1)
        Jtt(i:new_len_Jtt,:)=[];
        new_len_Jtt = i - 1;
    end

    i=2;

    while(i<=new_len_Jtt)
        if(Jtt(i,rankc2)>=Jtt(i-1,rankc2))
            Jtt(i,:)=[];
            new_len_Jtt = new_len_Jtt - 1;
            continue;
        end
        i=i+1;
    end

    if(ismember(id_new_point, Jtt(:,3)))
        % the new point is under the Pareto front, has been reserved
        decision_add_point = 1;
    else
        % the new point is over the Pareto front, has been deleted
        decision_add_point = 0;
    end


