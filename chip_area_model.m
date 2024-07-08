function [ total_area, G_A ] = chip_area_model(devicemodel_allocation, area_model, X,optimizing_stage,optimizing_stage_weight )

    total_area = 0; %%initial total area to 0
    G_A = [];
    for i_stage = 1: length(optimizing_stage)
        if( isempty( find(devicemodel_allocation(1,:) == i_stage)) == 0)
            G_A(1:3,i_stage) = 0;
            % fprintf('1.opt_stg: %i\n',optimizing_stage(i_stage)*2);
            for i_gate = 4 : 2 : optimizing_stage(i_stage)*2 + 2
                % fprintf('1.i_gate: %i\n',i_gate);
                total_area = total_area + ( X(i_gate, i_stage) * area_model(1) + area_model(2) ) *optimizing_stage_weight(i_stage) ;
                G_A(i_gate,i_stage) = area_model(1) * optimizing_stage_weight(i_stage);
            end
            for i_gate = 5 : 2 : optimizing_stage(i_stage)*2 + 3
                % fprintf('2.i_gate: %i\n',i_gate);
                total_area = total_area + ( X(i_gate, i_stage) * area_model(3) + area_model(4) ) *optimizing_stage_weight(i_stage) ;
                G_A(i_gate,i_stage) = area_model(3) * optimizing_stage_weight(i_stage);
            end
        else
            G_A(1:3,i_stage) = 0;
            % fprintf('2.opt_stg: %i\n',optimizing_stage(i_stage)*2);
            for i_gate = 4 : 2 : optimizing_stage(i_stage)*2 + 2;
                % fprintf('1.i_gate: %i\n',i_gate);
                total_area = total_area + ( X(i_gate, i_stage) * area_model(5) + area_model(6) ) *optimizing_stage_weight(i_stage) ;
                G_A(i_gate,i_stage) = area_model(5) * optimizing_stage_weight(i_stage);
            end
            for i_gate = 5 : 2 : optimizing_stage(i_stage)*2 + 3;
                % fprintf('2.i_gate: %i\n',i_gate);
                total_area = total_area + ( X(i_gate, i_stage) * area_model(7) + area_model(8) ) *optimizing_stage_weight(i_stage) ;
                G_A(i_gate,i_stage) = area_model(7) * optimizing_stage_weight(i_stage);
            end
        end
    end
    % pause(100000);