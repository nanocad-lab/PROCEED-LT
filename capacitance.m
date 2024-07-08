function capacitancy = capacitance(wire_len_layer, layer_below, layer, layer_after)
    % capacitancy library for GF12LP
    % cap lib col: 1) layer before, 2) layer level, 3) layer above
    % cap lib col: 4) resistance targ = fF/um
    cap_ref = 'gf12lp_metal_capacitance_ref.csv';
    input_file_folder = 'input_files';
    cap_data = readcell(fullfile(input_file_folder,cap_ref));
    % cap_data
    % strcmp(layer_below,cap_data{1,1})

    layer_below_mark = 0;
    layer_chosen_mark = 0;
    cap_data_chosen = 0;

    for i=1:1:height(cap_data)
        if (strcmp(layer_below,cap_data{i,1}) == 1)
            layer_below_mark = i;
            cap_err = 0;
            break;
        else
            cap_err = 1;
        end
    end

    for i=layer_below_mark:1:height(cap_data)
        if (strcmp(layer,cap_data{i,2}) == 1)
            layer_chosen_mark = i;
            cap_err = 0;
            break;
        else
            cap_err = 1;
        end
    end

    for i=layer_chosen_mark:1:height(cap_data)
        if (strcmp(layer_after,cap_data{i,3}) == 1)
            cap_data_chosen = i;
            cap_err = 0;
            break;
        else
            cap_err = 1;
        end
    end
    

    if (cap_err == 0)
        % the output now is in pF
        capacitancy = cap_data{cap_data_chosen,7}*(10^(-3))*wire_len_layer;
    else
        fprintf ('Capacitance for layer did not exist in %s.\nPlease consult library documentation!\nReverting back to capacitance = 0 to prevent error...\n', cap_ref);
        capacitancy = 0;
    end

    fprintf('len:%.3f - capacitancy: %.3f\n',wire_len_layer,capacitancy);
    % fprintf('tempe: %.3f - tempe_ref: %i\n',tempe,tempe_ref);
end
