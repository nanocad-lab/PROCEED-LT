function resistivity = resistive(tempe,mode,wire_len_per_net,wireload)
    tempe_ref = 25;
    if mode == 0
        total_resistivity = 0;
        net_count = 9999999;
    elseif (mode == 1 || mode == 3)
        % r total from spef file
        total_resistivity = 7.72090903E+06;
        % number of distributed nets from spef file
        d_net_count = 32151;
        % just conversion for final resistance output
        net_count = d_net_count;
    elseif mode == 2
        % mode 2 is accurate per temperature. calibrated at 25c.
        % load design data inside input_files folder.
        wire_data = 'm3_wires_vias.txt';
        metal_ref = 'gf12lp_metal_resistance_ref.csv';
        via_ref = 'gf12lp_via_resistance_ref.csv';
        input_file_folder = 'input_files';

        fid = fopen(fullfile(input_file_folder,wire_data));
        metal_resistance_data = readcell(fullfile(input_file_folder,metal_ref));
        via_resistance_data = readcell(fullfile(input_file_folder,via_ref));

        layer_metal = [];
        via_metal = [];

        tline = fgetl(fid);
        tline_found = 0;
        % 1: extract total number of nets
        while ischar(tline)
            [store rem] = strtok(tline);
            if (strcmp(store,'#Number') == 1)
                [store rem] = strtok(rem);
                if (strcmp(store,'of') == 1)
                    [store rem] = strtok(rem);
                    if (strcmp(store,'nets') == 1)
                        tline_found = 1;
                        break;
                    end
                end
            end
            tline = fgetl(fid);
        end
        % fprintf('whats here: %i\n',net_count);
        % pause(235423);

        if (tline_found == 1)
            for i = 1:1:2
                [store rem] = strtok(rem);
            end
        end
        net_count = str2num(store);     % in 'm3_wires_vias.txt' file, the number of nets is 40239
        % fprintf('whats here: %i\n',net_count);
        % pause(235423);

        frewind(fid)
        tline = fgetl(fid);
        tline_found = 0;
        % 2: extract per-layer metal length
        while ischar(tline)
            [store rem] = strtok(tline);
            if (strcmp(store,'#Total') == 1)
                [store rem] = strtok(rem);
                if (strcmp(store,'wire') == 1)
                    [store rem] = strtok(rem);
                    if (strcmp(store,'length') == 1)
                        tline_found = 1;
                        break;
                    end
                end
            end
            tline = fgetl(fid);
        end
        % fprintf('whats here: %s\n',tline);
        % pause(235423);

        if (tline_found == 1)
            for i = 1:1:2
                tline = fgetl(fid);
            end

            layer_num = 0;
            tline_metal_found = 0;
            while ischar(tline)
                layer_num = layer_num + 1;
                [store rem] = strtok(tline);
                if (strcmp(store,'#Total') == 1)
                    [store rem] = strtok(rem);
                    if (strcmp(store,'wire') == 1)
                        [store rem] = strtok(rem);
                        if (strcmp(store,'length') == 1)
                            tline_metal_found = 1;
                        end
                    end
                end
                
                if (tline_metal_found == 1)
                    for i = 1:1:3
                        [store rem] = strtok(rem);
                    end
                    % fprintf('whats here: %s\n',store);
                    % pause(235423);

                    % store layer name
                    layer_metal{layer_num,1} = store;

                    for i = 1:1:2
                        [store rem] = strtok(rem);
                    end
                    % fprintf('whats here: %s\n',store);
                    % fprintf('whats here: %i\n',str2num(store));
                    % pause(235423);

                    % store wire length of that particular layer in um
                    layer_metal{layer_num,2} = str2num(store);
                    % fprintf('whats here: %i\n',layer_metal{layer_num,2});
                    % pause(235423);

                    % calc resistance
                    metal_notfound = 1;
                    for i=1:1:height(metal_resistance_data)
                        if (strcmp(layer_metal{layer_num,1},metal_resistance_data{i,1}) == 1)
                            metal_notfound = 0;
                            resistance_targ = metal_resistance_data{i,5};
                            tcr_orig = metal_resistance_data{i,7};
                            break;
                        else
                            metal_notfound = 1;
                        end
                    end
                    % fprintf('whats here: %d - %d\n',resistance_targ,tcr_orig);
                    % pause(235423);
                    % fprintf('metal notfound: %d\n',metal_notfound);

                    if (metal_notfound == 0)
                        layer_metal{layer_num,3} = resistivity_calc(resistance_targ,layer_metal{layer_num,2},tcr_orig,tempe,tempe_ref);
                        % fprintf('res_targ: %d - tcl: %d\n',metal_resistance_data{i,5},metal_resistance_data{i,7});
                        % function resistivity = resistivity_calc(resistance_targ,value,tcr_orig,tempe,tempe_ref)
                    else
                        fprintf('Metal resistance data is not right!\nConsult the documentation of this library\n')
                    end
                    % fprintf('whats here: %d\n',layer_metal{layer_num,3});
                    % pause(235423);

                    tline_metal_found = 0;
                end
                tline = fgetl(fid);
            end
        else
            fprintf('Metal layer report is not right!\nPlease regenerate data using this command:\n>> pdi report_design -summary\n');
        end

        % layer_metal
        % pause(2352356);

        % 3: extract vias
        frewind(fid)
        tline = fgetl(fid);
        tline_found = 0;
        while ischar(tline)
            [store rem] = strtok(tline);
            if (strcmp(store,'#Vias') == 1)
                [store rem] = strtok(rem);
                if (strcmp(store,'used') == 1)
                    [store rem] = strtok(rem);
                    if (strcmp(store,'for') == 1)
                        [store rem] = strtok(rem);
                        if (strcmp(store,'rule') == 1)
                            tline_found = 1;
                            break;
                        end
                    end
                end
            end
            tline = fgetl(fid);
        end
        % fprintf('whats tline: %i\n',tline_found);
        % fprintf('whats here BoT: %s\n',tline);
        % pause(235423);

        if (tline_found == 1)
            via_num = 0;
            tline = fgetl(fid);
            % fprintf('whats here BoT: %s\n',tline);
            % pause(235423);

            while ischar(tline)
                via_num = via_num + 1;
                [store rem] = strtok(tline);
                if (isempty(rem) == 1)
                    break;
                else
                    for i=1:1:1
                        [store rem] = strtok(rem);
                    end
                    % fprintf('whats here: %s\n',store);
                    % pause(235423);
                    via_name = store;

                    for i=1:1:1
                        [store rem] = strtok(rem);
                    end
                    % fprintf('whats here: %s\n',store);
                    % pause(235423);
                    via_cnt = str2num(store);
                    
                    [store rem] = strtok(via_name,'_');
                    % for i=1:1:5
                    %     [store rem] = strtok(rem,'_');
                    % end
                    % fprintf('whats here: %s\n',store);
                    % pause(235423);

                    % do search
                    via_notfound_1 = 1;
                    while (via_notfound_1 == 1)
                        % fprintf('via nf top:%i\n',via_notfound_1);
                        [store rem] = strtok(rem,'_');
                        for i=1:1:height(via_resistance_data)
                            if (strcmp(store,via_resistance_data{i,1}) == 1)
                                via_notfound_1 = 0;
                                via_metal{via_num,1} = via_name;
                                via_metal{via_num,2} = store;
                                via_metal{via_num,3} = via_cnt;
                                break;
                            else
                                via_notfound_1 = 1;
                            end
                        end
                        if (via_notfound_1 == 1 && (isempty(store) || isempty(rem)))
                            break;
                        end
                        % fprintf('via nf bot:%i\n',via_notfound_1);
                        % pause(2352356);
                    end
                    % fprintf('via nf super bot:%i\n',via_notfound_1);
                    % pause(235235);

                    if (via_notfound_1 == 0)
                        resistance_targ = via_resistance_data{i,5};
                        tcr_orig = via_resistance_data{i,7};
                        % fprintf('res_targ: %d - tcl: %d\n',via_resistance_data{i,5},via_resistance_data{i,7});
                        via_metal{via_num,4} = resistivity_calc(resistance_targ,via_metal{via_num,3},tcr_orig,tempe,tempe_ref);
                    else
                        fprintf('Via resistance data is not right!\nConsult the documentation of this library\n')
                    end
                    tline = fgetl(fid);
                end
            end
        else
            fprintf('Via layer report is not right!\nPlease regenerate data using this command:\n>> pdi report_design -summary\n');
        end

        % via_metal

        total_resistivity = sum([layer_metal{:,3}]) + sum([via_metal{:,4}]);
        % fprintf('ToT:%.3f\n',total_resistivity);
        % resistivity = total_resistivity/(net_count*wire_len_per_net) * wireload;
        fclose(fid);
    end
    resistivity = total_resistivity/(net_count*wire_len_per_net) * wireload;
end

function resistivity = resistivity_calc(resistance_targ_calc,value,tcr_calc,tempe_calc,tempe_ref_calc)
    resistivity = resistance_targ_calc*value*(1+(tcr_calc/100)*(tempe_calc-tempe_ref_calc));
    % fprintf('tcr:%d\n', tcr_calc);
    fprintf('resistivity:%d\n', resistivity);
end