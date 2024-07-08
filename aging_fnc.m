% AGING CALCUL WITH LOOKUP TABLE
% SUPPORTED TEMPERATURE: -196C, 25C, 85C

function [vtn,vtp] = aging_fnc(tempe,vd,vtn_inp,vtp_inp)
    vtn = 0;
    vtp = 0;
    input_file_folder = 'input_files';
    aging_csv = ['aging_',num2str(tempe),'.csv'];
    aging_data = readmatrix(fullfile(input_file_folder,aging_csv),'Range',1);

    valid_lookup = 0;
    tempe_chk = 0;
    check_col = 1;
    % CONFIRM tempe from the first row

    if (aging_data(1,check_col) == tempe)
        tempe_chk = 1;
        check_col = check_col + 1; %2
    else
        fprintf('Temperature is wrong!\n Lookup temperature is %i. VTH will be set to 0\n',str2num(aging_data(1,check_col)));
    end

    if (tempe_chk == 1)
        % check vtn
        vtn_val = aging_data(1,check_col);
        check_col = check_col + 1; %3
        
        % check vtp
        vtp_val = aging_data(1,check_col);
        check_col = check_col + 1; %4

        % check vdmax
        vd_max_val = aging_data(1,check_col);       % 1V
        % check_col = check_col + 1; %5

        % notify valid lookup
        valid_lookup = 1;
    end
    %------ LOOKUP CHECK END -------

    if (valid_lookup == 0)
        fprintf('Lookup error. Check values on the first row.\n')
    else
        % Check VTH from lookup table
        num_row_vtn = length(aging_data(4,:))-1;
        num_row_vtp = length(aging_data(5,:))-1;
        vdd_min_val = min(vtn_val,vtp_val);
        
        vdd_n = (vd - vtn_inp);
        vdd_p = (vd - vtp_inp);

        if (vdd_n < vtn_val)
            vtn = vtn_inp;
            % fprintf('cs1\n');
        elseif (vdd_n > vd_max_val)
            vtn = vtn_inp + (aging_data(4,num_row_vtn));
            % fprintf('cs2\n');
        else
            idx = int16((vdd_n-vdd_min_val)*1000+1);
            vtn = vtn_inp + (aging_data(4,idx));
            % fprintf('Vtn_shift:%.12f\n',aging_data(4,idx));
            % fprintf('cs3\n');
            % fprintf('vtn idx:%.4f - ',idx);
        end
        % fprintf('VTN_inp:%4f - vdd_n:%.4f - vtn:%.4f\n',vtn_inp, vdd_n, vtn);

        if (vdd_p < vtp_val)
            vtp = vtp_inp;
        elseif (vdd_p > vd_max_val)
            vtp = vtp_inp + (aging_data(5,num_row_vtp));
        else
            idx = int16((vdd_p-vdd_min_val)*1000+1);
            vtp = vtp_inp + (aging_data(5,idx));
            % fprintf('Vtp_shift:%.12f\n',aging_data(5,idx));
            % fprintf('vtp idx:%.4f -',idx);
        end
        % fprintf('VTP_inp:%4f - vdd_p: %.4f - vtp:%.4f\n',vtp_inp, vdd_p,vtp);
    end
    
end