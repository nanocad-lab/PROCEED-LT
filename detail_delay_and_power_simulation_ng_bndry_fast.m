fprintf('Detail Delay and Power Simulation for Multiple bins - NG\n')
% This file will generate a netlist including all bins of circuit that need to be simulated
% And run the simulation then return the delay and power of the circuit
% ============================================================
% Input is the 'X': vdd, vtn, vtp, width data of all bins of circuit.
% Output is the delay and power of the all bins of circuit
% ============================================================

% include the library
libname = ['YOUR_LIBRARY_NAME'];

% ============================================================
% Netlist Generation for Leakage Power Calculation
% ============================================================
leak_arr = zeros(4,n_bins);
% generate the netlist
fid = fopen(spname,'w');
fprintf(fid,'\n.lib ''%s'' TT\n',libname);
fprintf(fid,'.option converge=1\n');

% temperature parameter
if (temp)   
    fprintf(fid,'.temp %d\n', temp_v);      
end


% write the netlist for every bin
for i_stage = 1:n_bins
    n = 2 * optimizing_stage(i_stage);
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        % Extract the vdd, vtn, vtp, width data of the bin
        Vd = Xkc(1,i_stage);
        Vtn = Xkc(2,i_stage);
        Vtp = Xkc(3,i_stage);
    end

    fprintf(fid,'\n\n.param sigma=0 l=%d pulse_b%i_leak=%4.4fn clock_b%i_leak=%4.4fn cl=%d vd_b%i_leak=%d Vtn_b%i_leak=%d Vtp_b%i_leak=%d rl=%d fanout=%3.2f Cout=%d ', ...
                                        L, i_stage, Delay1*1e9*5, i_stage, Delay1*2e9*5, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, fanout-1, Cout);
    for width_ind = 1:i_stage*2
        fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(Xkc(3+width_ind,i_stage)/width_eff_mult));
    end
    fprintf(fid,'\n');

    % NAND Gate Netlist FOR LEAK TEST
    fprintf(fid,['.subckt inv1_b%i_leak 1 2 3 4\n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_leak p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 4 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_leak p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 4 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_leak p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_leak p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n'], i_stage, i_stage, i_stage, i_stage, i_stage);
    % INV Gate Netlistzc FOR LEAK TEST
    fprintf(fid,['.subckt inv2_b%i_leak 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_leak p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_leak p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n'], i_stage, i_stage, i_stage);
    fprintf(fid, '\n');
    
    % write the netlist for the bin
    for leak_iter = 1:1:4
        % Specify the connection of the gates
        for i = 1:1:i_stage*2
            if (mod(i,2))
                fprintf(fid, 'xinv_b%i_%i_leak%i in_b%i_1_leak%i %i_b%i_leak%i vdd_b%i_leak%i in_b%i_2_leak%i inv1_b%i_leak w=w_b%i_%i\n', i_stage, i, leak_iter, i_stage, leak_iter, i+1, i_stage, leak_iter, i_stage, leak_iter, i_stage, leak_iter, i_stage, i_stage, i);
            else
                fprintf(fid, 'xinv_b%i_%i_leak%i %i_b%i_leak%i %i_b%i_leak%i vdd_b%i_leak%i inv2_b%i_leak w=w_b%i_%i\n', i_stage, i, leak_iter, i, i_stage, leak_iter, 100+i-1, i_stage, leak_iter, i_stage, leak_iter, i_stage, i_stage, i);
            end
        end

        % Specify all situations for leakage power calculation
        if (leak_iter == 1) % 00
            fprintf(fid,'Vd_b%i_leak%i vdd_b%i_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
            fprintf(fid,'Vin_b%i_1_leak%i in_b%i_1_leak%i 0 DC 0\n', i_stage, leak_iter, i_stage, leak_iter);
            fprintf(fid,'Vin_b%i_2_leak%i in_b%i_2_leak%i 0 DC 0\n', i_stage, leak_iter, i_stage, leak_iter);
        elseif (leak_iter == 2) % 01
            fprintf(fid,'Vd_b%i_leak%i vdd_b%i_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
            fprintf(fid,'Vin_b%i_1_leak%i in_b%i_1_leak%i 0 DC 0\n', i_stage, leak_iter, i_stage, leak_iter);
            fprintf(fid,'Vin_b%i_2_leak%i in_b%i_2_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
        elseif (leak_iter == 3) % 10
            fprintf(fid,'Vd_b%i_leak%i vdd_b%i_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
            fprintf(fid,'Vin_b%i_1_leak%i in_b%i_1_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
            fprintf(fid,'Vin_b%i_2_leak%i in_b%i_2_leak%i 0 DC 0\n', i_stage, leak_iter, i_stage, leak_iter);
        elseif (leak_iter == 4) %11
            fprintf(fid,'Vd_b%i_leak%i vdd_b%i_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
            fprintf(fid,'Vin_b%i_1_leak%i in_b%i_1_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
            fprintf(fid,'Vin_b%i_2_leak%i in_b%i_2_leak%i 0 DC vd_b%i_leak\n', i_stage, leak_iter, i_stage, leak_iter, i_stage);
        end
    end
    fprintf(fid,'.DC sweep clock_b%i_leak start=1n stop=2n step=1n\n', i_stage);
    fprintf(fid,'.option ingold=1 brief=1 \n');
    for leak_iter = 1:1:4
        fprintf(fid,'.print par''V(vdd_b%i_leak%i)*I(Vd_b%i_leak%i)''\n', i_stage, leak_iter, i_stage, leak_iter);
    end
    
end



% ============================================================
% Netlist Generation for Delay and Power Calculation
% ============================================================

% ============================================================
% ============= Aging Netlist Generation =============
% ============================================================
% Calculate the power and latency of the circuit separately
if (aging || vt_var_mode ~= 0 || flag_IR_drop)
    fprintf(fid, '\n\n\n');
    for i_stage = 1:n_bins
        % ===================================================================================
        % ================== Frist calculate power (this is average power) ==================
        % Should use non-aging vt to calculate the power here for the worst case
        if fixed_V
            Vd = 0.7;
            Vtn = 0;
            Vtp = 0;
        else
            % Extract the vdd, vtn, vtp, width data of the bin
            Vd = Xkc(1,i_stage);
            Vtn = Xkc(2,i_stage);
            Vtp = Xkc(3,i_stage);
        end
        if(vt_var_mode ~= 0)
            % fprintf('Doing negative vth variation now...\n');
            % fprintf('Before variation: Vtn = %.4f, Vtp = %.4f\n', Vtn, Vtp);
            Vtn = Vtn - vtn_var_val - vtn_percent * (Vtn0 + Vtn);   % Here Vtn is the shift of NFET Vth
            Vtp = Vtp - vtp_var_val - vtp_percent * (Vtp0 + Vtp);   % Here Vtp is the shift of PFET Vth
            % fprintf('After variation: Vtn = %.4f, Vtp = %.4f\n', Vtn, Vtp);
        end
    
        fprintf(fid,'\n\n.param sigma=0 l=%d pulse_b%i_pwr=%4.4fn clock_b%i_pwr=%4.4fn cl=%d vd_b%i_pwr=%d Vtn_b%i_pwr=%d Vtp_b%i_pwr=%d rl=%d fanout=%3.2f Cout=%d ', ...
                                            L, i_stage, Delay1*1e9*5, i_stage, Delay1*2e9*5, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, fanout-1, Cout);
        for width_ind = 1:i_stage*2
            fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(Xkc(3+width_ind,i_stage)/width_eff_mult));
        end
        fprintf(fid,'\n');
    
        % NAND Gate Netlist 
        fprintf(fid,['.subckt inv1_b%i_pwr 1 2 3 \n', ...
        'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_pwr p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_pwr p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_pwr p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_pwr p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        '.ends\n'], i_stage, i_stage, i_stage, i_stage, i_stage);
        % INV Gate Netlistzc 
        fprintf(fid,['.subckt inv2_b%i_pwr 1 2 3 \n', ...
        'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_pwr p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_pwr p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Rl 6 2 rl\n', ...
        'Cl 2 0 cl\n', ...
        '.ends\n'], i_stage, i_stage, i_stage);
        fprintf(fid, '\n');

        % write the netlist for the driver
        fprintf(fid, 'xinv100_b%i_pwr 100_b%i_pwr 1_b%i_pwr vdr_b%i_pwr inv2_b%i_pwr w=w_b%i_1\n', i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);
        for i=1:1:i_stage*2
            if(mod(i,2))
                fprintf(fid,'xinv_b%i_%i_pwr %i_b%i_pwr %i_b%i_pwr vdd_b%i_pwr inv1_b%i_pwr w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
            else
                fprintf(fid,'xinv_b%i_%i_pwr %i_b%i_pwr %i_b%i_pwr vdd_b%i_pwr inv2_b%i_pwr w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
                if(i<i_stage*2)
                    fprintf(fid,'xinv_b%i_%i_pwr %i_b%i_pwr %i_b%i_pwr vdr_b%i_pwr inv1_b%i_pwr w=''fanout*w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1);
                else
                    fprintf(fid,'xinv_b%i_%i_pwr %i_b%i_pwr %i_b%i_pwr vdr_b%i_pwr inv1_b%i_pwr w=''fanout*w_b%i_%i+w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1);
                end
            end
        end

        if (flag_var)
            fprintf(fid,'.subckt inv1_s_b%i_pwr 1 2 3 dvt=0 w_b%i=1u\n', i_stage);
            fprintf(fid,'Xn1 7 4 0 0 d_nfet W=w_b%i L=l\n', i_stage);
            fprintf(fid,'Xn2 2 8 7 0 d_nfet W=w_b%i L=l\n');
            fprintf(fid,'Xp1 2 9 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Xp2 5 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
            fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
            fprintf(fid,'Vtno 8 3 ''dvt+vshift'' \n');
            fprintf(fid,'Vtpo 9 3 ''-dvt-vshift''\n.ends\n');            
            fprintf(fid,'.subckt inv2_s_b%i_pwr 1 2 3 dvt=0 w_b%i=1u\n');
            fprintf(fid,'Xn 6 4 0 0 d_nfet W=w_b%i L=l\n');
            fprintf(fid,'Xp 6 5 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
            fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
            fprintf(fid,'Rl 6 2 rl\nCl 2 0 cl\n.ends\n');

            fprintf(fid,'xinv100_s_b%i_pwr 100_s_b%i_pwr 1_s_b%i_pwr vdr_s_b%i_pwr inv2_s_b%i_pwr w=''w_b%i_1-%4.4f'' dvt=%3.3f\n', ...
                                        i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, dw1, dvt);
            for i=1:1:i_stage*2
                if(mod(i,2))
                    fprintf(fid,'xinv%i_s_b%i_pwr %i_s_b%i_pwr %i_s_b%i_pwr vdd_s_b%i_pwr inv1_s_b%i_pwr w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                        i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                else
                    fprintf(fid,'xinv%i_s_b%i_pwr %i_s_b%i_pwr %i_s_b%i_pwr vdd_s_b%i_pwr inv2_s_b%i_pwr w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                        i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                    if(i<i_stage*2)
                        fprintf(fid,'xinv%i_s_b%i_pwr %i_s_b%i_pwr %i_s_b%i_pwr vdr_s_b%i_pwr inv1_s_b%i_pwr w=''fanout*w_b%i_%i-fanout*%4.4f'' dvt=%3.3f\n', ...
                                        100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, dw1, dvt);
                    else
                        fprintf(fid,'xinv%i_s_b%i_pwr %i_s_b%i_pwr %i_S_b%i_pwr vdr_s_b%i_pwr inv1_s_b%i_pwr w=''fanout*w_b%i_%i+w_b%i_%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n', ...
                                        100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1, dw1, dw1, dvt);
                    end
                    
                end
            end
            fprintf(fid,'Vdr_s_b%i_pwr vdr_s_b%i_pwr 0 DC=''vd_b%i_pwr*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
            fprintf(fid,'Vd_s_b%i_pwr vdd_s_b%i_pwr 0 DC=''vd_b%i_pwr*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
            fprintf(fid,'Vpulse_s_b%i_pwr 100_s_b%i_pwr 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock_b%i_pwr/2'',clock_b%i_pwr)\n', ...
                                i_stage, i_stage, p_slow, i_stage, i_stage);
        end
        
        
        fprintf(fid,['Vdr_b%i_pwr vdr_b%i_pwr 0 DC vd_b%i_pwr\n', ...
                    'Vd_b%i_pwr vdd_b%i_pwr 0 DC vd_b%i_pwr\n', ...
                    'Vpulse_b%i_pwr 100_b%i_pwr 0 PULSE(0,vd_b%i_pwr,10n,10p,10p,''clock_b%i_pwr/2'',clock_b%i_pwr)\n', ...
                    '.tran ''clock_b%i_pwr/10000'' ''clock_b%i_pwr*5+10n'' \n'], i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);

        if (flag_var)
            % DO NOTHING
        else
            % Only measure the power here
            fprintf(fid,'.measure tran avgpower_b%i avg par''(-1)*abs(p(Vd_b%i_pwr))'' FROM=''clock_b%i_pwr*2+10n'' TO=''clock_b%i_pwr*5+10n''\n', i_stage, i_stage, i_stage, i_stage);
        end

        fprintf(fid,'.ic ');
        fprintf(fid,'v(100_b%i_pwr)=0 ', i_stage);
        for i=1:i_stage*2
            if(mod(i,2)==0)
                fprintf(fid,'v(%i_b%i_pwr)=0 ', i, i_stage);
            else
                fprintf(fid,'v(%i_b%i_pwr)=''vd_b%i_pwr'' ', i, i_stage, i_stage);
            end
        end
        if (flag_var)
            fprintf(fid,'v(100_s_b%i_pwr)=0 ', i_stage);
            for i=1:i_stage*2
                if(mod(i,2)==0)
                    fprintf(fid,'v(%i_s_b%i_pwr)=0 ', i, i_stage);
                else
                    fprintf(fid,'v(%i_s_b%i_pwr)=''vd*%2.2f'' ', i, i_stage, p_slow);
                end
            end
        end


        % ==============================================================
        % ================== Then calculate the delay ==================
        if (aging)
            [Vtn_age, Vtp_age] = aging_fnc(temp_v, Vd, Vtn, Vtp);
            % fprintf('\nBefore Aging: Vtn:%d, Vtp:%d; \n After Aging: Vtn_age:%d, Vtp_age:%d\n', Vtn, Vtp, Vtn_age, Vtp_age);
            Vtn = Vtn_age;
            Vtp = Vtp_age;
        end
        if (vt_var_mode ~= 0)
            % fprintf('Doing positive vth variation now...\n');
            % fprintf('\nBefore adding variation, Vtn:%d, Vtp:%d\n', Vtn, Vtp);
            % Add Vth variation impact
            Vtn = Vtn + vtn_var_val + vtn_percent * (Vtn0 + Vtn);   % Here Vtn is the shift of NFET Vth
            Vtp = Vtp + vtp_var_val + vtp_percent * (Vtp0 + Vtp);   % Here Vtp is the shift of PFET Vth
            % fprintf('\nAfter adding variation, Vtn:%d, Vtp:%d\n', Vtn, Vtp);
        end
        % if (flag_IR_drop)
        %     Vd = Vd - 13.435 * ;
        % end

        % Reset the Vtn and Vtp parameters for the aging calculation
        fprintf(fid,'\n\n.param sigma=0 l=%d pulse_b%i_dly=%4.4fn clock_b%i_dly=%4.4fn cl=%d vd_b%i_dly=%d Vtn_b%i_dly=%d Vtp_b%i_dly=%d rl=%d fanout=%3.2f Cout=%d ', ...
                                            L, i_stage, Delay1*1e9*5, i_stage, Delay1*2e9*5, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, fanout-1, Cout);
        
        for width_ind = 1:i_stage*2
            fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(Xkc(3+width_ind,i_stage)/width_eff_mult));
        end
        fprintf(fid,'\n');
    
        % NAND Gate Netlist
        fprintf(fid,['.subckt inv1_b%i_dly 1 2 3\n', ...
        'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_dly p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_dly p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_dly p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_dly p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        '.ends\n'], i_stage, i_stage, i_stage, i_stage, i_stage);
        % INV Gate Netlistzc
        fprintf(fid,['.subckt inv2_b%i_dly 1 2 3 \n', ...
        'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_dly p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_dly p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Rl 6 2 rl\n', ...
        'Cl 2 0 cl\n', ...
        '.ends\n'], i_stage, i_stage, i_stage);
        fprintf(fid, '\n');

        % write the netlist for the driver
        fprintf(fid, 'xinv100_b%i_dly 100_b%i_dly 1_b%i_dly vdr_b%i_dly inv2_b%i_dly w=w_b%i_1\n', i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);
        for i=1:1:i_stage*2
            if(mod(i,2))
                fprintf(fid,'xinv_b%i_%i_dly %i_b%i_dly %i_b%i_dly vdd_b%i_dly inv1_b%i_dly w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
            else
                fprintf(fid,'xinv_b%i_%i_dly %i_b%i_dly %i_b%i_dly vdd_b%i_dly inv2_b%i_dly w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
                if(i<i_stage*2)
                    fprintf(fid,'xinv_b%i_%i_dly %i_b%i_dly %i_b%i_dly vdr_b%i_dly inv1_b%i_dly w=''fanout*w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1);
                else
                    fprintf(fid,'xinv_b%i_%i_dly %i_b%i_dly %i_b%i_dly vdr_b%i_dly inv1_b%i_dly w=''fanout*w_b%i_%i+w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1);
                end
            end
        end

        if (flag_var)
            fprintf(fid,'.subckt inv1_s_b%i_dly 1 2 3 dvt=0 w_b%i=1u\n', i_stage);
            fprintf(fid,'Xn1 7 4 0 0 d_nfet W=w_b%i L=l\n', i_stage);
            fprintf(fid,'Xn2 2 8 7 0 d_nfet W=w_b%i L=l\n');
            fprintf(fid,'Xp1 2 9 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Xp2 5 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
            fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
            fprintf(fid,'Vtno 8 3 ''dvt+vshift'' \n');
            fprintf(fid,'Vtpo 9 3 ''-dvt-vshift''\n.ends\n');            
            fprintf(fid,'.subckt inv2_s_b%i_dly 1 2 3 dvt=0 w_b%i=1u\n');
            fprintf(fid,'Xn 6 4 0 0 d_nfet W=w_b%i L=l\n');
            fprintf(fid,'Xp 6 5 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
            fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
            fprintf(fid,'Rl 6 2 rl\nCl 2 0 cl\n.ends\n');

            fprintf(fid,'xinv100_s_b%i_dly 100_s_b%i_dly 1_s_b%i_dly vdr_s_b%i_dly inv2_s_b%i_dly w=''w_b%i_1-%4.4f'' dvt=%3.3f\n', ...
                                        i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, dw1, dvt);
            for i=1:1:i_stage*2
                if(mod(i,2))
                    fprintf(fid,'xinv%i_s_b%i_dly %i_s_b%i_dly %i_s_b%i_dly vdd_s_b%i_dly inv1_s_b%i_dly w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                        i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                else
                    fprintf(fid,'xinv%i_s_b%i_dly %i_s_b%i_dly %i_s_b%i_dly vdd_s_b%i_dly inv2_s_b%i_dly w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                        i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                    if(i<i_stage*2)
                        fprintf(fid,'xinv%i_s_b%i_dly %i_s_b%i_dly %i_s_b%i_dly vdr_s_b%i_dly inv1_s_b%i_dly w=''fanout*w_b%i_%i-fanout*%4.4f'' dvt=%3.3f\n', ...
                                        100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, dw1, dvt);
                    else
                        fprintf(fid,'xinv%i_s_b%i_dly %i_s_b%i_dly %i_S_b%i_dly vdr_s_b%i_dly inv1_s_b%i_dly w=''fanout*w_b%i_%i+w_b%i_%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n', ...
                                        100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1, dw1, dw1, dvt);
                    end
                    
                end
            end
            fprintf(fid,'Vdr_s_b%i_dly vdr_s_b%i_dly 0 DC=''vd_b%i_dly*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
            fprintf(fid,'Vd_s_b%i_dly vdd_s_b%i_dly 0 DC=''vd_b%i_dly*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
            fprintf(fid,'Vpulse_s_b%i_dly 100_s_b%i_dly 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock_b%i_dly/2'',clock_b%i_dly)\n', i_stage, i_stage, p_slow, i_stage, i_stage);
        end
        
        
        fprintf(fid,['Vdr_b%i_dly vdr_b%i_dly 0 DC vd_b%i_dly\n',  ...
                    'Vd_b%i_dly vdd_b%i_dly 0 DC vd_b%i_dly\n',  ...
                    'Vpulse_b%i_dly 100_b%i_dly 0 PULSE(0,vd_b%i_dly,10n,10p,10p,''clock_b%i_dly/2'',clock_b%i_dly)\n', ...
                    '.tran ''clock_b%i_dly/10000'' ''clock_b%i_dly*5+10n'' \n'], i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);

        if (flag_var)
            fprintf(fid,'.measure tran t1_b%i trig V(1_s_b%i_dly) val=''%2.2f*vd_b%i_dly/2'' TD=''10n+clock_b%i_dly/4'' rise=2 targ V(%i_s_b%i_dly) val=''%2.2f*vd_b%i_dly/2'' TD=''10n+clock_b%i_dly/4'' %s=2\n', ...
                                            i_stage, i_stage, p_slow, i_stage, i_stage, i_stage*2+1, i_stage, p_slow, i_stage, i_stage, edge);
        else
            % Only measure the delay here
            fprintf(fid,'.measure avgdelay_b%i param=''(t1_b%i+t2_b%i)*0.5''\n', i_stage, i_stage, i_stage);
            fprintf(fid,'.measure tran t1_b%i trig V(1_b%i_dly) val=''vd_b%i_dly*0.5'' TD=0 rise=3 targ V(%i_b%i_dly) val=''vd_b%i_dly*0.5'' TD=0 rise=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
            fprintf(fid,'.measure tran t2_b%i trig V(1_b%i_dly) val=''vd_b%i_dly*0.5'' TD=0 fall=3 targ V(%i_b%i_dly) val=''vd_b%i_dly*0.5'' TD=0 fall=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
        end

        fprintf(fid,'.ic ');
        fprintf(fid,'v(100_b%i_dly)=0 ', i_stage);
        for i=1:i_stage*2
            if(mod(i,2)==0)
                fprintf(fid,'v(%i_b%i_dly)=0 ', i, i_stage);
            else
                fprintf(fid,'v(%i_b%i_dly)=''vd_b%i_dly'' ', i, i_stage, i_stage);
            end
        end
        if (flag_var)
            fprintf(fid,'v(100_s_b%i_dly)=0 ', i_stage);
            for i=1:n
                if(mod(i,2)==0)
                    fprintf(fid,'v(%i_s_b%i_dly)=0 ', i, i_stage);
                else
                    fprintf(fid,'v(%i_s_b%i_dly)=''vd_b%i_dly*%2.2f'' ', i, i_stage, i_stage, p_slow);
                end
            end
        end

    end


% ============================================================
% ============= non-Aging Netlist Generation =============
% ============================================================
% Calculate the power and latency of the circuit together
else
    fprintf(fid, '\n\n\n');
    for i_stage = 1:n_bins
        if fixed_V
            Vd = 0.7;
            Vtn = 0;
            Vtp = 0;
        else
            % Extract the vdd, vtn, vtp, width data of the bin
            Vd = Xkc(1,i_stage);
            Vtn = Xkc(2,i_stage);
            Vtp = Xkc(3,i_stage);
        end
    
        fprintf(fid,'\n\n.param sigma=0 l=%d pulse_b%i=%4.4fn clock_b%i=%4.4fn cl=%d vd_b%i=%d Vtn_b%i=%d Vtp_b%i=%d rl=%d fanout=%3.2f Cout=%d ', ...
                                            L, i_stage, Delay1*1e9*5, i_stage, Delay1*2e9*5, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, fanout-1, Cout);
        for width_ind = 1:i_stage*2
            fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(Xkc(3+width_ind,i_stage)/width_eff_mult));
        end
        fprintf(fid,'\n');

        % NAND Gate Netlist
        fprintf(fid,['.subckt inv1_b%i 1 2 3 \n', ...
        'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        '.ends\n'], i_stage, i_stage, i_stage, i_stage, i_stage);
        % INV Gate Netlistzc
        fprintf(fid,['.subckt inv2_b%i 1 2 3 \n', ...
        'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Rl 6 2 rl\n', ...
        'Cl 2 0 cl\n', ...
        '.ends\n'], i_stage, i_stage, i_stage);
        fprintf(fid, '\n');

        % write the netlist for the driver
        fprintf(fid, 'xinv100_b%i 100_b%i 1_b%i vdr_b%i inv2_b%i w=w_b%i_1\n', i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);
        for i=1:1:i_stage*2
            if(mod(i,2))
                fprintf(fid,'xinv_b%i_%i %i_b%i %i_b%i vdd_b%i inv1_b%i w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
            else
                fprintf(fid,'xinv_b%i_%i %i_b%i %i_b%i vdd_b%i inv2_b%i w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
                if(i<i_stage*2)
                    fprintf(fid,'xinv_b%i_%i %i_b%i %i_b%i vdr_b%i inv1_b%i w=''fanout*w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1);
                else
                    fprintf(fid,'xinv_b%i_%i %i_b%i %i_b%i vdr_b%i inv1_b%i w=''fanout*w_b%i_%i+w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1);
                end
            end
        end

        if (flag_var)
            fprintf(fid,'.subckt inv1_s_b%i 1 2 3 dvt=0 w_b%i=1u\n', i_stage);
            fprintf(fid,'Xn1 7 4 0 0 d_nfet W=w_b%i L=l\n', i_stage);
            fprintf(fid,'Xn2 2 8 7 0 d_nfet W=w_b%i L=l\n');
            fprintf(fid,'Xp1 2 9 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Xp2 5 3 3 d_pfet W=w_b%i L=l\n');
            fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
            fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
            fprintf(fid,'Vtno 8 3 ''dvt+vshift'' \n');
            fprintf(fid,'Vtpo 9 3 ''-dvt-vshift''\n.ends\n');            
            fprintf(fid,'.subckt inv2_s_b%i 1 2 3 dvt=0 w_b%i=1u\n', i_stage, i_stage);
            fprintf(fid,'Xn 6 4 0 0 d_nfet W=w_b%i L=l\n', i_stage);
            fprintf(fid,'Xp 6 5 3 3 d_pfet W=w_b%i L=l\n', i_stage);
            fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
            fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
            fprintf(fid,'Rl 6 2 rl\nCl 2 0 cl\n.ends\n');

            fprintf(fid,'xinv100_s_b%i 100_s_b%i 1_s_b%i vdr_s_b%i inv2_s_b%i w=''w_b%i_1-%4.4f'' dvt=%3.3f\n', ...
                                        i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, dw1, dvt);
            for i=1:1:i_stage*2
                if(mod(i,2))
                    fprintf(fid,'xinv%i_s_b%i %i_s_b%i %i_s_b%i vdd_s_b%i inv1_s_b%i w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                        i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                else
                    fprintf(fid,'xinv%i_s_b%i %i_s_b%i %i_s_b%i vdd_s_b%i inv2_s_b%i w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                        i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                    if(i<i_stage*2)
                        fprintf(fid,'xinv%i_s_b%i %i_s_b%i %i_s_b%i vdr_s_b%i inv1_s_b%i w=''fanout*w_b%i_%i-fanout*%4.4f'' dvt=%3.3f\n', ...
                                        100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, dw1, dvt);
                    else
                        fprintf(fid,'xinv%i_s_b%i %i_s_b%i %i_S_b%i vdr_s_b%i inv1_s_b%i w=''fanout*w_b%i_%i+w_b%i_%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n', ...
                                        100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1, dw1, dw1, dvt);
                    end
                    
                end
            end
            fprintf(fid,'Vdr_s_b%i vdr_s_b%i 0 DC=''vd_b%i*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
            fprintf(fid,'Vd_s_b%i vdd_s_b%i 0 DC=''vd_b%i*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
            fprintf(fid,'Vpulse_s_b%i 100_s_b%i 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock_b%i/2'',clock_b%i)\n', i_stage, i_stage, p_slow, i_stage, i_stage);
        end
        
        
        fprintf(fid,['Vdr_b%i vdr_b%i 0 DC vd_b%i\n', ...
                    'Vd_b%i vdd_b%i 0 DC vd_b%i\n',  ...
                    'Vpulse_b%i 100_b%i 0 PULSE(0,vd_b%i,10n,10p,10p,''clock_b%i/2'',clock_b%i)\n', ...
                    '.tran ''clock_b%i/10000'' ''clock_b%i*5+10n'' \n'], ...
                    i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);

        if (flag_var)
            fprintf(fid,'.measure tran t1_b%i trig V(1_s_b%i) val=''%2.2f*vd_b%i/2'' TD=''10n+clock_b%i/4'' rise=2 targ V(%i_s_b%i) val=''%2.2f*vd_b%i/2'' TD=''10n+clock_b%i/4'' %s=2\n', ...
                                            i_stage, i_stage, p_slow, i_stage, i_stage, i_stage*2+1, i_stage, p_slow, i_stage, i_stage, edge);
        else
            % Only measure the delay here
            fprintf(fid,'.measure tran avgpower_b%i avg par''(-1)*abs(p(Vd_b%i))'' FROM=''clock_b%i*2+10n'' TO=''clock_b%i*5+10n''\n', i_stage, i_stage, i_stage, i_stage);
            fprintf(fid,'.measure avgdelay_b%i param=''(t1_b%i+t2_b%i)*0.5''\n', i_stage, i_stage, i_stage);
            fprintf(fid,'.measure tran t1_b%i trig V(1_b%i) val=''vd_b%i*0.5'' TD=0 rise=3 targ V(%i_b%i) val=''vd_b%i*0.5'' TD=0 rise=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
            fprintf(fid,'.measure tran t2_b%i trig V(1_b%i) val=''vd_b%i*0.5'' TD=0 fall=3 targ V(%i_b%i) val=''vd_b%i*0.5'' TD=0 fall=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
        end

        fprintf(fid,'.ic ');
        fprintf(fid,'v(100_b%i)=0 ', i_stage);
        for i=1:i_stage*2
            if(mod(i,2)==0)
                fprintf(fid,'v(%i_b%i)=0 ', i, i_stage);
            else
                fprintf(fid,'v(%i_b%i)=''vd_b%i'' ', i, i_stage, i_stage);
            end
        end
        if (flag_var)
            fprintf(fid,'v(100_s_b%i)=0 ', i_stage);
            for i=1:n
                if(mod(i,2)==0)
                    fprintf(fid,'v(%i_s_b%i)=0 ', i, i_stage);
                else
                    fprintf(fid,'v(%i_s_b%i)=''vd_b%i*%2.2f'' ', i, i_stage, i_stage, p_slow);
                end
            end
        end

    end

end






% end the netlist
fprintf(fid,'\n\n.end');
% Run the simulation
fclose(fid);
[hspicestatus, cmdout] = unix(hspicerun);
hspice_running = 1;
pause_second = 5;
pause_count = 0;
while(hspice_running == 1)
    pause(pause_second);
    pause_count = pause_count + 1;
    found = check_hspice_conclude(logname, 6, 'hspice job concluded');
    fprintf('Hspice running... Consuming %d sec\n', pause_count*pause_second);
    if(pause_second * pause_count > 1000)
        error('Hspice running time out!\n');
    end
    if(found)
        hspice_running = 0;
        fprintf('Hspice simulation done!\n');
    end
end


% Read the output file and extract the leakage power
leakage2temp=38;
fid=fopen(logname,"r");
f3 = 1; % a flag that will be '1' if reading the leakage power unsuccessfully
str2search = 'v\(vdd_b(\d+)_leak(\d+)\)\*i\(vd_b(\d+)_leak(\d+)\)';      % regular expression to search for the leakage power
while (~feof(fid))
    tline = fgetl(fid);
    [t r] = strtok(tline);
    matches = regexp(t, str2search, 'tokens');
    if (~isempty(matches))
        % found where the leakage power is
        numericValues = str2double(matches{1});
        bin_ind = numericValues(1);
        leak_iter_ind = numericValues(2);
        tline = fgetl(fid);
        [t r] = strtok(tline);
        [x status] = str2num(r);  % 'status' show that if it transform from str to num successfully?
        if (status ~= 0)
            leakage2temp = x;
            f3 = 0;
            leak_arr(leak_iter_ind, bin_ind) = (-1)*abs(leakage2temp);
        end
        numericValues = [];
        if (bin_ind == n_bins && leak_iter_ind == 4)
            break;
        end
    end
end
fclose(fid);


% Calculate the mean value of leakage power of each bin
Leakage2 = mean(leak_arr, 1);


%  Read the output file and extract the delay and average power
fid=fopen(resultname, "r");
Delay2 = [];
Energy2 = [];
if (fid ~= -1)
    for (j=1:(2+n_bins+1+1))
        tline=fgetl(fid);
    end
    for i=1:1:n_bins
        [avgpwr_temp r] = strtok(tline);
        [avgdly_temp r] = strtok(r);
        avgpwr_temp = str2num(avgpwr_temp);
        avgdly_temp = str2num(avgdly_temp);

        if(avgpwr_temp>=0  || avgpwr_temp>Leakage2(i)/2)
            f4=1;
            disp('energy2 >=0 or energy2 < leakage');
            disp(avgpwr_temp);
            disp(Leakage2(i));
        else
            if(avgpwr_temp>Leakage2(i))
                avgpwr_temp=Leakage2(i);
            end
        end
        if(isempty(avgdly_temp))
            disp('delay simulation failed in this boundary allocation scheme (maybe vd too low...)')
            avgdly_temp=38;
        end
        % Here 'Energy2' is the average power
        % But later it will be energy through calculation
        Energy2(i) = avgpwr_temp;
        Delay2(i) = avgdly_temp;

        tline=fgetl(fid);
    end
    
else
    error('Error: cannot open the result file');
end
fclose(fid);

if (f4 == 1)
    error('Error: energy2 >=0 or energy2 < leakage');
end



Dynamic2 = (-Energy2 + Leakage2) / calib_dyn;
Leakage2 = Leakage2 / calib_leak;

Energy2 = Dynamic2 .* (10 * Delay1) * 0.5 * activity;
fprintf('Detail Delay2 for each bin: \n');
disp(Delay2);
fprintf("VDD1: %d, VDD2: %d\n", Xkc(1, min(optimizing_stage)), Xkc(1, max(optimizing_stage)));
fprintf("Vtn1: %d, Vtn2: %d\n", Xkc(2, min(optimizing_stage)), Xkc(2, max(optimizing_stage)));
fprintf("Vtp1: %d, Vtp2: %d\n", Xkc(3, min(optimizing_stage)), Xkc(3, max(optimizing_stage)));


