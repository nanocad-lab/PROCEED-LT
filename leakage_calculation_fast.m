fprintf('Partial Delay and Power Sim NEW (FAST)\n');

libname = ['YOUR_LIBRARY_NAME'];

% Begin to write the netlist
fid=fopen(spname,'w');
fprintf(fid,'\n.lib ''%s'' TT\n',libname);
fprintf(fid,'\n.option converge=1 cptime=%i\n', hspice_time_limit);
% Set the temperature
if (temp)
    fprintf(fid,'.temp %d\n', temp_v);
end

for i_stage = 1:n_bins
    n=2*optimizing_stage(i_stage);
    % the number of partial term for this bin
    n_partial=1+n+flag_vdd+flag_vtn+flag_vtp+(n+flag_vdd+flag_vtn+flag_vtp+1)*(n+flag_vdd+flag_vtn+flag_vtp)/2;
    f2=0;
    % You don't need to choose the model in this experiment
    pulse=5e9*Delay_i(i_stage);
    clock1=10e9*Delay_i(i_stage);
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xc(1,i_stage); % Xc is the starting point for the optimization
        Vtn = Xc(2,i_stage);
        Vtp = Xc(3,i_stage);
    end
    
    % Here Calculate the Leakage Power
    fprintf(fid,'.param sigma=0 l=%d clock_b%i_leak=%4.4fn cl=%d vd_b%i_leak=%d Vtn_b%i_leak=%d Vtp_b%i_leak=%d rl=%d pulse_b%i_leak=%4.4fn fanout=%3.2f Cout=%d ',...
                                    L, i_stage, clock1, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, i_stage, pulse, fanout-1, Cout);
    for width_ind = 1:n
        fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(X(3+width_ind,i_stage)/width_eff_mult));
    end
    fprintf(fid,'\n'); 
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%
    write_derivation7_fast;
    %%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(spname, 'a');

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
        for i = 1:1:n
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
    fprintf(fid,'.DC sweep data=cv_leak_b%i\n', i_stage);
    fprintf(fid,'.option ingold=1 brief=1 \n');
    for leak_iter = 1:1:4
        fprintf(fid,'.print par''-V(vdd_b%i_leak%i)*I(Vd_b%i_leak%i)''\n', i_stage, leak_iter, i_stage, leak_iter);
    end
   
end

fprintf(fid,'.end');
fclose(fid);

unix(hspicerun);
% pause(10);

% To be continued
% Read the Delay and Average Power simulation results
fid=fopen(resultname, 'r');

% Read the Leakage Power simulation results
fid=fopen(logname,'r');



