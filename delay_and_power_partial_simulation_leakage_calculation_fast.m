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
    
    fprintf(fid,'\n.param sigma=0 l=%d clock_b%i_part=%4.4fn cl=%d vd_b%i_part=%d Vtn_b%i_part=%d Vtp_b%i_part=%d rl=%d pulse_b%i_part=%4.4fn fanout=%3.2f Cout=%d ', ...
                                    L, i_stage, clock1, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, i_stage, pulse, fanout-1, Cout);
    for width_ind = 1:n
        fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(X(3+width_ind,i_stage)/width_eff_mult));
    end
    fprintf(fid,'\n');
    fclose(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    write_derivation6_fast;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fid=fopen(spname,'a');

    % NAND Gate Netlist
    fprintf(fid,['.subckt inv1_b%i_part 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_part p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_part p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_part p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_part p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n'], i_stage, i_stage, i_stage, i_stage, i_stage);
    % INV Gate Netlistzc 
    fprintf(fid,['.subckt inv2_b%i_part 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_part p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_part p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n'], i_stage, i_stage, i_stage);
    fprintf(fid, '\n');

    % write the netlist for the driver
    fprintf(fid, 'xinv100_b%i_part 100_b%i_part 1_b%i_part vdr_b%i_part inv2_b%i_part w=w_b%i_1\n', i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);
    for i=1:1:i_stage*2
        if(mod(i,2))
            fprintf(fid,'xinv_b%i_%i_part %i_b%i_part %i_b%i_part vdd_b%i_part inv1_b%i_part w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
        else
            fprintf(fid,'xinv_b%i_%i_part %i_b%i_part %i_b%i_part vdd_b%i_part inv2_b%i_part w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
            if(i<i_stage*2)
                fprintf(fid,'xinv_b%i_%i_part %i_b%i_part %i_b%i_part vdr_b%i_part inv1_b%i_part w=''fanout*w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1);
            else
                fprintf(fid,'xinv_b%i_%i_part %i_b%i_part %i_b%i_part vdr_b%i_part inv1_b%i_part w=''fanout*w_b%i_%i+w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1);
            end
        end
    end

    if (flag_var)
        fprintf(fid,'.subckt inv1_s_b%i_part 1 2 3 dvt=0 w_b%i=1u\n', i_stage);
        fprintf(fid,'Xn1 7 4 0 0 d_nfet W=w_b%i L=l\n', i_stage);
        fprintf(fid,'Xn2 2 8 7 0 d_nfet W=w_b%i L=l\n');
        fprintf(fid,'Xp1 2 9 3 3 d_pfet W=w_b%i L=l\n');
        fprintf(fid,'Xp2 5 3 3 d_pfet W=w_b%i L=l\n');
        fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
        fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
        fprintf(fid,'Vtno 8 3 ''dvt+vshift'' \n');
        fprintf(fid,'Vtpo 9 3 ''-dvt-vshift''\n.ends\n');            
        fprintf(fid,'.subckt inv2_s_b%i_part 1 2 3 dvt=0 w_b%i=1u\n');
        fprintf(fid,'Xn 6 4 0 0 d_nfet W=w_b%i L=l\n');
        fprintf(fid,'Xp 6 5 3 3 d_pfet W=w_b%i L=l\n');
        fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
        fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
        fprintf(fid,'Rl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s_b%i_part 100_s_b%i_part 1_s_b%i_part vdr_s_b%i_part inv2_s_b%i_part w=''w_b%i_1-%4.4f'' dvt=%3.3f\n', ...
                                    i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, dw1, dvt);
        for i=1:1:i_stage*2
            if(mod(i,2))
                fprintf(fid,'xinv%i_s_b%i_part %i_s_b%i_part %i_s_b%i_part vdd_s_b%i_part inv1_s_b%i_part w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                    i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
            else
                fprintf(fid,'xinv%i_s_b%i_part %i_s_b%i_part %i_s_b%i_part vdd_s_b%i_part inv2_s_b%i_part w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                    i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                if(i<i_stage*2)
                    fprintf(fid,'xinv%i_s_b%i_part %i_s_b%i_part %i_s_b%i_part vdr_s_b%i_part inv1_s_b%i_part w=''fanout*w_b%i_%i-fanout*%4.4f'' dvt=%3.3f\n', ...
                                    100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, dw1, dvt);
                else
                    fprintf(fid,'xinv%i_s_b%i_part %i_s_b%i_part %i_S_b%i_part vdr_s_b%i_part inv1_s_b%i_part w=''fanout*w_b%i_%i+w_b%i_%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n', ...
                                    100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1, dw1, dw1, dvt);
                end
                
            end
        end
        fprintf(fid,'Vdr_s_b%i_part vdr_s_b%i_part 0 DC=''vd_b%i_part*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
        fprintf(fid,'Vd_s_b%i_part vdd_s_b%i_part 0 DC=''vd_b%i_part*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
        fprintf(fid,'Vpulse_s_b%i_part 100_s_b%i_part 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock_b%i_part/2'',clock_b%i_part)\n', i_stage, i_stage, p_slow, i_stage, i_stage);
    end
    
    
    fprintf(fid,['Vdr_b%i_part vdr_b%i_part 0 DC vd_b%i_part\n',  ...
                'Vd_b%i_part vdd_b%i_part 0 DC vd_b%i_part\n',  ...
                'Vpulse_b%i_part 100_b%i_part 0 PULSE(0,vd_b%i_part,10n,10p,10p,pulse_b%i_part,clock_b%i_part)\n', ...
                '.tran ''clock_b%i_part/10000'' ''clock_b%i_part*5+10n'' sweep data=cv_b%i\n'], ...
                i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);

    if (flag_var)
        fprintf(fid,'.measure tran t1_b%i trig V(1_s_b%i_part) val=''%2.2f*vd_b%i_part/2'' TD=''10n+pulse_b%i_part*3/4'' rise=2 targ V(%i_s_b%i_part) val=''%2.2f*vd_b%i_part/2'' TD=''10n+pulse_b%i_part*3/4'' %s=2\n', ...
                                        i_stage, i_stage, p_slow, i_stage, i_stage, n+1, i_stage, p_slow, i_stage, i_stage, edge);
    else
        fprintf(fid,'.measure tran avgpower_b%i avg par''(-1)*abs(p(Vd_b%i_part))'' FROM=''clock_b%i_part*2+10n'' TO=''clock_b%i_part*5+10n''\n', i_stage, i_stage, i_stage, i_stage);
        fprintf(fid,'.measure avgdelay_b%i param=''(t1_b%i+t2_b%i)*0.5''\n', i_stage, i_stage, i_stage);
        fprintf(fid,'.measure tran t1_b%i trig V(1_b%i_part) val=''vd_b%i_part*0.5'' TD=0 rise=3 targ V(%i_b%i_part) val=''vd_b%i_part*0.5'' TD=0 rise=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
        fprintf(fid,'.measure tran t2_b%i trig V(1_b%i_part) val=''vd_b%i_part*0.5'' TD=0 fall=3 targ V(%i_b%i_part) val=''vd_b%i_part*0.5'' TD=0 fall=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
    end

    fprintf(fid,'.ic ');
    fprintf(fid,'v(100_b%i_part)=0 ', i_stage);
    for i=1:n+1
        if(mod(i,2)==0)
            fprintf(fid,'v(%i_b%i_part)=0 ', i, i_stage);
        else
            fprintf(fid,'v(%i_b%i_part)=''vd_b%i_part'' ', i, i_stage, i_stage);
        end
    end
    if (flag_var)
        fprintf(fid,'v(100_s_b%i_part)=0 ', i_stage);
        for i=1:n
            if(mod(i,2)==0)
                fprintf(fid,'v(%i_s_b%i_part)=0 ', i, i_stage);
            else
                fprintf(fid,'v(%i_s_b%i_part)=''vd_b%i_part*%2.2f'' ', i, i_stage, i_stage, p_slow);
            end
        end
    end

    fprintf(fid,'\n\n\n');



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



