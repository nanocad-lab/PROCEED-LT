fprintf('Xk Init Delay Sim NEW - NG (FAST)\n');

libname = ['YOUR_LIBRARY_NAME'];

%%end
fid=fopen(spname,'w');
%fprintf(fid,'\n.inc ''%s'' \n',subcktname);
fprintf(fid,'\n.lib ''%s'' TT\n',libname);
fprintf(fid,'.option converge=1\n');
%%temperature
if (temp)
    fprintf(fid,'.temp %d\n', temp_v);
end


for i_stage=1:n_bins
    n = 2 * optimizing_stage(i_stage);
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xkc(1,i_stage);
        Vtn = Xkc(2,i_stage);
        Vtp = Xkc(3,i_stage);
    end

    if aging
        [Vtn_age,Vtp_age] = aging_fnc(temp_v,Vd,Vtn,Vtp);
        Vtn = Vtn_age;
        Vtp = Vtp_age;
    end

    if (vt_var_mode ~= 0)
        Vtn = Vtn + vtn_var_val + vtn_percent * Vtn;
        Vtp = Vtp + vtp_var_val + vtp_percent * Vtp;
    end

    fprintf(fid,'.param sigma=0 l=%d clock_b%i_ini=%4.4fn cl=%d vd_b%i_ini=%d Vtn_b%i_ini=%d Vtp_b%i_ini=%d rl=%d fanout=%3.2f fc_global_lib=''(mc_global==0)'' Cout=%d ', ...
                                    L, i_stage, 30000, C0, i_stage, Vd, i_stage, Vtn, i_stage, Vtp, R, fanout-1, Cout);
    for width_ind = 1:n
        fprintf(fid, 'w_b%i_%i=%d ', i_stage, width_ind, ceil(Xkc(3+width_ind,i_stage)/width_eff_mult));
    end
    fprintf(fid, '\n');


    %%%% 注意这里是Xkc！！！！！！！！！
    if (flag_var)
        fprintf(fid,'.subckt inv1_s_b%i_ini 1 2 3 dvt=0 w_b%i=1u\n', i_stage);
        fprintf(fid,'Xn1 7 4 0 0 d_nfet W=w_b%i L=l\n', i_stage);
        fprintf(fid,'Xn2 2 8 7 0 d_nfet W=w_b%i L=l\n');
        fprintf(fid,'Xp1 2 9 3 3 d_pfet W=w_b%i L=l\n');
        fprintf(fid,'Xp2 5 3 3 d_pfet W=w_b%i L=l\n');
        fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
        fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
        fprintf(fid,'Vtno 8 3 ''dvt+vshift'' \n');
        fprintf(fid,'Vtpo 9 3 ''-dvt-vshift''\n.ends\n');            
        fprintf(fid,'.subckt inv2_s_b%i_ini 1 2 3 dvt=0 w_b%i=1u\n');
        fprintf(fid,'Xn 6 4 0 0 d_nfet W=w_b%i L=l\n');
        fprintf(fid,'Xp 6 5 3 3 d_pfet W=w_b%i L=l\n');
        fprintf(fid,'Vtn 4 1 ''dvt+vshift''\n');
        fprintf(fid,'Vtp 5 1 ''-dvt-vshift'' \n');
        fprintf(fid,'Rl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s_b%i_ini 100_s_b%i_ini 1_s_b%i_ini vdr_s_b%i_ini inv2_s_b%i_ini w=''w_b%i_1-%4.4f'' dvt=%3.3f\n', ...
                                    i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, dw1, dvt);
        for i=1:1:n
            if(mod(i,2))
                fprintf(fid,'xinv%i_s_b%i_ini %i_s_b%i_ini %i_s_b%i_ini vdd_s_b%i_ini inv1_s_b%i_ini w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                    i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
            else
                fprintf(fid,'xinv%i_s_b%i_ini %i_s_b%i_ini %i_s_b%i_ini vdd_s_b%i_ini inv2_s_b%i_ini w=''w_b%i_%i-%4.4f'' dvt=%3.3f\n', ...
                                    i, i_stage, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i, dw1, dvt);
                if(i<i_stage*2)
                    fprintf(fid,'xinv%i_s_b%i_ini %i_s_b%i_ini %i_s_b%i_ini vdr_s_b%i_ini inv1_s_b%i_ini w=''fanout*w_b%i_%i-fanout*%4.4f'' dvt=%3.3f\n', ...
                                    100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, dw1, dvt);
                else
                    fprintf(fid,'xinv%i_s_b%i_ini %i_s_b%i_ini %i_S_b%i_ini vdr_s_b%i_ini inv1_s_b%i_ini w=''fanout*w_b%i_%i+w_b%i_%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n', ...
                                    100+i-1, i_stage, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1, dw1, dw1, dvt);
                end
                
            end
        end
        fprintf(fid,'Vdr_s_b%i_ini vdr_s_b%i_ini 0 DC=''vd_b%i_ini*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
        fprintf(fid,'Vd_s_b%i_ini vdd_s_b%i_ini 0 DC=''vd_b%i_ini*%2.2f''\n', i_stage, i_stage, i_stage, p_slow);
        fprintf(fid,'Vpulse_s_b%i_ini 100_s_b%i_ini 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock_b%i_ini/2'',clock_b%i_ini)\n', i_stage, i_stage, p_slow, i_stage, i_stage);
    else
        % NAND Gate Netlist FOR INITIAL SIMULATION
        fprintf(fid,['.subckt inv1_b%i_ini 1 2 3 \n', ...
        'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_ini p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_ini p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_ini p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_ini p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        '.ends\n'], i_stage, i_stage, i_stage, i_stage, i_stage);
        % INV Gate Netlistzc FOR INITIAL SIMULATION
        fprintf(fid,['.subckt inv2_b%i_ini 1 2 3 \n', ...
        'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn_b%i_ini p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp_b%i_ini p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Rl 6 2 rl\n', ...
        'Cl 2 0 cl\n', ...
        '.ends\n'], i_stage, i_stage, i_stage);

        % write the netlist for the driver
        fprintf(fid, 'xinv100_b%i_ini 100_b%i_ini 1_b%i_ini vdr_b%i_ini inv2_b%i_ini w=w_b%i_1\n', i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);
        for i=1:1:n
            if(mod(i,2))
                fprintf(fid,'xinv_b%i_%i_ini %i_b%i_ini %i_b%i_ini vdd_b%i_ini inv1_b%i_ini w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
            else
                fprintf(fid,'xinv_b%i_%i_ini %i_b%i_ini %i_b%i_ini vdd_b%i_ini inv2_b%i_ini w=w_b%i_%i\n', i_stage, i, i, i_stage, i+1, i_stage, i_stage, i_stage, i_stage, i);
                if(i<i_stage*2)
                    fprintf(fid,'xinv_b%i_%i_ini %i_b%i_ini %i_b%i_ini vdr_b%i_ini inv1_b%i_ini w=''fanout*w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1);
                else
                    fprintf(fid,'xinv_b%i_%i_ini %i_b%i_ini %i_b%i_ini vdr_b%i_ini inv1_b%i_ini w=''fanout*w_b%i_%i+w_b%i_%i''\n', i_stage, 100+i-1, i+1, i_stage, 200+i, i_stage, i_stage, i_stage, i_stage, i-1, i_stage, i-1);
                end
            end
        end

        fprintf(fid,['Vdr_b%i_ini vdr_b%i_ini 0 DC vd_b%i_ini\n', ...
                    'Vd_b%i_ini vdd_b%i_ini 0 DC vd_b%i_ini\n', ...
                    'Vpulse_b%i_ini 100_b%i_ini 0 PULSE(0,vd_b%i_ini,10n,10p,10p,''clock_b%i_ini/2'',clock_b%i_ini)\n', ...
                    '.tran ''clock_b%i_ini/10000'' ''clock_b%i_ini*5+10n'' \n'], ...
                        i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage, i_stage);

        if (flag_var)
            fprintf(fid,'.measure tran t1_b%i trig V(1_s_b%i_ini) val=''%2.2f*vd_b%i_ini/2'' TD=''10n+clock_b%i_ini/4'' rise=2 targ V(%i_s_b%i_ini) val=''%2.2f*vd_b%i_ini/2'' TD=''10n+clock_b%i_ini/4'' %s=2\n', ...
                                            i_stage, i_stage, p_slow, i_stage, i_stage, n+1, i_stage, p_slow, i_stage, i_stage, edge);
        else
            % Only measure the delay here
            fprintf(fid,'.measure avgdelay_b%i param=''(t1_b%i+t2_b%i)*0.5''\n', i_stage, i_stage, i_stage);
            fprintf(fid,'.measure b%i_holder param=0\n', i_stage);
            fprintf(fid,'.measure tran t1_b%i trig V(1_b%i_ini) val=''vd_b%i_ini*0.5'' TD=0 rise=3 targ V(%i_b%i_ini) val=''vd_b%i_ini*0.5'' TD=0 rise=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
            fprintf(fid,'.measure tran t2_b%i trig V(1_b%i_ini) val=''vd_b%i_ini*0.5'' TD=0 fall=3 targ V(%i_b%i_ini) val=''vd_b%i_ini*0.5'' TD=0 fall=3\n', i_stage, i_stage, i_stage, i_stage*2+1, i_stage, i_stage);
        end

        fprintf(fid,'.ic ');

        if (flag_var)
            fprintf(fid,'v(100_s_b%i_ini)=0 ', i_stage);
            for i=1:n
                if(mod(i,2)==0)
                    fprintf(fid,'v(%i_s_b%i_ini)=0 ', i, i_stage);
                else
                    fprintf(fid,'v(%i_s_b%i_ini)=''vd_b%i_ini*%2.2f'' ', i, i_stage, i_stage, p_slow);
                end
            end
        else
            fprintf(fid,'v(100_b%i_ini)=0 ', i_stage);
            for i=1:n
                if(mod(i,2)==0)
                    fprintf(fid,'v(%i_b%i_ini)=0 ', i, i_stage);
                else
                    fprintf(fid,'v(%i_b%i_ini)=''vd_b%i_ini'' ', i, i_stage, i_stage);
                end
            end
        end
    end
    fprintf(fid,'\n\n\n');
end

fprintf(fid,'.end');
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
    if(found)
        hspice_running = 0;
        fprintf('Hspice simulation done!');
    end
end


% Read the simulation results
fid=fopen(resultname, 'r');
if fid == -1
    error('File cannot be opened');
end

for (j=1:(2+n_bins+1+1))
    tline=fgetl(fid);
end
for i=1:1:n_bins
    [avgdly r] = strtok(tline);
    if(ischar(avgdly))
        [Delay1(i), status] = str2num(avgdly);
    else
        status = 0;
    end
    if( status==0 || (Delay1(i)>1e-4) || (Delay1(i)<0))
        f2=1;
        error('Xk: Delay1 is not correct');
    end
    tline=fgetl(fid);
end
fclose(fid);
fprintf('Estimated Delay1 for each bin: \n');
disp(Delay1);
fprintf("VDD1: %d, VDD2: %d\n", Xkc(1, min(optimizing_stage)), Xkc(1, max(optimizing_stage)));
fprintf("Vtn1: %d, Vtn2: %d\n", Xkc(2, min(optimizing_stage)), Xkc(2, max(optimizing_stage)));
fprintf("Vtp1: %d, Vtp2: %d\n", Xkc(3, min(optimizing_stage)), Xkc(3, max(optimizing_stage)));

% delay calibration factor
Delay1 = Delay1./calib_delay;
