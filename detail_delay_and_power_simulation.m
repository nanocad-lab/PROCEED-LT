
fprintf('Detail Delay and Power Sim NEW - ');
fprintf('Delay1 Value: %d\n', Delay1);
% choose model for this stage
% error_detect = 1;
% for r_allocation = 1: length(devicemodel_allocation(:,1))
%     if( length( find (devicemodel_allocation(r_allocation,:) == i_stage ) ) == 1 )
%         subcktname = devicemodel(r_allocation).subcktname;
%         error_detect = 0;
%         break;
%     end
% end
% if( error_detect)
%     print('no model found for this bin stage');
%     %break;
% end

libname = ['YOUR_LIBRARY_NAME'];
%%end

% LEAKAGE POWER CALCULATION
leak_arr = zeros(2,1);
for leak_iter = 1:1:4
    fid=fopen(spname,'w');
    % fprintf(fid,'\n.inc ''%s'' \n',subcktname);
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=1\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end

    if aging
        Vdd_n = (Vd - Vtn);
        Vdd_p = (Vd - Vtp);
        % Vdd_n check
        if Vdd_n<0.305
            Vtn = Vtn;
        elseif(Vdd_n>=0.95)
            Vtn = Vtn + deltav(3,650);
            fprintf('deltav_N: %d -',deltav(3,650));
        else
            index = int16((Vdd_n-0.3)*1000); 
            Vtn = Vtn + deltav(3,index);
            fprintf('deltav_N: %d -',deltav(3,index));
        end
        % Vdd_p check
        if Vdd_p<0.305
            Vtp = Vtp;
        elseif(Vdd_p>=0.95)
            Vtp = Vtp + deltav(4,650);
            fprintf('deltav_P: %d\n',deltav(4,650));
        else
            index = int16((Vdd_p-0.3)*1000); 
            Vtp = Vtp + deltav(4,index);
            fprintf('deltav_P: %d\n',deltav(4,index));
        end
    else
        Vtn = Vtn;
        Vtp = Vtp;
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);

    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    end
    fprintf(fid,'\n');

    % NAND Netlist FOR LEAK TEST
    fprintf(fid,['.subckt inv1 1 2 3 4\n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 4 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 4 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    %INV Netlist FOR LEAK TEST
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
    for i=1:1:n
        if(mod(i,2))
            % fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
            fprintf(fid,'xinv%i in1 %i vdd in2 inv1 w=w%i\n',i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,100+i-1,i);
            % if(i<n)
            %     fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
            % else
            %     fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
            % end
        %  fprintf(fid,'c%i %i 0 0.1fF\n',100+i-1,200+i);
        end
    end

    % OLD METHOD --> LOOK THROUGH LOG FILE
    if (leak_iter == 1) %00
        fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC 0\n Vin2 in2 0 DC 0\n');
    elseif (leak_iter == 2) %01
        % fprintf(fid,'Vd vdd 0 DC vd\n');
        fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC 0\n Vin2 in2 0 DC vd\n');
    elseif (leak_iter == 3) %10
        fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC vd\n Vin2 in2 0 DC 0\n');
    elseif (leak_iter == 4) %11
        fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC vd\n Vin2 in2 0 DC vd\n');
    end

    fprintf(fid,'.DC sweep clock start=1n stop=2n step=1n\n');
    fprintf(fid,'.option ingold=1 brief=1 \n');
    fprintf(fid,'.print par''V(vdd)*I(Vd)''\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    leakage2temp=38;
    fid=fopen(logname,'r');
    f3=1;
    while(~feof(fid))
        tline=fgetl(fid);
        [t r]=strtok(tline);
        if(strcmp(t,'v(vdd)*i(vd)'))
            tline=fgetl(fid);
            [t r]=strtok(tline);
            [t r]=strtok(r);
            [x status]=str2num(t);
            if(status~=0)
                leakage2temp=x;
                f3=0;
            end
            break;
        end
    end
    fclose(fid);

    if (f3 == 0)
        leak_arr(leak_iter,1) = (-1)*abs(leakage2temp);
    end
    % fprintf('Leak Power %i: %d\n', leak_iter, leak_arr(leak_iter,1));
end
Leakage2 = mean(leak_arr);
% fprintf('Leak Power Method 1 Stored: %d\n', Leakage2);

% ---------------------------
% delay and power calculation
%%if aging, calculate power and latency seperately
if aging
    fid=fopen(spname,'w');
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=1\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end

    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        end
    end
    fprintf(fid,'\n');

    % fprintf(fid,'.subckt inv1 1 2 3 \nXn1 7 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXn2 2 3 7 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp1 2 3 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nXp2 2 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\n.ends\n');
    % fprintf(fid,'.subckt inv2 1 2 3 \nXn 6 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp 6 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nRl 6 2 rl\nCl 2 0 cl\n.ends\n');
    % NAND Netlist
    fprintf(fid,['.subckt inv1 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    % INVERTER Netlist
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
    for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,i+1,i);
            if(i<n)
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
                % for (j = 1:1:ceil(fanout-1))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            else
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
                % for (j = 1:1:ceil(fanout))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            end
        end
    end

    if(flag_var)
        fprintf(fid,'.subckt inv1_s 1 2 3 dvt=0 w=1u\nXn1 7 4 0 0 d_nfet W=w L=l\nXn2 2 8 7 0 d_nfet W=w L=l\nXp1 2 9 3 3 d_pfet W=w L=l\nXp2 2 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nVtno 8 3 ''dvt+vshift'' \nVtpo 9 3 ''-dvt-vshift''\n.ends\n');
        fprintf(fid,'.subckt inv2_s 1 2 3 dvt=0 w=1u\nXn 6 4 0 0 d_nfet W=w L=l\nXp 6 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nRl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s 100_s 1_s vdr_s inv2_s w=''w1-%4.4f'' dvt=%3.3f\n',dw1,dvt);
        for i=1:1:n
            if(mod(i,2))
                fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv1_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
            else
                fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv2_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
                if(i<n)
                    fprintf(fid,'xinv%i_s %i_s %i_s vdr_s inv1_s w=''fanout*w%i-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,dw1,dvt);
                else
                    fprintf(fid,'xinv%i_s %i_s %i_S vdr_s inv1_s w=''fanout*w%i+w%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,i-1,dw1,dw1,dvt);
                end
                
            end
        end
        fprintf(fid,'Vdr_s vdr_s 0 DC=''vd*%2.2f''\nVd_s vdd_s 0 DC=''vd*%2.2f''\nVpulse_s 100_s 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock/2'',clock)\n',p_slow,p_slow,p_slow);
    end

    % fprintf(fid,'Vdr vdr 0 DC vd\nVd vdd 0 DC vd\nVpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n.tran ''clock/1000'' ''clock*3+10n'' \n');
    fprintf(fid,['Vdr vdr 0 DC vd\n', ...
    'Vd vdd 0 DC vd\n', ...
    'Vpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n', ...
    '.tran ''clock/1000'' ''clock*5+10n'' \n']);
    
    if(flag_var)
        fprintf(fid,'.measure tran t1 trig V(1_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' %s=2\n',p_slow,n+1,p_slow,edge);
    else
        fprintf(fid,'.measure tran t1 trig V(1) val=''vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i) val=''vd/2'' TD=''10n+clock/4'' %s=2\n',n+1,edge);
        fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))'' FROM=''clock*2+10n'' TO=''clock*5+10n''\n');
        fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.2+1)'' TO=''10n+clock*(1.4+1)''\n');
        fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.7+1)'' TO=''10n+clock*(1.9+1)''\n');
        
        % fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        % fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10n'' TO=''clock*3+10n''\n');
        % fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        % fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        % fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        % fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.07n+3*clock+500n'' TO=''10.07n+1*10000n+3*clock-500n''\n');
        % fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.08n+1*10000n+3*clock+500n'' TO=''10.08n+2*10000n+3*clock-500n''\n');
    end

    % fprintf(fid,'\n.measure tran avgpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*0.3+10n'' TO=''clock*2.3+10n''\n');
    % fprintf(fid,'.measure tran leakpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*1.3+10n+10p'' TO=''clock*1.4+10n+10p''\n');
    fprintf(fid,'.ic ');
    fprintf(fid,'v(100)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i)=0 ',i);
        else
            fprintf(fid,'v(%i)=''vd'' ',i);
        end
    end
    if(flag_var)
        fprintf(fid,'v(100_s)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i_s)=0 ',i);
        else
            fprintf(fid,'v(%i_s)=''vd*%2.2f'' ',i,p_slow);
        end
    end
    end
    fprintf(fid,'\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    fid=fopen(resultname, 'r');
    % delay3 is just a placeholder so the code won't break
    if (fid ~= -1)
        x = [];
        for (j=1:(5+1))
            tline=fgetl(fid);
        end
        [t r] = strtok(tline);
        if(~ischar(t))
            f4=1;
            Delay3=38;
        else
            [Delay3,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Delay3=38;
        end
        if(Delay3<=0)
            f4=1;
        end
        
        [t r] = strtok(r);
        if(~ischar(t))
            f4=1;
        else
            [Energy2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Energy2=38;
        end
                
        if(Energy2>=0  || Energy2>Leakage2/2)
            f4=1;
            disp('energy2 >=0 or energy2 < leakage');
            disp(Energy2);
            disp(Leakage2);
        else
            if(Energy2>Leakage2)
                Energy2=Leakage2;
            end
        end
    else
        f4=1;
    end
    fclose(fid);

    %second round for latency
    fid=fopen(spname,'w');
    % fprintf(fid,'\n.inc ''%s'' \n',subcktname);
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=1\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end
    if aging
        Vdd_n = (Vd - Vtn);
        Vdd_p = (Vd - Vtp);
        % Vdd_n check
        if Vdd_n<0.305
            Vtn = Vtn;
        elseif(Vdd_n>=0.95)
            Vtn = Vtn + deltav(3,650);
            fprintf('deltav_N: %d -',deltav(3,650));
        else
            index = int16((Vdd_n-0.3)*1000); 
            Vtn = Vtn + deltav(3,index);
            fprintf('deltav_N: %d -',deltav(3,index));
        end
        % Vdd_p check
        if Vdd_p<0.305
            Vtp = Vtp;
        elseif(Vdd_p>=0.95)
            Vtp = Vtp + deltav(4,650);
            fprintf('deltav_P: %d\n',deltav(4,650));
        else
            index = int16((Vdd_p-0.3)*1000); 
            Vtp = Vtp + deltav(4,index);
            fprintf('deltav_P: %d\n',deltav(4,index));
        end
    else
        Vtn = Vtn;
        Vtp = Vtp;
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        end
    end
    fprintf(fid,'\n');

    % fprintf(fid,'.subckt inv1 1 2 3 \nXn1 7 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXn2 2 3 7 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp1 2 3 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nXp2 2 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\n.ends\n');
    % fprintf(fid,'.subckt inv2 1 2 3 \nXn 6 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp 6 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nRl 6 2 rl\nCl 2 0 cl\n.ends\n');
    % NAND Netlist
    fprintf(fid,['.subckt inv1 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    % INVERTER Netlist
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
    for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,i+1,i);
            if(i<n)
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
                % for (j = 1:1:ceil(fanout-1))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            else
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
                % for (j = 1:1:ceil(fanout))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            end
        end
    end

    if(flag_var)
        fprintf(fid,'.subckt inv1_s 1 2 3 dvt=0 w=1u\nXn1 7 4 0 0 d_nfet W=w L=l\nXn2 2 8 7 0 d_nfet W=w L=l\nXp1 2 9 3 3 d_pfet W=w L=l\nXp2 2 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nVtno 8 3 ''dvt+vshift'' \nVtpo 9 3 ''-dvt-vshift''\n.ends\n');
        fprintf(fid,'.subckt inv2_s 1 2 3 dvt=0 w=1u\nXn 6 4 0 0 d_nfet W=w L=l\nXp 6 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nRl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s 100_s 1_s vdr_s inv2_s w=''w1-%4.4f'' dvt=%3.3f\n',dw1,dvt);
        for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv1_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
        else
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv2_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
            if(i<n)
                fprintf(fid,'xinv%i_s %i_s %i_s vdr_s inv1_s w=''fanout*w%i-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,dw1,dvt);
            else
                fprintf(fid,'xinv%i_s %i_s %i_S vdr_s inv1_s w=''fanout*w%i+w%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,i-1,dw1,dw1,dvt);
            end
            
        end
        end
        fprintf(fid,'Vdr_s vdr_s 0 DC=''vd*%2.2f''\nVd_s vdd_s 0 DC=''vd*%2.2f''\nVpulse_s 100_s 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock/2'',clock)\n',p_slow,p_slow,p_slow);
    end

    % fprintf(fid,'Vdr vdr 0 DC vd\nVd vdd 0 DC vd\nVpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n.tran ''clock/1000'' ''clock*3+10n'' \n');
    fprintf(fid,['Vdr vdr 0 DC vd\n', ...
    'Vd vdd 0 DC vd\n', ...
    'Vpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n', ...
    '.tran ''clock/1000'' ''clock*5+10n'' \n']);

    if(flag_var)
        fprintf(fid,'.measure tran t1 trig V(1_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' %s=2\n',p_slow,n+1,p_slow,edge);
    else
        fprintf(fid,'.measure tran t1 trig V(1) val=''vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i) val=''vd/2'' TD=''10n+clock/4'' %s=2\n',n+1,edge);
        fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))'' FROM=''clock*2+10n'' TO=''clock*5+10n''\n');
        fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.2+1)'' TO=''10n+clock*(1.4+1)''\n');
        fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.7+1)'' TO=''10n+clock*(1.9+1)''\n');

        % fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        % fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10n'' TO=''clock*3+10n''\n');
        % fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        % fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        % fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        % fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.07n+3*clock+500n'' TO=''10.07n+1*10000n+3*clock-500n''\n');
        % fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.08n+1*10000n+3*clock+500n'' TO=''10.08n+2*10000n+3*clock-500n''\n');
    end

    % fprintf(fid,'\n.measure tran avgpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*0.3+10n'' TO=''clock*2.3+10n''\n');
    % fprintf(fid,'.measure tran leakpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*1.3+10n+10p'' TO=''clock*1.4+10n+10p''\n');
    fprintf(fid,'.ic ');
    fprintf(fid,'v(100)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i)=0 ',i);
        else
            fprintf(fid,'v(%i)=''vd'' ',i);
        end
    end
    if(flag_var)
        fprintf(fid,'v(100_s)=0 ');
        for i=1:n
            if(mod(i,2)==0)
                fprintf(fid,'v(%i_s)=0 ',i);
            else
                fprintf(fid,'v(%i_s)=''vd*%2.2f'' ',i,p_slow);
            end
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    % energy3 is just a placeholder so the code won't break
    fid=fopen(resultname, 'r');
    if (fid ~= -1)
        x = [];
        for (j=1:(5+1))
            tline=fgetl(fid);
        end
        [t r] = strtok(tline);
        if(~ischar(t))
            f4=1;
            Delay2=38;
        else
            [Delay2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Delay2=38;
        end
        if(Delay2<=0)
            f4=1;
        end
        [t r] = strtok(r);
        if(~ischar(t))
            f4=1;
        else
            [Energy3,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Energy3=38;
        end
                
    else
        f4=1;
    end
    fclose(fid);

%no aging
else
    fid=fopen(spname,'w');
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=1\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    end
    fprintf(fid,'\n');

    % NAND Netlist
    fprintf(fid,['.subckt inv1 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    % INVERTER Netlist
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,i+1,i);
            if(i<n)
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
            else
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
            end
        end
    end

    if(flag_var)
        fprintf(fid,'.subckt inv1_s 1 2 3 dvt=0 w=1u\nXn1 7 4 0 0 d_nfet W=w L=l\nXn2 2 8 7 0 d_nfet W=w L=l\nXp1 2 9 3 3 d_pfet W=w L=l\nXp2 2 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nVtno 8 3 ''dvt+vshift'' \nVtpo 9 3 ''-dvt-vshift''\n.ends\n');
        fprintf(fid,'.subckt inv2_s 1 2 3 dvt=0 w=1u\nXn 6 4 0 0 d_nfet W=w L=l\nXp 6 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nRl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s 100_s 1_s vdr_s inv2_s w=''w1-%4.4f'' dvt=%3.3f\n',dw1,dvt);
        for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv1_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
        else
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv2_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
            if(i<n)
                fprintf(fid,'xinv%i_s %i_s %i_s vdr_s inv1_s w=''fanout*w%i-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,dw1,dvt);
            else
                fprintf(fid,'xinv%i_s %i_s %i_S vdr_s inv1_s w=''fanout*w%i+w%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,i-1,dw1,dw1,dvt);
            end
            
        end
        end
        fprintf(fid,'Vdr_s vdr_s 0 DC=''vd*%2.2f''\nVd_s vdd_s 0 DC=''vd*%2.2f''\nVpulse_s 100_s 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock/2'',clock)\n',p_slow,p_slow,p_slow);
    end

    % fprintf(fid,'Vdr vdr 0 DC vd\nVd vdd 0 DC vd\nVpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n.tran ''clock/1000'' ''clock*3+10n'' \n');
    fprintf(fid,['Vdr vdr 0 DC vd\n', ...
    'Vd vdd 0 DC vd\n', ...
    'Vpulse 100 0 PULSE(0,vd,10n,10p,10p,''clock/2'',clock)\n', ...
    '.tran ''clock/10000'' ''clock*5+10n'' \n']);

    if(flag_var)
        fprintf(fid,'.measure tran t1 trig V(1_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' %s=2\n',p_slow,n+1,p_slow,edge);
    else
        fprintf(fid,'.measure tran t1 trig V(1) val=''vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i) val=''vd/2'' TD=''10n+clock/4'' %s=2\n',n+1,edge);
        fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))'' FROM=''clock*2+10n'' TO=''clock*5+10n''\n');
        fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.2+1)'' TO=''10n+clock*(1.4+1)''\n');
        fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.7+1)'' TO=''10n+clock*(1.9+1)''\n');

        % fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        % fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10n'' TO=''clock*3+10n''\n');
        % fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        % fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        % fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        % fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.07n+3*clock+500n'' TO=''10.07n+1*10000n+3*clock-500n''\n');
        % fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.08n+1*10000n+3*clock+500n'' TO=''10.08n+2*10000n+3*clock-500n''\n');
    end

    % fprintf(fid,'\n.measure tran avgpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*0.3+10n'' TO=''clock*2.3+10n''\n');
    % fprintf(fid,'.measure tran leakpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*1.3+10n+10p'' TO=''clock*1.4+10n+10p''\n');
    fprintf(fid,'.ic ');
    fprintf(fid,'v(100)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i)=0 ',i);
        else
            fprintf(fid,'v(%i)=''vd'' ',i);
        end
    end
    if(flag_var)
        fprintf(fid,'v(100_s)=0 ');
        for i=1:n
            if(mod(i,2)==0)
                fprintf(fid,'v(%i_s)=0 ',i);
            else
                fprintf(fid,'v(%i_s)=''vd*%2.2f'' ',i,p_slow);
            end
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    fid=fopen(resultname, 'r');
    if (fid ~= -1)
        x = [];
        for (j=1:(5+1))
            tline=fgetl(fid);
        end
        [t r] = strtok(tline);
        if(~ischar(t))
            f4=1;
            Delay2=38;
        else
            [Delay2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Delay2=38;
        end
        if(Delay2<=0)
            f4=1;
        end
        
        [t r] = strtok(r);
        if(~ischar(t))
            f4=1;
        else
            [Energy2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Energy2=38;
        end
                
        if(Energy2>=0  || Energy2>Leakage2/2)
            f4=1;
            disp('energy2 >=0 or energy2 < leakage');
            disp(Energy2);
            disp(Leakage2);
        else
            if(Energy2>Leakage2)
                Energy2=Leakage2;
            end
        end
    else
        f4=1;
    end
    fclose(fid);
end


% if f4
if(f4==1)
    f4=0;
    leak_arr = zeros(2,1);

    for (leak_iter = 1:1:4)
        fid=fopen(spname,'w');
        % fprintf(fid,'\n.inc ''%s'' \n',subcktname);
        fprintf(fid,'\n.lib ''%s'' TT\n',libname);
        fprintf(fid,'.option converge=3\n');
        %%temperature
        if (temp)
            fprintf(fid,'.temp %d\n', temp_v);
        end
        if fixed_V
            Vd = 0.7;
            Vtn = 0;
            Vtp = 0;
        else
            Vd = Xk(1);
            Vtn = Xk(2);
            Vtp = Xk(3);
        end

        if aging
            Vdd_n = (Vd - Vtn);
            Vdd_p = (Vd - Vtp);
            % Vdd_n check
            if Vdd_n<0.305
                Vtn = Vtn;
            elseif(Vdd_n>=0.95)
                Vtn = Vtn + deltav(3,650);
                fprintf('deltav_N: %d -',deltav(3,650));
            else
                index = int16((Vdd_n-0.3)*1000); 
                Vtn = Vtn + deltav(3,index);
                fprintf('deltav_N: %d -',deltav(3,index));
            end
            % Vdd_p check
            if Vdd_p<0.305
                Vtp = Vtp;
            elseif(Vdd_p>=0.95)
                Vtp = Vtp + deltav(4,650);
                fprintf('deltav_P: %d\n',deltav(4,650));
            else
                index = int16((Vdd_p-0.3)*1000); 
                Vtp = Vtp + deltav(4,index);
                fprintf('deltav_P: %d\n',deltav(4,index));
            end
        else
            Vtn = Vtn;
            Vtp = Vtp;
        end

        fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1,Cout);
        % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1,Cout);
        % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
        if flag_w
            for i=1:n
                fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
            %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
            end
        else
            for i=1:n
                fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
            %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
            end
        end
        fprintf(fid,'\n');

        % fprintf(fid,'.subckt inv1 1 2 3 \nXn1 7 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXn2 2 3 7 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp1 2 3 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nXp2 2 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\n.ends\n');
        %     fprintf(fid,'.subckt inv2 1 2 3 \nXn 6 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp 6 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nRl 6 2 rl\nCl 2 0 cl\n.ends\n');
        % NAND Netlist FOR LEAK TEST
        fprintf(fid,['.subckt inv1 1 2 3 4\n', ...
        'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xn2 2 4 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp1 2 4 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        '.ends\n']);
        % INV Netlist FOR LEAK TEST
        fprintf(fid,['.subckt inv2 1 2 3 \n', ...
        'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
        'Rl 6 2 rl\n', ...
        'Cl 2 0 cl\n', ...
        '.ends\n']);

        % fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
        % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
        for i=1:1:n
            if(mod(i,2))
                % fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
                fprintf(fid,'xinv%i in1 %i vdd in2 inv1 w=w%i\n',i,i+1,i);
            else
                fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,100+i-1,i);
                % if(i<n)
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
                % else
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
                % end
            %  fprintf(fid,'c%i %i 0 0.1fF\n',100+i-1,200+i);
            end
        end

        if (leak_iter == 1) %00
            fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC 0\n Vin2 in2 0 DC 0\n');
        elseif (leak_iter == 2) %01
            % fprintf(fid,'Vd vdd 0 DC vd\n');
            fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC 0\n Vin2 in2 0 DC vd\n');
        elseif (leak_iter == 3) %10
            fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC vd\n Vin2 in2 0 DC 0\n');
        elseif (leak_iter == 4) %11
            fprintf(fid,'Vd vdd 0 DC vd\n Vin1 in1 0 DC vd\n Vin2 in2 0 DC vd\n');
        end

        fprintf(fid,'.DC sweep clock start=1n stop=2n step=1n\n');
        fprintf(fid,'.option ingold=1 brief=1 \n');
        fprintf(fid,'.print par''V(vdd)*I(Vd)'' \n');
        fprintf(fid,'.end');
        fclose(fid);
        unix(hspicerun);

        leakage2temp=38;
        fid=fopen(logname,'r');
        f3=1;
        while(~feof(fid))
            tline=fgetl(fid);
            [t r]=strtok(tline);
            if(strcmp(t,'v(vdd)*i(vd)'))
                tline=fgetl(fid);
                    [t r]=strtok(tline);
                    [t r]=strtok(r);
                    [x status]=str2num(t);
                    if(status==0)
                        x=-1;
                    end
                    leakage2temp=x;
                    break;
            end
        end
        fclose(fid);
        if (f3 == 0)
            leak_arr(leak_iter,1) = (-1)*abs(leakage2temp);
        end
    end
    Leakage2 = mean(leak_arr);
    % fprintf('Leak Power Method 2 Stored: %d\n', Leakage2);

    if aging
    fid=fopen(spname,'w');
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=3\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    end
    fprintf(fid,'\n');

    % fprintf(fid,'.subckt inv1 1 2 3 \nXn1 7 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXn2 2 3 7 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp1 2 3 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nXp2 2 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\n.ends\n');
    %     fprintf(fid,'.subckt inv2 1 2 3 \nXn 6 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp 6 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nRl 6 2 rl\nCl 2 0 cl\n.ends\n');
    % NAND Netlist
    fprintf(fid,['.subckt inv1 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    % INVERTER Netlist
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
    for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,i+1,i);
            if(i<n)
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
                % for (j = 1:1:ceil(fanout-1))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            else
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
                % for (j = 1:1:ceil(fanout))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            end
        end
    end

    if(flag_var)
        fprintf(fid,'.subckt inv1_s 1 2 3 dvt=0 w=1u\nXn1 7 4 0 0 d_nfet W=w L=l\nXn2 2 8 7 0 d_nfet W=w L=l\nXp1 2 9 3 3 d_pfet W=w L=l\nXp2 2 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nVtno 8 3 ''dvt+vshift'' \nVtpo 9 3 ''-dvt-vshift''\n.ends\n');
        fprintf(fid,'.subckt inv2_s 1 2 3 dvt=0 w=1u\nXn 6 4 0 0 d_nfet W=w L=l\nXp 6 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nRl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s 100_s 1_s vdr_s inv2_s w=''w1-%4.4f'' dvt=%3.3f\n',dw1,dvt);
        for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv1_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
        else
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv2_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
            if(i<n)
                fprintf(fid,'xinv%i_s %i_s %i_s vdr_s inv1_s w=''fanout*w%i-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,dw1,dvt);
            else
                fprintf(fid,'xinv%i_s %i_s %i_S vdr_s inv1_s w=''fanout*w%i+w%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,i-1,dw1,dw1,dvt);
            end
            
        end
        end
        fprintf(fid,'Vdr_s vdr_s 0 DC=''vd*%2.2f''\nVd_s vdd_s 0 DC=''vd*%2.2f''\nVpulse_s 100_s 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock/2'',clock)\n',p_slow,p_slow,p_slow);
    end 

    % fprintf(fid,'Vdr vdr 0 DC vd\nVd vdd 0 DC vd\nVpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n.tran ''clock/1000'' ''clock*3+10n'' \n');
    fprintf(fid,['Vdr vdr 0 DC vd\n', ...
    'Vd vdd 0 DC vd\n', ...
    'Vpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n', ...
    '.tran ''clock/1000'' ''clock*5+10n'' \n']);

    if(flag_var)
        fprintf(fid,'.measure tran t1 trig V(1_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' %s=2\n',p_slow,n+1,p_slow,edge);
    else
        fprintf(fid,'.measure tran t1 trig V(1) val=''vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i) val=''vd/2'' TD=''10n+clock/4'' %s=2\n',n+1,edge);
        fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))'' FROM=''clock*2+10n'' TO=''clock*5+10n''\n');
        fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.2+1)'' TO=''10n+clock*(1.4+1)''\n');
        fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.7+1)'' TO=''10n+clock*(1.9+1)''\n');

        % fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        % fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10n'' TO=''clock*3+10n''\n');
        % fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        % fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        % fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        % fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.07n+3*clock+500n'' TO=''10.07n+1*10000n+3*clock-500n''\n');
        % fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.08n+1*10000n+3*clock+500n'' TO=''10.08n+2*10000n+3*clock-500n''\n');
    end
    % fprintf(fid,'\n.measure tran avgpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*0.3+10n'' TO=''clock*2.3+10n''\n');
    % fprintf(fid,'.measure tran leakpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*1.3+10n+10p'' TO=''clock*1.4+10n+10p''\n');
    % fprintf(fid,'.measure VTH\n');
    fprintf(fid,'.ic ');
    fprintf(fid,'v(100)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i)=0 ',i);
        else
            fprintf(fid,'v(%i)=''vd'' ',i);
        end
    end
    if(flag_var)
        fprintf(fid,'v(100_s)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i_s)=0 ',i);
        else
            fprintf(fid,'v(%i_s)=''vd*%2.2f'' ',i,p_slow);
        end
    end
    end
    fprintf(fid,'\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    fid=fopen(resultname, 'r');
    if (fid ~= -1)
        x = [];
        for (j=1:(5+1))
            tline=fgetl(fid);
        end
        [t r] = strtok(tline);
        if(~ischar(t))
            f4=1;
            Delay3=38;
        else
            [Delay3,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Delay3=38;
        end
        if(Delay3<=0)
            f4=1;
        end
        
        [t r] = strtok(r);
        if(~ischar(t))
            f4=1;
        else
            [Energy2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Energy2=38;
        end
                
        if(Energy2>=0  || Energy2>Leakage2/2)
            f4=1;
            disp('energy2 >=0 or energy2 < leakage');
            disp(Energy2);
            disp(Leakage2);
        else
            if(Energy2>Leakage2)
                Energy2=Leakage2;
            end
        end
    else
        f4=1;
    end
    fclose(fid);

    %second round for latency
    fid=fopen(spname,'w');
    %fprintf(fid,'\n.inc ''%s'' \n',subcktname);
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=3\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end

    if aging
        Vdd_n = (Vd - Vtn);
        Vdd_p = (Vd - Vtp);
        % Vdd_n check
        if Vdd_n<0.305
            Vtn = Vtn;
        elseif(Vdd_n>=0.95)
            Vtn = Vtn + deltav(3,650);
        else
            index = int16((Vdd_n-0.3)*1000); 
            Vtn = Vtn + deltav(3,index);
        end
        % Vdd_p check
        if Vdd_p<0.305
            Vtp = Vtp;
        elseif(Vdd_p>=0.95)
            Vtp = Vtp + deltav(4,650);
        else
            index = int16((Vdd_p-0.3)*1000); 
            Vtp = Vtp + deltav(4,index);
        end
    else
        Vtn = Vtn;
        Vtp = Vtp;
    end

    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    end
    fprintf(fid,'\n');

    % fprintf(fid,'.subckt inv1 1 2 3 \nXn1 7 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXn2 2 3 7 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp1 2 3 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nXp2 2 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\n.ends\n');
    %     fprintf(fid,'.subckt inv2 1 2 3 \nXn 6 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp 6 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nRl 6 2 rl\nCl 2 0 cl\n.ends\n');
    % NAND Netlist
    fprintf(fid,['.subckt inv1 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    % INVERTER Netlist
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
    for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,i+1,i);
            if(i<n)
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
                % for (j = 1:1:ceil(fanout-1))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            else
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
                % for (j = 1:1:ceil(fanout))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            end
        end
    end

    if(flag_var)
        fprintf(fid,'.subckt inv1_s 1 2 3 dvt=0 w=1u\nXn1 7 4 0 0 d_nfet W=w L=l\nXn2 2 8 7 0 d_nfet W=w L=l\nXp1 2 9 3 3 d_pfet W=w L=l\nXp2 2 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nVtno 8 3 ''dvt+vshift'' \nVtpo 9 3 ''-dvt-vshift''\n.ends\n');
        fprintf(fid,'.subckt inv2_s 1 2 3 dvt=0 w=1u\nXn 6 4 0 0 d_nfet W=w L=l\nXp 6 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nRl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s 100_s 1_s vdr_s inv2_s w=''w1-%4.4f'' dvt=%3.3f\n',dw1,dvt);
        for i=1:1:n
            if(mod(i,2))
                fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv1_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
            else
                fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv2_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
                if(i<n)
                    fprintf(fid,'xinv%i_s %i_s %i_s vdr_s inv1_s w=''fanout*w%i-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,dw1,dvt);
                else
                    fprintf(fid,'xinv%i_s %i_s %i_S vdr_s inv1_s w=''fanout*w%i+w%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,i-1,dw1,dw1,dvt);
                end
                
            end
        end
        fprintf(fid,'Vdr_s vdr_s 0 DC=''vd*%2.2f''\nVd_s vdd_s 0 DC=''vd*%2.2f''\nVpulse_s 100_s 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock/2'',clock)\n',p_slow,p_slow,p_slow);
    end 

    % fprintf(fid,'Vdr vdr 0 DC vd\nVd vdd 0 DC vd\nVpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n.tran ''clock/1000'' ''clock*3+10n'' \n');
    fprintf(fid,['Vdr vdr 0 DC vd\n', ...
    'Vd vdd 0 DC vd\n', ...
    'Vpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n', ...
    '.tran ''clock/1000'' ''clock*5+10n'' \n']);

    if(flag_var)
        fprintf(fid,'.measure tran t1 trig V(1_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' %s=2\n',p_slow,n+1,p_slow,edge);
    else
        fprintf(fid,'.measure tran t1 trig V(1) val=''vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i) val=''vd/2'' TD=''10n+clock/4'' %s=2\n',n+1,edge);
        fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))'' FROM=''clock*2+10n'' TO=''clock*5+10n''\n');
        fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.2+1)'' TO=''10n+clock*(1.4+1)''\n');
        fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.7+1)'' TO=''10n+clock*(1.9+1)''\n');

        % fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        % fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10n'' TO=''clock*3+10n''\n');
        % fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        % fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        % fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        % fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.07n+3*clock+500n'' TO=''10.07n+1*10000n+3*clock-500n''\n');
        % fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.08n+1*10000n+3*clock+500n'' TO=''10.08n+2*10000n+3*clock-500n''\n');
    end
    % fprintf(fid,'\n.measure tran avgpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*0.3+10n'' TO=''clock*2.3+10n''\n');
    % fprintf(fid,'.measure tran leakpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*1.3+10n+10p'' TO=''clock*1.4+10n+10p''\n');
    % fprintf(fid,'.measure VTH\n');
    fprintf(fid,'.ic ');
    fprintf(fid,'v(100)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i)=0 ',i);
        else
            fprintf(fid,'v(%i)=''vd'' ',i);
        end
    end
    if(flag_var)
        fprintf(fid,'v(100_s)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i_s)=0 ',i);
        else
            fprintf(fid,'v(%i_s)=''vd*%2.2f'' ',i,p_slow);
        end
    end
    end
    fprintf(fid,'\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    fid=fopen(resultname, 'r');
    if (fid ~= -1)
        x = [];
        for (j=1:(5+1))
            tline=fgetl(fid);
        end
        [t r] = strtok(tline);
        if(~ischar(t))
            f4=1;
            Delay2=38;
        else
            [Delay2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Delay2=38;
        end
        if(Delay2<=0)
            f4=1;
        end
        [t r] = strtok(r);
        if(~ischar(t))
            f4=1;
        else
            [Energy3,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Energy3=38;
        end
    else
        f4=1;
    end
    fclose(fid);

    else

    fid=fopen(spname,'w');
    %fprintf(fid,'\n.inc ''%s'' \n',subcktname);
    fprintf(fid,'\n.lib ''%s'' TT\n',libname);
    fprintf(fid,'.option converge=3\n');
    %%temperature
    if (temp)
        fprintf(fid,'.temp %d\n', temp_v);
    end
    if fixed_V
        Vd = 0.7;
        Vtn = 0;
        Vtp = 0;
    else
        Vd = Xk(1);
        Vtn = Xk(2);
        Vtp = Xk(3);
    end

    fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*5,Delay1*2e9*5,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d fanout=%3.2f Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,fanout-1, Cout);
    % fprintf(fid,'.param sigma=0 l=%d pulse=%4.4fn clock=%4.4fn cl=%d vd=%d Vtn=%d Vtp=%d rl=%d Cout=%d ',L,Delay1*1e9*3,Delay1*2e9*3,C0,Vd,Vtn,Vtp,R,Cout);
    if flag_w
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(Xk(3+i)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    else
        for i=1:n
            fprintf(fid,'w%i=%d ',i,ceil(X(3+i,i_stage)/width_eff_mult));
        %     fprintf(fid,'w%i=%d ',i,X(2+i,i_stage)*1e-6);
        end
    end
    fprintf(fid,'\n');

    % fprintf(fid,'.subckt inv1 1 2 3 \nXn1 7 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXn2 2 3 7 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp1 2 3 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nXp2 2 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\n.ends\n');
    %     fprintf(fid,'.subckt inv2 1 2 3 \nXn 6 1 0 0 lvtnfet L=l nfin=w p_vta=Vtn\nXp 6 1 3 3 lvtpfet L=l nfin=w p_vta=Vtp\nRl 6 2 rl\nCl 2 0 cl\n.ends\n');
    fprintf(fid,['.subckt inv1 1 2 3 \n', ...
    'Xn1 7 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xn2 2 3 7 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp1 2 3 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp2 2 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    '.ends\n']);
    % INVERTER Netlist
    fprintf(fid,['.subckt inv2 1 2 3 \n', ...
    'Xn 6 1 0 0 lvtnfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtn p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Xp 6 1 3 3 lvtpfet L=l nfin=w m=1 nf=1 par=1 par_nf=1 asej=''594e-18*w+(w*int(500e-3))*594e-18'' adej=''(w*(1+int(0.0)))*594e-18'' psej=''(((1*w)*2)*54e-9+11e-9)+(((w*int(500e-3))*1)*2)*54e-9'' pdej=''(((w*(1+int(0.0)))*1)*2)*54e-9'' pdevdops=1 pdevlgeos=1 pdevwgeos=1 psw_acv_sign=1 plnest=1 pldist=1 plorient=0 cpp=78e-9 fpitch=48e-9 xpos=-99 ypos=-99 ptwell=0 sca=0 scb=0 scc=0 pre_layout_local=1 ngcon=1 p_vta=Vtp p_la=0 u0mult_fet=1 lle_sa=71e-9 lle_sb=71e-9 lle_rxrxa=78e-9 lle_rxrxb=78e-9 lle_rxrxn=192e-9 lle_rxrxs=192e-9 lle_pcrxn=65e-9 lle_pcrxs=65e-9 lle_nwa=2e-6 lle_nwb=2e-6 lle_nwn=192e-9 lle_nws=192e-9 lle_ctne=0 lle_ctnw=0 lle_ctse=0 lle_ctsw=0 lle_sctne=0 lle_sctnw=0 lle_sctse=0 lle_sctsw=0 lrsd=27e-9 dtemp=0 l_shape=0 l_shape_s=0 nsig_dop1=0 nsig_dop2=0 nsig_dibl=0 nsig_pc=0 nsig_rx=0 fc_index=0 fc_sigma=3 analog=-1 nf_pex=1\n', ...
    'Rl 6 2 rl\n', ...
    'Cl 2 0 cl\n', ...
    '.ends\n']);

    fprintf(fid,'xinv100 100 1 vdr inv2 w=w1\n');
    % fprintf(fid,'xinv100 100 1 vdr inv2 w=w2\n');
    for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i %i %i vdd inv1 w=w%i\n',i,i,i+1,i);
        else
            fprintf(fid,'xinv%i %i %i vdd inv2 w=w%i\n',i,i,i+1,i);
            if(i<n)
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i''\n',100+i-1,i+1,200+i,i-1);
                % for (j = 1:1:ceil(fanout-1))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            else
                fprintf(fid,'xinv%i %i %i vdr inv1 w=''fanout*w%i+w%i''\n',100+i-1,i+1,200+i,i-1,i-1);
                % for (j = 1:1:ceil(fanout))
                %     fprintf(fid,'xinv%i %i %i vdr inv1 w=w%i\n',1000+(i*10)+j-1,i+1,200+i,i-1);
                % end
            end
        end
    end

    if(flag_var)
        fprintf(fid,'.subckt inv1_s 1 2 3 dvt=0 w=1u\nXn1 7 4 0 0 d_nfet W=w L=l\nXn2 2 8 7 0 d_nfet W=w L=l\nXp1 2 9 3 3 d_pfet W=w L=l\nXp2 2 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nVtno 8 3 ''dvt+vshift'' \nVtpo 9 3 ''-dvt-vshift''\n.ends\n');
        fprintf(fid,'.subckt inv2_s 1 2 3 dvt=0 w=1u\nXn 6 4 0 0 d_nfet W=w L=l\nXp 6 5 3 3 d_pfet W=w L=l\nVtn 4 1 ''dvt+vshift''\nVtp 5 1 ''-dvt-vshift'' \nRl 6 2 rl\nCl 2 0 cl\n.ends\n');

        fprintf(fid,'xinv100_s 100_s 1_s vdr_s inv2_s w=''w1-%4.4f'' dvt=%3.3f\n',dw1,dvt);
        for i=1:1:n
        if(mod(i,2))
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv1_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
        else
            fprintf(fid,'xinv%i_s %i_s %i_s vdd_s inv2_s w=''w%i-%4.4f'' dvt=%3.3f\n',i,i,i+1,i,dw1,dvt);
            if(i<n)
                fprintf(fid,'xinv%i_s %i_s %i_s vdr_s inv1_s w=''fanout*w%i-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,dw1,dvt);
            else
                fprintf(fid,'xinv%i_s %i_s %i_S vdr_s inv1_s w=''fanout*w%i+w%i-%4.4f-fanout*%4.4f'' dvt=%3.3f\n',100+i-1,i+1,200+i,i-1,i-1,dw1,dw1,dvt);
            end
            
        end
        end
        fprintf(fid,'Vdr_s vdr_s 0 DC=''vd*%2.2f''\nVd_s vdd_s 0 DC=''vd*%2.2f''\nVpulse_s 100_s 0 PULSE(0,''vd*%2.2f'',10n,10p,10p,''clock/2'',clock)\n',p_slow,p_slow,p_slow);
    end 

    % fprintf(fid,'Vdr vdr 0 DC vd\nVd vdd 0 DC vd\nVpulse 100 0 PULSE(0,vd,10n,10p,10p,pulse,clock)\n.tran ''clock/1000'' ''clock*3+10n'' \n');
    fprintf(fid,['Vdr vdr 0 DC vd\n', ...
    'Vd vdd 0 DC vd\n', ...
    'Vpulse 100 0 PULSE(0,vd,10n,10p,10p,''clock/2'',clock)\n', ...
    '.tran ''clock/10000'' ''clock*5+10n'' \n']);
    % fprintf(fid,'Vdr vdr 0 DC vd\n Vd vdd 0 DC vd\n');
    % fprintf(fid, ['Vpwl 100 0 pwl(''0n'',0V \\ \n', ...
    % '''10.00n+0*clock'',0V ''10.01n+0*clock'',''vd'' ''10.01n+1*pulse+0*clock'',''vd'' ''10.02n+1*pulse+0*clock'',0V \\ \n', ...
    % '''10.02n+1*clock'',0V ''10.03n+1*clock'',''vd'' ''10.03n+1*pulse+1*clock'',''vd'' ''10.04n+1*pulse+1*clock'',0V \\ \n', ...
    % '''10.04n+2*clock'',0V ''10.05n+2*clock'',''vd'' ''10.05n+1*pulse+2*clock'',''vd'' ''10.06n+1*pulse+2*clock'',0V \\ \n', ...
    % '''10.06n+3*clock'',0V ''10.07n+3*clock'',''vd'' ''10.07n+1*10000n+3*clock'',''vd'' ''10.08n+1*10000n+3*clock'',0V \\ \n', ...
    % '''10.08n+2*10000n+3*clock'',0V ''10.09n+2*10000n+3*clock'',0V \n']);
    % fprintf(fid,'.tran ''clock/1000'' ''10.09n+2*10000n+3*clock'' \n');

    if(flag_var)
        fprintf(fid,'.measure tran t1 trig V(1_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i_s) val=''%2.2f*vd/2'' TD=''10n+clock/4'' %s=2\n',p_slow,n+1,p_slow,edge);
    else
        fprintf(fid,'.measure tran t1 trig V(1) val=''vd/2'' TD=''10n+clock/4'' rise=2 targ V(%i) val=''vd/2'' TD=''10n+clock/4'' %s=2\n',n+1,edge);
        fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))'' FROM=''clock*2+10n'' TO=''clock*5+10n''\n');
        fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.2+1)'' TO=''10n+clock*(1.4+1)''\n');
        fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))'' FROM=''10n+clock*(1.7+1)'' TO=''10n+clock*(1.9+1)''\n');

        % fprintf(fid,'.measure avgdelay param=''(t1+t2)*0.5''\n');
        % fprintf(fid,'.measure tran avgpower avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10n'' TO=''clock*3+10n''\n');
        % fprintf(fid,'.measure leakpower param=''(leakpower1+leakpower2)*0.5''\n\n');

        % fprintf(fid,'.measure tran t1 trig V(1) val=''vd*0.5'' TD=0 rise=3 targ V(%i) val=''vd*0.5'' TD=0 rise=3\n',n+1);
        % fprintf(fid,'.measure tran t2 trig V(1) val=''vd*0.5'' TD=0 fall=3 targ V(%i) val=''vd*0.5'' TD=0 fall=3\n',n+1);
        % fprintf(fid,'.measure tran leakpower1 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.07n+3*clock+500n'' TO=''10.07n+1*10000n+3*clock-500n''\n');
        % fprintf(fid,'.measure tran leakpower2 avg par''(-1)*abs(p(Vd))+(-1)*abs(p(Vpwl))'' FROM=''10.08n+1*10000n+3*clock+500n'' TO=''10.08n+2*10000n+3*clock-500n''\n');
    end
    % fprintf(fid,'\n.measure tran avgpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*0.3+10n'' TO=''clock*2.3+10n''\n');
    % fprintf(fid,'.measure tran leakpower avg par''(-1)*abs(p(Vd))+p(Vpulse)'' FROM=''clock*1.3+10n+10p'' TO=''clock*1.4+10n+10p''\n');
    % fprintf(fid,'.measure VTH\n');
    fprintf(fid,'.ic ');
    fprintf(fid,'v(100)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i)=0 ',i);
        else
            fprintf(fid,'v(%i)=''vd'' ',i);
        end
    end
    if(flag_var)
        fprintf(fid,'v(100_s)=0 ');
    for i=1:n
        if(mod(i,2)==0)
            fprintf(fid,'v(%i_s)=0 ',i);
        else
            fprintf(fid,'v(%i_s)=''vd*%2.2f'' ',i,p_slow);
        end
    end
    end
    fprintf(fid,'\n');
    fprintf(fid,'.end');
    fclose(fid);
    unix(hspicerun);

    fid=fopen(resultname, 'r');
    if (fid ~= -1)
        x = [];
        for (j=1:(5+1))
            tline=fgetl(fid);
        end
        [t r] = strtok(tline);
        if(~ischar(t))
            f4=1;
            Delay2=38;
        else
            [Delay2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Delay2=38;
        end
        if(Delay2<=0)
            f4=1;
        end
        
        [t r] = strtok(r);
        if(~ischar(t))
            f4=1;
        else
            [Energy2,status]=str2num(t);
        end
        if(status==0)
            f4=1;
            Energy2=38;
        end
                
        if(Energy2>=0  || Energy2>Leakage2/2)
            f4=1;
            disp('energy2 >=0 or energy2 < leakage');
            disp(Energy2);
            disp(Leakage2);
        else
            if(Energy2>Leakage2)
                Energy2=Leakage2;
            end
        end
    else
        f4=1;
    end
    fclose(fid);
    end
end

% delay calibration factor
Delay2 = Delay2/calib_delay;

% calibrating leak and dyn pwr results before sending back to main prog
Dynamic2 = (-Energy2+Leakage2)/calib_dyn;       % Energy2: P_average in slide(20230316 - Rhesa_s LTLT Progress Update)
Leakage2 = Leakage2/calib_leak;
% Energy2=Dynamic*6*Delay1/2*activity; 
Energy2=Dynamic2*10*Delay1*0.5*activity;       % Delay1*0.5: pulse width in slide(20230316 - Rhesa_s LTLT Progress Update)
fprintf('VDD: %4.4f - VTN: %4.4f - VTP: %4.4f\n',Vd, Vtn, Vtp);