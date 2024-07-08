fclose all;
clearvars;
tic;
fprintf('-----SEPVT PROGRAM----- CALIB\n');

%%%%% TUNING SWITCHES
% activity factor
acti        = 0.1;
% temperature setting: lo (-196.15) or md (+25C) or hi (+80)
tempe_set   = 'lo';
% iterations
kval        = 30;
k_extend = 0;
testingcalib = 0;
trynum = 3;

%%%% AGING
% aging = 0 disable 1 enabled
aging       = 1;

%%%% FIXED VOLTAGE
% fixed voltage
fixed_V     = 0;

%%%% SIZING OPTION
%0: constant vdd or 1: optimally picking vdd
flag_vdd    = 1;
%0:constant vtn or 1: optimally picking vt
flag_vtn    = 1; 
%0:constant vtp or 1: optimally picking vt
flag_vtp    = 1;
%0 : no plots, 1: show optimizing plot
flag_plot   = 1;
%0: constant width or 1: optimally picking widths
flag_w = 1;
%0：not consider vth variation, 1: consider vth variation same values as RT for LT, 2: consider vth variation same percentage (3σ/Vth0) as RT for LT
vt_var_mode = 0;
flag_IR_drop = 0;


if(vt_var_mode == 0)
    vtn_var_val = 0;
    vtp_var_val = 0;
    vtn_percent = 0;
    vtp_percent = 0;
elseif(vt_var_mode == 1)
    vtn_var_val = 0.13;
    vtp_var_val = 0.13;
    vtn_percent = 0;
    vtp_percent = 0;
elseif(vt_var_mode == 2)
    vtn_var_val = 0;
    vtp_var_val = 0;
    vtn_percent = 0.05;
    vtp_percent = 0.05;
end
%0: new start, 1: continue from old Pareto points
flag_continue = 0;

%%%% VARIATION SWITCHES
% using temperature = 1, not using = 0
temp = 1; 
flag
%%%% WIRE MODEL
% 1: enable wire model, 0: Constant wire model
flag_wire   = 0;
% temp wire
temp_wire   = 1;
% wire accuracy mode: 0) RC = 0, 1) std RC wire with average resistance from design
% wire accuracy mode: 2) R accurate model with constant C, 3) C accurate model with constant R
accurate_wire = 2;
wire_width_scale = 1;

% VTSELECT MVT or SVT
% vtSelect = '1vt';
vtSelect = '1vt';
% VTSELECT MVD or SVD
vdSelect = '1vd';
% vdSelect = '2vd';

if (flag_wire == 0) %nowire model still didn't use wire per bin,,, (diff path 0)
    flagWireStat = 'cwm';
else
    flagWireStat = 'ewm';
end

if (aging == 0)
    agingStat = 'na';
else
    agingStat = 'a';
end

% Temperature selector
if (tempe_set == 'hi')
    tempe   = 85;
elseif (tempe_set == 'lo')
    tempe   = -196;
elseif (tempe_set == 'md')
    tempe = 25;
end


% PRINT FOR FINAL CHECK
fprintf('> ACTIV: %.2f - TEMP SET: %i - kval: %i\n',acti, tempe, kval);
fprintf('> WIRE ACCURACY MODE: %i - WIRE MODEL: %s - AGING: %s\n',accurate_wire, flagWireStat, agingStat);
% fprintf('> WIRE MODEL: %s - AGING: %s\n',flagWireStat, agingStat);
fprintf('> VT SELECTOR: %s\n',vtSelect);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters for Pareto optimization%%%%%%%%%%%%%%%%%%%%%
% alpha per degree celsius
alpha       = 0.004041; 
% fin pitch
width_eff_mult   = 0.075;
% optimizing selections 1) delay and power, 2) delay and area, and 3) power and area
mode        = 1;
%the third metric limitation: [max araa] in mode 1, [max power] in mode 2, [max delay] in mode 3
max_limit_metric = 1e12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start_t=cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read logic depth histogram 
LDH             = dlmread('LDH_M3_1ns_prelayout1.txt'); 
% Read aging information
deltav          = readmatrix('aging.csv');
% Total number of gates in design
total_gates     = 68603;
% Critical path depth
depth_critical  = max(LDH);
% Number of bins divided in LDH (Must be sequency number and start from 1)
total_stage = 10;
optimizing_stage = 1:total_stage;
% optimizing_stage = [1,2,3,4,5,6,7,8,9,10];
% optimizing_stage = [1,2,3,4,5];
% Divide multiple vdds into different bins, e.g. [1 2; 3 4] bin 1&2 use Vdd1, bin3&4 use vdd2
if (strcmp (vdSelect, '1vd') == 1)
    vdd_allocation = 1:total_stage;
    % vdd_allocation = [1,2,3,4,5,6,7,8,9,10];
elseif (strcmp (vdSelect, '2vd') == 1)
    vdd_allocation = [1:round(total_stage/2); round(total_stage/2)+1:total_stage, zeros(mod(total_stage,2)==1)];
    % vdd_allocation = [1,2,3,4,5,6,7; 8,9,10,0,0,0,0];
end
% Divide mutiple vts into different bins
% SVT
if (strcmp (vtSelect, '1vt') == 1)
    vtn_allocation = 1:total_stage;
    % vtn_allocation = [1,2,3,4,5,6,7,8,9,10];
elseif (strcmp (vtSelect, '2vt') == 1)
    vtn_allocation = [1:round(total_stage/2); round(total_stage/2)+1:total_stage, zeros(mod(total_stage,2)==1)];
    % vtn_allocation = [1,2,3,4,5,6,7; 8,9,10,0,0,0,0];
end
vtp_allocation = vtn_allocation;

% initial boundry for 2vt and 2vd
boundry_ini = 1;

%%% divide multiple HGI devices into bins[1 2;3 0]bin 1 and bin 2 use model 1, bin3 use model 2
devicemodel_allocation = 1:total_stage;
% devicemodel_allocation = [1,2,3,4,5];
if(length( devicemodel_allocation(:,1)) > 1)
    flag_HGI = 1;
else
    flag_HGI = 0;
end
n_bins=length(optimizing_stage);

% define hspice time limit in seconds
hspice_time_limit = ceil((n_bins/4)*360);


stage_weight;
if(mode == 3)
    max_limit_metric = max_limit_metric / delay_weight; %consistent with delay calculation
end

n_vdd= length(vdd_allocation(:,1));
n_vtn= length(vtn_allocation(:,1));
n_vtp= length(vtp_allocation(:,1));
fprintf('> VT COUNT: %i - NUMBER OF BINS: %i\n' , n_vtn, n_bins);

%%variation knobs
% Capacitance out
Cout        = 100e-15;%
% Temperature pass variable; tempe var to temp_v
temp_v      = tempe;
%0: variation is not enabled, 1: Enable variation
flag_var    = 0;
% percentage of low vdd (suffered from variation) in nominal vdd
p_slow      = 0.95; 
% threshold voltage increasing due to variation, unit [V]
dvt         = 0.05;
% gate sizes variation in percentage
dw1         = 0.001;

%%%%%wire model*
%%% channel length m --> that's why we should use c14 instead of c16
L           = 14e-9; 
%%% activity factor pass variable
activity    = acti;
%% multiplier for wire load 
wireload    = 1;
% avg wire length per net in um
wire_len_per_net= 6.5300;
% total net in the design
% net_count= 40239;
%% DATA from design
% R total from spef (OHM)
% R_total = 7.72090903E+06;
% C total from spef (PF)
C_total = 5.02262046E+01;

total_length = 262758.0910;

if temp_wire
    tempeWire = tempe;
else
    tempeWire = 25;
end

if (accurate_wire == 0)
    % no resistance
    Ri = resistive(tempeWire,accurate_wire,wire_len_per_net,wireload);
    Ci = (C_total*(10^(-12))/(32151*wire_len_per_net)) * wireload; 
elseif (accurate_wire == 1 || accurate_wire == 2)
    Ri = resistive(tempeWire,accurate_wire,wire_len_per_net,wireload);
    % capacitance for 1um wire (avg capacitance per net/avg wire length per net)
    Ci = (C_total*(10^(-12))/(32151*wire_len_per_net)) * wireload; 
elseif (accurate_wire == 3)
    Ri = resistive(tempeWire,accurate_wire,wire_len_per_net,wireload);
    layer_c(1) = capacitance(522, 'PC', 'M1', 'M2');
    layer_c(2) = capacitance(46184, 'M1', 'M2', 'M3');
    layer_c(3) = capacitance(61255, 'M2', 'M3', 'C4');
    layer_c(4) = capacitance(64544, 'M3', 'C4', 'C5');
    layer_c(5) = capacitance(68716, 'C4', 'C5', 'K1');
    layer_c(6) = capacitance(21554, 'C5', 'K1', 'K2');
    Ci = (sum(layer_c,"all")*(10^(-12))/(40239*wire_len_per_net)) * wireload;
% elseif (accurate_wire == 4)
%     % determine R first with mode 2
%     Ri = resistive(tempeWire,2,wire_len_per_net,wireload);
%     Ci = (C_total*(10^(-12))/(32151*wire_len_per_net)) * wireload;
end

% WIRE SCALING: Scale R by R/wire_width_scale and C by C/wire_width_scale
Ri = Ri*wire_width_scale;
Ci = Ci/wire_width_scale;


fprintf('Tau Value: %d\n', Ri*Ci);
fprintf('Ri Value: %d\n', Ri);
fprintf('Ci Value: %d\n', Ci);


%%% design area um^2
input_total_area    = 0.1452 * total_gates;
%Fan-out (w/o -1, raw from spef)
fanout              = ceil(3.81636652);

%%% calibration factors
% calib factor for dynamic power; calibrated pwr = raw_dyn_pwr/calib_dyn
calib_dyn = 1.7392 / 1.7707
% calib factor for leakage power; calibrated pwr = raw_leak_pwr/calib_leak
calib_leak = 1.4765 / 1.7707
% calib factor for delay
calib_delay = 1.1733

%timing upper bound
delay_upper_bound   = 0.7e-9;
%% timing lower bound
delay_lower_bound   = 100e-9;
% design power upper bound
power_upper_bound   = 1e-1;
% design power lower bound 
power_lower_bound   = 1e-8;

devicemodel(1).subcktname=['/app/design/puneet/projects/Shaodiwang/model/standard_model/TFET.sp'];%% model for device1 should include both NMOS and PMOS
devicemodel(2).subcktname=['/w/design/puneet/projects/Shaodiwang/model/standard_model/SOI.sp']; % model for device 2

% parameter_initial, upper bound, lower bound -> [Vdd, VtN, VtP, Gate Size (nm)]
% sweep starts at the middle values of each lower and upper bound.

if (testingcalib == 1)
    % params for testing
    ub=[0.81, 0, 0, 0.150];
    lb=[0.81, 0, 0, 0.150];
    parameter_initial = [0.9 0 0 0.150]
elseif (testingcalib == 0)
    if (tempe == 85)
        ub=[1.0, 0.66, 0.66, 4.500];
        lb=[0.05, -0.5, -0.5, 0.001];
        if (trynum == 1)
            parameter_initial = [0.6 -0.18 -0.2 0.150]
        elseif (trynum == 2)
            parameter_initial = [0.5 -0.14 -0.16 0.150]
        elseif(trynum == 3)
            parameter_initial = [0.42 -0.09 -0.13 0.150]
        elseif(trynum == 4)
            parameter_initial = [0.8 0 0 0.150]
        end
    end
    if (tempe == -196)
        ub=[1.0, -0.2400, -0.2400, 4.500];
        lb=[0.08, -0.3580, -0.3580, 0.100];
        if (trynum == 1)
            parameter_initial = [0.2500 -0.2800 -0.3200 0.150]
        elseif (trynum == 2)
            parameter_initial = [0.2500 -0.2800 -0.3200 0.150]
        elseif(trynum == 3)
            parameter_initial = [0.2500 -0.2800 -0.3200 0.150]
        elseif(trynum == 4)
            parameter_initial = [0.2500 -0.2800 -0.3200 0.150]
        end
    end
end

devicename=['Bulk']; % desgn number
number_sim = 1; % number of simulations performed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% You may not need to change the following codes%%%%%%%%
save('Design_information.mat','devicemodel_allocation','optimizing_stage','Ri','Ci','input_total_area','flag_wire','fanout','wireload','activity','delay_weight','optimizing_stage_weight','flag_var','p_slow','dvt','dw1','vdd_allocation','vtn_allocation','vtp_allocation','LDH','L','devicemodel','devicename','mode');
dwli=0.3;  %%% partial step in times unit
dwui=0.5;
dvi=0.002;  %%% partial voltage
Rti=20;  %%%% Trust region radius
Rt_th=4;
pd=0.5; %%% decreation speed on Rt
glb_pd=0.99; %%% global decreation speed on Rt
istart=0;
ivt=0;
    if(flag_HGI)
        area_model = fitting_HGI_interconnect_model('./nangate_TFET.data','./inv_TFET.data','./nangate.data','./inv.data' );
    else
        area_model = fitting_interconnect_model('./nangate.data','./inv.data');
    end
[min_area,temp_scale] = chip_area_model(devicemodel_allocation, area_model, lb(4)*ones(3+ 2*length(optimizing_stage), length(optimizing_stage)),optimizing_stage,optimizing_stage_weight );
[max_area,temp_scale] = chip_area_model(devicemodel_allocation, area_model, ub(4)*ones(3+ 2*length(optimizing_stage) , length(optimizing_stage)),optimizing_stage,optimizing_stage_weight );
area_scale = (max_area - min_area) /5; %%%This scale is used to scale area to the same range of delay and power
if(mode == 1)
    max_limit_metric = max_limit_metric / area_scale;
end

for i_sim = 1: number_sim
    for ivd=1:1
        while( 1)
            if(ivt==1)
                ivt=0;
                break;
            end
            ivt=ivt+1;
            for iw=1:1
                X=[];
                for i_model = 1:length(devicemodel_allocation(:,1))
                    for i_ss = 1:length(devicemodel_allocation(i_model,:))
                        i_stage = devicemodel_allocation(i_model,i_ss);
                        if(i_stage ==0)
                            break;
                        end
                        n=2*optimizing_stage(i_stage);
                        X(1,i_stage)=parameter_initial(i_model,1)+0.2*(ivd-1);
                        X(2,i_stage)=parameter_initial(i_model,2)+0.05*(ivt-1); %    this is for VtN
                        X(3,i_stage)=parameter_initial(i_model,3)+0.05*(ivt-1); %    new --> this is for VtP
                        ub_i(1:3,i_stage)=ub(1:3);
                        lb_i(1:3,i_stage)=lb(1:3);
                        for iwtt=1:n/2
                                X(3+iwtt*2-1,i_stage)=parameter_initial(i_model,4);
                                X(3+iwtt*2,i_stage)=(1.0+0.5*iw)*X(3+iwtt*2-1,i_stage);
                                % X(3+iwtt*2,i_stage)=X(3+iwtt*2-1,i_stage)*((fanout-1)/3.59);
                                ub_i(3+2*iwtt-1:3+2*iwtt,i_stage)=ub(4);
                                lb_i(3+2*iwtt-1:3+2*iwtt,i_stage)=lb(4);
                        end
                    end
                end



            % You can specify the initial point here to accelerate the optimization 
            if (i_sim == 1)
                if (tempe_set == 'hi')
                    if (length(vdd_allocation(:,1)) == 1)
                        X(1,vdd_allocation(1,:)) = 0.9500;      % vdd1
                    else
                        % for 2vdd,
                        X(1,vdd_allocation(1,:)) = 0.4500;      % vdd1
                        X(1,vdd_allocation(2,1:max(optimizing_stage)-max(vdd_allocation(1,:)))) = 0.9500;      % vdd2
                    end
                    if (length(vtn_allocation(:,1)) == 1)
                        X(2,vtn_allocation(1,:)) = -0.1800;     % vtn1
                    else
                        % for 2vtn,
                        X(2,vtn_allocation(1,:)) = -0.1500;     % vtn1
                        X(2,vtn_allocation(2,1:max(optimizing_stage)-max(vtn_allocation(1,:)))) = -0.1800;     % vtn2
                    end
                    % for 2vtn,
                    if (length(vtp_allocation(:,1)) == 1)
                        X(3,vtp_allocation(1,:)) = -0.2280;     % vtp1
                    else
                        % for 2vtp,
                        X(3,vtp_allocation(1,:)) = -0.1880;     % vtp1
                        X(3,vtp_allocation(2,1:max(optimizing_stage)-max(vtp_allocation(1,:)))) = -0.2280;     % vtp2
                    end
                else % tempe_set == 'lo'
                    if (length(vdd_allocation(:,1)) == 1)
                        X(1,vdd_allocation(1,:)) = 0.2672;      % vdd1
                    else
                        % for 2vdd,
                        X(1,vdd_allocation(1,:)) = 0.2672;      % vdd1
                        X(1,vdd_allocation(2,1:max(optimizing_stage)-max(vdd_allocation(1,:)))) = 0.2672;      % vdd2
                    end
                    if (length(vtn_allocation(:,1)) == 1)
                        X(2,vtn_allocation(1,:)) = -0.3030;     % vtn1
                    else
                        % for 2vtn,
                        X(2,vtn_allocation(1,:)) = -0.3030;     % vtn1
                        X(2,vtn_allocation(2,1:max(optimizing_stage)-max(vtn_allocation(1,:)))) = -0.3030;     % vtn2
                    end
                    % for 2vtn,
                    if (length(vtp_allocation(:,1)) == 1)
                        X(3,vtp_allocation(1,:)) = -0.3030;     % vtp1
                    else
                        % for 2vtp,
                        X(3,vtp_allocation(1,:)) = -0.3030;     % vtp1
                        X(3,vtp_allocation(2,1:max(optimizing_stage)-max(vtp_allocation(1,:)))) = -0.3030;     % vtp2
                    end
                end
            end
            

            if(i_sim >1)
                for i_stage=1:n_bins
                    X(1,i_stage) = next_X.X(1,i_stage);
                    X(2,i_stage) = next_X.X(2,i_stage);
                    X(3,i_stage) = next_X.X(3,i_stage);
                    for i=4:2*optimizing_stage(i_stage)+3
                        X(i,i_stage) = next_X.X(i,i_stage);
                    end
                end
            end
            X(1:3,n_bins+1)=0;
            %    X
            DB(1).X=X;
            istart=istart+1;
            Xc=[];
            J=[];
            nof=[];
            DB(1).X=X;

%%% enviroment specification

Rt=Rti;  %%%% Trust region radius
DB(1).X(1,n_bins+1)=Rt;
sigma=1.5;%%%%minimum Rt
Resolution=0.2;
% fprintf('mark here\n')''

%%%%parameter calculationr
spname      = ['acti',num2str(activity),'_',num2str(n_vtn),'vt_',num2str(temp_v),'C_k',num2str(kval),'_wireacc',num2str(accurate_wire),'_',flagWireStat,'_',agingStat,'_',num2str(n_bins),'_try',num2str(trynum),'_wirescale',num2str(wire_width_scale),'.sp'];
resultname  = ['acti',num2str(activity),'_',num2str(n_vtn),'vt_',num2str(temp_v),'C_k',num2str(kval),'_wireacc',num2str(accurate_wire),'_',flagWireStat,'_',agingStat,'_',num2str(n_bins),'_try',num2str(trynum),'_wirescale',num2str(wire_width_scale),'.mt0'];
matname     = ['acti',num2str(activity),'_',num2str(n_vtn),'vt_',num2str(temp_v),'C_k',num2str(kval),'_wireacc',num2str(accurate_wire),'_',flagWireStat,'_',agingStat,'_',num2str(n_bins),'_try',num2str(trynum),'_wirescale',num2str(wire_width_scale),'.mat'];
logname     = ['acti',num2str(activity),'_',num2str(n_vtn),'vt_',num2str(temp_v),'C_k',num2str(kval),'_wireacc',num2str(accurate_wire),'_',flagWireStat,'_',agingStat,'_',num2str(n_bins),'_try',num2str(trynum),'_wirescale',num2str(wire_width_scale),'.log'];
hspicerun   = ['hspice -i ',spname,' > ',logname];
txtname     = [num2str(istart),'_','acti',num2str(activity),'_',num2str(n_vtn),'vt_',num2str(temp_v),'C_k',num2str(kval),'_wireacc',num2str(accurate_wire),'_',flagWireStat,'_',agingStat,'_',num2str(n_bins),'_try',num2str(trynum),'_wirescale',num2str(wire_width_scale),'.txt'];
Xtxtname    = [num2str(istart),'_','X_','acti',num2str(activity),'_',num2str(n_vtn),'vt_',num2str(temp_v),'C_k',num2str(kval),'_wireacc',num2str(accurate_wire),'_',flagWireStat,'_',agingStat,'_',num2str(n_bins),'_try',num2str(trynum),'_wirescale',num2str(wire_width_scale),'.txt'];

edge=['rise'];

if (flag_continue)
    % rename the previous file
    oldname = [matname(1:end-4), '_old.mat'];
    status = movefile(matname, oldname);
    % check if the renaming is successful
    if status
        disp('rename success.');
    else
        disp('rename failed.');
    end
end

%%%%%wireload
%%%%%calculate the approximate delay

[Area_i,G_A] = chip_area_model(devicemodel_allocation, area_model, lb(4)*ones(3+ 2*length(optimizing_stage), length(optimizing_stage)),optimizing_stage,optimizing_stage_weight );
if(flag_wire)
    [R,C0]=interconnect_model( devicemodel_allocation, area_model,X,optimizing_stage,optimizing_stage_weight,input_total_area,Ri,Ci);
else
    C0=Ci;
    R=Ri;
end


%%%%initial point
color_ng=['r','y'];
k=1;
nof(k)=1;

f2=0;
f4=0;
f5=0;
current_E=0;
current_L=0;
if (flag_continue)
    [DB, J, J_DB, X] = start_from_pre_results(oldname, n_bins);
else
    initial_measure_simulation_fast;
    % This simulation will calculate the detail delay, dynamic power and leakage power.
    detail_delay_and_power_simulation_starting_point_fast;
    % J is the reported matrix
    %Delay J(:,1), power (total) J(:,2), leakage J(:,4), area J(:,5) and id is J(:,3)
    J(1,1)=log10(max_D);
    J(1,2)=log10(Power_total);
    J_DB(1,1:n_bins)=Delay_each_stage; %Stage delay
    J(1,3)=1;
    J(1,4)=log10(current_L);
    J(1,5) = Area_i / area_scale;
    fprintf('%.2f\n',current_L/Power_total*100);
end
nof(k) = size(J,1);

f=0;%%%flag 1 
%%%%pareto search starts
metric_scale = 5; % scale max delay or power up in gradient descent to avoid model error, this value will be automatically updated

% plot all current Pareto points
if(flag_plot)
    switch mode
        case 1 
            plot_x = log10(delay_weight)+J(:,1);
            plot_y = J(:,2);
        case 2
            plot_x = log10(delay_weight)+J(:,1);
            plot_y = J(:,5) ;
        case 3
            plot_x = J(:,5);
            plot_y = J(:,2);
    end
    h=plot(plot_x,plot_y,'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
    title(vtSelect + " " + vdSelect);
    drawnow;
end
for i=1:nof(k)
    DB(i).X(1,n_bins+1)=Rti;
    DB(i).X(2,n_bins+1)=0;
end

while (1)
    updateinformation=['k=',num2str(k),' nof(k)=',num2str(nof(k)),' metric scale is ',num2str(metric_scale)];
    disp(updateinformation);
    if(k>kval)
        break;
    end

    for i=1:nof(k)
        DB(i).X(1,n_bins+1)=DB(i).X(1,n_bins+1)*glb_pd;
    end

    disp(updateinformation);
    % --- below here is still unedited ---
    if (nof(k)>2)
        nt=nof(k);
        d=zeros(nt,2);
        for i=2:nt-1
            switch mode
                case 1
                    index_com = [1 2];
                case 2
                    index_com = [1 5];
                case 3
                    index_com = [2 5];
            end
            d(i,1) = sqrt( (J(i,index_com)-J(i-1,index_com))*((J(i,index_com)-J(i-1,index_com))'))+sqrt( (J(i+1,index_com)-J(i,index_com))*((J(i+1,index_com)-J(i,index_com))'));
            d(i,2) = i;
            % The slope of the line with the left and right points
            d(i,3) = (J(i-1,index_com(2))-J(i,index_com(2)))/(J(i-1,index_com(1))-J(i,index_com(1)) - 1e-12) - (J(i,index_com(2))-J(i+1,index_com(2)))/(J(i,index_com(1))-J(i+1,index_com(1)) - 1e-12);
        end
        d(1,1) = 2*sqrt( (J(2,index_com)-J(1,index_com))*((J(2,index_com)-J(1,index_com))'));
        d(1,2) = 1;
        d(1,3) =  0 - ((J(1,index_com(2))-J(2,index_com(2)))/(J(1,index_com(1))-J(2,index_com(1)) - 1e-12));
        d(nt,1) = 2*sqrt( (J(nt,index_com)-J(nt-1,index_com))*((J(nt,index_com)-J(nt-1,index_com))'));
        d(nt,2) = nt;
        d(nt,3) = (J(nt-1,index_com(2))-J(nt,index_com(2)))/(J(nt-1,index_com(1))-J(nt,index_com(1)) - 1e-12);
        
        if (k <= k_extend)
            if(mod(k,2)==0)
                l=d(nt,2);            % First 10 iterations, choose the most right point as the starting point to extend the Pareto front
            else
                l=d(1,2);           % First 10 iterations, choose the most left point as the starting point to extend the Pareto front
            end
        else
            d=sortrows(d,1);
            l=d(nt,2);              % Choose the one with the largest distance to its neighbors as the starting point
        end
        nt_temp = nt;
        while (d(nt_temp,3)>0 || J(d(nt_temp,2), 1) + log10(delay_weight) > -8.85 || J(d(nt_temp,2), 1) + log10(delay_weight) < -8.95)
            nt_temp=nt_temp-1;
            l=d(nt_temp,2);
            fprintf('Switch the starting point to the neighbor point\n');
        end
        if(k==1)
            l=646;
        end


        Xc=DB(l).X;
        while (Xc(2,n_bins+1)==1)
            nt=nt-1;
            if(nt==0)
                f=1;
                disp('all are selected as Xc');                 
                break;                
            end
            l=d(nt,2);
            Xc=DB(l).X;
        end
        if(f==1)
            break;
        end
        
        switch mode
            case 1 
            if( J(nof(k),1)>log10(delay_lower_bound/delay_weight) && (J(1,1)<log10(delay_upper_bound/delay_weight) || J(1,5) > (max_limit_metric/1.5) /area_scale )&& nt>80 && d(nt,1)<Resolution)
                f=1;
                disp('finish after resolution is satisfied');
                break;
            end
            case 2
                if( J(nof(k),5) < min_area / area_scale * 1.5  && (J(1,1) > max_area / area_scale/2 || J(1,2) > log10(max_limit_metric/2) ) && nt>80 && d(nt,1)<Resolution)
                    f=1;
                    disp('finish after resolution is satisfied');
                    break;
                end
            case 3
            if( (J(nof(k),2)>log10(power_lower_bound/delay_weight) || J(nof(k),1) > log10(max_limit_metric/1.5 ) )&& J(1,2)<log10(power_upper_bound/delay_weight) && nt>80 && d(nt,1)<Resolution)
                f=1;
                disp('finish after resolution is satisfied');
                break;
            end
        end

        Delay_i(1:n_bins)=J_DB(l,1:n_bins);
        Rt=DB(l).X(1,n_bins+1);
        DB(l).X(1,n_bins+1)=Rt*pd;
        
        if(Rt<sigma)
            DB(l).X(2,n_bins+1)=1;
        end
    else 
        l=1;
        if(DB(l).X(2,n_bins+1)==1)
            l=2;
            if(nof(k)==1)
                disp('break when there is no front can be chosen as Xc');
                break;
            end
            if(DB(l).X(2,n_bins)==1)
                disp('break when there is no front can be chosen as Xc~');
                break;
            end
        end
        Xc=DB(l).X;
        
        Delay_i(1:n_bins)=J_DB(l,1:n_bins);
        Rt = DB(l).X(1,n_bins+1);
        DB(l).X(1,n_bins+1)=Rt*pd;
        
        if(Rt<sigma)
            DB(l).X(2,n_bins+1)=1;
        end
    end
    
    Delay1 = max(J_DB(l,1:n_bins));
    % Reselect allocation boundry for multi-voltage configuration
    if (strcmp(vdSelect, '1vd') == 1 && strcmp(vtSelect, '1vt') == 1)
        fprintf('1vd 1vt, no need to reselect the allocation boundry\n');
    else
        reselect_allocation_boundry_vd_vt_fast;
    end
    


    if(flag_plot)
        hold on;
        switch mode
            case 1 
                plot_x = log10(delay_weight)+J(l,1);
                plot_y = J(l,2);
            case 2
                plot_x = log10(delay_weight)+J(l,1);
                plot_y = J(l,5) ;
            case 3
                plot_x = J(l,5);
                plot_y = J(l,2);
        end
        plot(plot_x,plot_y,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',15);
        drawnow;
        title(vtSelect + " " + vdSelect);
    end

    disp_rt=['Rt= ',num2str(Rt)];
    disp(disp_rt);
    %%%%%%%%%%%%%%%%%%%Update derivation step
    dv=dvi*(Rt/Rti)^(2/3);
    dwu=dwui*(Rt/Rti)^(2/3);
    dwl=dwli*(Rt/Rti)^(2/3);
    %dwu = 0;
    %dwl = 0;
    
    %%%%%%%delay and dynamic power partial calculation
    if(flag_wire)
        [R,C0]=interconnect_model( devicemodel_allocation, area_model,Xc,optimizing_stage,optimizing_stage_weight,input_total_area,Ri,Ci);
    else
        C0=Ci;
        R=Ri;
    end

    for i_stage=1:n_bins
        fprintf ('i_stage: %i\n', i_stage);
        n=2*optimizing_stage(i_stage);
        n_partial=1+n+flag_vdd+flag_vtn+flag_vtp+(n+flag_vdd+flag_vtn+flag_vtp+1)*(n+flag_vdd+flag_vtn+flag_vtp)/2;
        f2=0;
        delay_and_power_partial_simulation;
        f3=0;
        leakage_calculation;
    
        if(f2==1)
            disp('faild in getting dynamic or delay partial');
            break;
        end
        if(f3==1)
            disp('faild in getting leakage partial');
            break;
        end
        Derivative(i_stage).D=D;
        Derivative(i_stage).DE=DE;
        Derivative(i_stage).LE=LE;
    end 

    if(f2||f3)
        continue;
    end
    
    fprintf('its here before derivative, f2 = %i - f3 = %i\n',f2,f3);

    fde0=0;
    [current_area, G_A] = chip_area_model(devicemodel_allocation, area_model, Xc,optimizing_stage,optimizing_stage_weight );
    for i_stage=1:n_bins
        n=2*optimizing_stage(i_stage);
        n_partial=1+n+flag_vdd+flag_vtn+flag_vtp+(n+flag_vdd+flag_vtn+flag_vtp+1)*(n+flag_vdd+flag_vtn+flag_vtp)/2;
        DE=[];
        LE=[];
        D=[];
        Derivative(i_stage).area = current_area/area_scale;
        DE=Derivative(i_stage).DE;
        D=Derivative(i_stage).D;
        LE=Derivative(i_stage).LE;
        derivation_calculation;
        Derivative(i_stage).D=D;
        Derivative(i_stage).E=E;
        Derivative(i_stage).d_D=d_D;
        Derivative(i_stage).d_D2=d_D2;
        Derivative(i_stage).d_E=d_E;
        Derivative(i_stage).d_E2=d_E2;
        Derivative(i_stage).d_A = G_A(:,i_stage)/area_scale;
    end
    fprintf('its here finished derivative\n');
    
    if(fde0==1)
        continue;
    end
    k=k+1;  
    model_and_descent;
    if(flag_plot)
        hold on;
        switch mode
            case 1 
                plot_x = log10(delay_weight)+Jk(:,1);
                plot_y = Jk(:,2);
            case 2
                plot_x = log10(delay_weight)+Jk(:,1);
                plot_y = Jk(:,5) ;
            case 3
                plot_x = Jk(:,5);
                plot_y = Jk(:,2);
        end
        plot(plot_x ,plot_y,'rs',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10);
        drawnow;
    end


    %%%% test descent results
    cout_max_metric_vio =0;
    for ij=1:nk
        fprintf('Checking new points (%i/%i) ...\n', ij, nk);
        if(flag_plot)
            hold on;
            switch mode
            case 1 
                plot_x = log10(delay_weight)+Jk(ij,1);
                plot_y = Jk(ij,2);
            case 2
                plot_x = log10(delay_weight)+Jk(ij,1);
                plot_y = Jk(ij,5) ;
            case 3
                plot_x = Jk(ij,5);
                plot_y = Jk(ij,2);
            end
            plot(plot_x,plot_y,'rs',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','m',...
                'MarkerSize',10);
            drawnow;
        end
        f4=0;
        f2=0;
        if(flag_wire)
            [R,C0]=interconnect_model( devicemodel_allocation, area_model,Xkt(ij).X,optimizing_stage,optimizing_stage_weight,input_total_area,Ri,Ci);
        else
            C0=Ci;
            R=Ri;
        end
        max_D=0;

        Xkc = Xkt(ij).X;
        width_mask = Xkc(3+1:3+2*n_bins, 1:n_bins) ~= 0;
        % new number of fins
        new_nf = round(Xkt(ij).X(3+1:3+2*n_bins, 1:n_bins) / width_eff_mult);
        Xkt(ij).X(3+1:3+2*n_bins, 1:n_bins) = max(2, new_nf) * (width_eff_mult) .* width_mask;
        Xkc = Xkt(ij).X;
        % xk_initial_delay_estimation_ng_fast;
        % max_D = max(Delay1);
        

        if(f2==1)
            continue;
        end

        Delay1 = 2 * max(J_DB(l,1:n_bins)); % Delay1 is the 10x delay of the current point (Then you don't need to do initial simulation!!! Save a lot of time.)

        current_E=0;
        current_L=0;
        Xkc = Xkt(ij).X;
        detail_delay_and_power_simulation_ng_bndry_fast;
        if(flag_IR_drop)
            fprintf('Finish the initial simulation for IR drop power estimation\n');
            J_DB(nof(k-1)+ij,:)=Delay2;
            current_E=sum(optimizing_stage_weight.*(Energy2));      % Here is still dynamic energy
            current_L=sum(optimizing_stage_weight.*(-Leakage2));
            fprintf ('Xk Detail Sim Done - Energy2:%d, Leakage2:%d\n',Energy2,Leakage2);
            max_D=max(J_DB(nof(k-1)+ij,:));
            current_E=current_E/(delay_weight*max_D)+current_L;     % (Then you don't need to do initial simulation!!! Save a lot of time.)
            Power_total_now = current_E;
    
            detail_delay_and_power_simulation_ng_bndry_fast_IRdrop;
            fprintf('Finish the simulation for IR drop\n');
        end

        if(f4==1 )
            disp('error in Xk');
            J(nof(k-1)+ij,1)=38;
            J(nof(k-1)+ij,2)=38;
            J(nof(k-1)+ij,3)=nof(k-1)+ij;
            J(nof(k-1)+ij,4)=38;
            J(nof(k-1)+ij,5)=38;
            DB(nof(k-1)+ij).X=Xkt(ij).X;
            J_DB(nof(k-1)+ij,1:n_bins)=38;
        end
        J_DB(nof(k-1)+ij,:)=Delay2;
        current_E=sum(optimizing_stage_weight.*(Energy2));      % Here is still dynamic energy
        current_L=sum(optimizing_stage_weight.*(-Leakage2));
        fprintf ('Xk Detail Sim Done - Energy2:%d, Leakage2:%d\n',Energy2,Leakage2);


        % fprintf ('FINAL currentE:%d - currentL:%d\n',current_E,current_L);
        if(f4==1)
            continue;
        end
        max_D=max(J_DB(nof(k-1)+ij,:));


        current_E=current_E/(delay_weight*max_D)+current_L;     % (Then you don't need to do initial simulation!!! Save a lot of time.)

        fprintf ('FINAL currentE:%d - currentL:%d\n\n',current_E,current_L);
        if(f4==1)  
            continue;
        end
        %%%Check max metric violation
        [Current_area, temp_G_A] = chip_area_model(devicemodel_allocation, area_model, Xkt(ij).X, optimizing_stage,optimizing_stage_weight );
        f5 =0 ;
        switch mode
            case 1
                if( Current_area /area_scale  > max_limit_metric )
                    f5=1;
                end
            case 2
                if( current_E > max_limit_metric )
                    f5=1;
                end
            case 3
                if( max(J_DB(nof(k-1)+ij,:)) > max_limit_metric )
                    f5=1;
                end
        end
        if (f5 ==1)
            cout_max_metric_vio = cout_max_metric_vio + 1;
            switch mode
                case 1
                    message=['violate max limit metric ', num2str(Current_area ),' > ', num2str(max_limit_metric * area_scale)];
                case 2
                    message=['violate max limit metric ', num2str(current_E ),' > ', num2str(max_limit_metric )];
                case 3
                    message=['violate max limit metric ', num2str( max(J_DB(nof(k-1)+ij,:)) * delay_weight ),' > ', num2str(max_limit_metric *delay_weight)];
            end
            disp(message);
            J(nof(k-1)+ij,1)=38;
            J(nof(k-1)+ij,2)=38;
            J(nof(k-1)+ij,3)=nof(k-1)+ij;
            J(nof(k-1)+ij,4)=38;
            J(nof(k-1)+ij,5)=38;
            DB(nof(k-1)+ij).X=Xkt(ij).X;
            J_DB(nof(k-1)+ij,1:n_bins)=38;
            continue;
        end
        %%%Pass all checks, now write valid data
        J(nof(k-1)+ij,1)=log10(max(J_DB(nof(k-1)+ij,:)));
        J(nof(k-1)+ij,2)=log10(current_E);
        J(nof(k-1)+ij,3)=nof(k-1)+ij;
        J(nof(k-1)+ij,4)=log10(current_L);
        J(nof(k-1)+ij,5) = Current_area /area_scale;
        DB(nof(k-1)+ij).X=Xkt(ij).X;
        current_time=clock;
        fprintf('current time is %i:%i\n',current_time(4),current_time(5));
        if(flag_plot)
            hold on;
            switch mode
            case 1 
                plot_x = log10(delay_weight)+J(nof(k-1)+ij,1);
                plot_y = J(nof(k-1)+ij,2);
            case 2
                plot_x = log10(delay_weight)+J(nof(k-1)+ij,1);
                plot_y = J(nof(k-1)+ij,5);
            case 3
                plot_x = J(nof(k-1)+ij,5);
                plot_y = J(nof(k-1)+ij,2);
            end
            plot(plot_x,plot_y,'rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',color_ng(mod(ij,2)+1),...
                'MarkerSize',10);
                drawnow;
                title(vtSelect + " " + vdSelect);
        end
    end
    %%%%update metric_scale 
    if( cout_max_metric_vio / nk > 0.4)
        metric_scale = metric_scale * 0.7;
    end
    if( cout_max_metric_vio / nk < 0.1)
        metric_scale = metric_scale * 1.3;
    end
    if(metric_scale <1 )
        metric_scale = 1;
    end
%     if(metric_scale >10 )
%         metric_scale = 10;
%     end
    %%%%% dead with current results
    nof(k)=nof(k-1)+nk;
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
    Jtt=sortrows(real(J),[rankc1 3]);       % sort by delay and then by id
    J=Jtt;
    f5=0;
    for i=1:1:nof(k)
        if(J(i,rankc1)==38)
            f5=1;
            break;
        end
    end
    if(f5==1)
        J(i:nof(k),:)=[];
        nof(k)=i-1;
    end

    i=2;

    while(i<=nof(k))
        if(J(i,rankc2)>=J(i-1,rankc2))
            J(i,:)=[];
            nof(k)=nof(k)-1;
            continue;
        end
        i=i+1;
    end
    DBtemp=[];
    J_DBtemp=[];
    for i=1:nof(k)
        DBtemp(i).X=DB(J(i,3)).X;                       % J(i,3) is the index of the X in DB
        J_DBtemp(i,1:n_bins)=J_DB(J(i,3),1:n_bins);     % J(i,3) is the index of the X in DB
    end

    
   

    DB=[];
    DB=DBtemp;
    J_DB=[];
    J_DB=J_DBtemp;
    J(1:nof(k),3)=1:nof(k);
    if(flag_plot)
        close;
    end
    save(matname,'J','J_DB','DB','txtname','Xtxtname');
    fclose all;
    if(flag_plot)
        switch mode
            case 1 
                plot_x = log10(delay_weight)+J(:,1);
                plot_y = J(:,2);
            case 2
                plot_x = log10(delay_weight)+J(:,1);
                plot_y = J(:,5) ;
            case 3
                plot_x = J(:,5);
                plot_y = J(:,2);
        end
        h=plot(plot_x,plot_y,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
        drawnow;
        title(vtSelect + " " + vdSelect);
    end
end

for l=1:1:length(J(:,1))
    % Reselect_allocation_boundry in the end for every existing point
    if (strcmp(vdSelect, '1vd') == 1 && strcmp(vtSelect, '1vt') == 1)
        fprintf('1vd 1vt, no need to reselect the allocation boundry\n');
    else
        Xc=DB(l).X;
        reselect_allocation_boundry_vd_vt_fast;
    end
end
save(matname,'J','J_DB','DB','txtname','Xtxtname');
fclose all;

% saving all the results in txt files --> for flexible postprocessing
fprintf('Saving TXT files...\n');
if(k>1)
    fid=fopen(txtname,'w');
    for j=1:1:length(J(:,1))
        fprintf(fid,'%d %d %d %d\n',log10(delay_weight)+J(j,1),J(j,2),J(j,4),J(j,5)*area_scale);
    end
    fclose(fid);
    stage_delay=['stage_delay_',txtname];
    dlmwrite(stage_delay,'J_DB');
    fclose all;

    Xtxtname=['X_',txtname];
    % 1) *width* in the xtxt file was saved in number of fins instead of width. 0) save width instead
    savefin = 1;
    fid=fopen(Xtxtname,'w');
    for j=1:1:length(J(:,1))
        for i_stage=1:n_bins
            for i=1:2*optimizing_stage(i_stage)+3
                if (i > 3 && savefin == 1)
                    fprintf(fid,'%d ',ceil(DB(j).X(i,i_stage)/width_eff_mult));
                else
                    fprintf(fid,'%d ',DB(j).X(i,i_stage));
                end
                % DB(j).X(i,i_stage)
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end
 
 
        end
    end
end
if(mod(i_sim,2) == 1)
    next_X = DB(2);
else
    next_X = DB(length(DB)-1);
end

fprintf('Voila!s\n');
end %% end sim
    
% end_t=cputime;
% runtime=['Run time is ',num2str(end_t-start_t),'\n'];
% fprintf(runtime);
toc;
% end
%data_deal;
