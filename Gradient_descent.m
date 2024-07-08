function [Xf,num_points_return,Jf]=Gradient_descent(w_1, w_2, Xc, Rt, trust_int_v, dwui, dwli, Derivative, num, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric, metric_scale)

%%%%%Debug testing
% Xc=[1 2 3 4 0 0;1 3 4 5 6 7]';
% w_1=1/2;
% w_2=1/2;
% ub_i=[2 3 4 5 0 0; 2 4 5 6 7 8]';
% lb_i=[0 1 2 3 0 0; 0 2 3 4 5 6]';
% Derivative(1).d_D=[2;1;3;4];
% Derivative(2).d_D=[1;1;2;2;2;2];
% Derivative(1).d_D2=[1 2 3 4;2 1 3 4;3 3 1 1; 4 4 1 1];
% Derivative(2).d_D2=ones(6,6);
% Derivative(1).d_E=[2;1;3;4];
% Derivative(2).d_E=[1;1;2;2;2;2];
% Derivative(1).d_E2=[1 2 3 4;2 1 3 4;3 3 1 1; 4 4 1 1];
% Derivative(2).d_E2=ones(6,6);
% Derivative(1).D(1)=1;
% Derivative(1).E(1)=1;
% Derivative(2).D(1)=2;
% Derivative(2).E(1)=2;
% optimizing_stage_weight=[1 2];
% vdd_allocation=[1 2];
% vt_allocation=[1;2];
% num=[2 4];
% flag_vdd=1;
% flag_vt=1;
%%%%%%%%%%%%% Set ub_i and lb_i to ensure gradient descent in model region
% fin pitch
width_eff_mult   = 0.075;

num_points=5;
num_points_return=2*num_points;
n_norm=200;
plus_delay=12;
a=0.1;
b=0.75;
M=0;
n_bins=length(num);
Df=-0.5;
n=n_norm;
for i_stage=1:n_bins
    D(num(i_stage)/2)=Derivative(i_stage).D(1);
    E(num(i_stage)/2)=Derivative(i_stage).E(1);
    W(num(i_stage)/2)=optimizing_stage_weight(i_stage); % the number of paths in each bins. - Zhichao
end
%%%%%%%%%%%%%%%%%%%Gradient and Hessian calculation
G_D=[];
G_P=[];
H_D=[];
H_P=[];
n_v=0;
delay_base=0;
power_base=0;
for i_stage=1:n_bins   %%% delay and power base are formular used in G and H calculations
    current_D=Derivative(i_stage).D(1);
    current_E=Derivative(i_stage).E(1);
    delay_base=delay_base+(plus_delay+current_D)^n_norm;        % plus delay is for changing the delay unit to pico second I guess? - Zhichao
    power_base=power_base+optimizing_stage_weight(i_stage)*10^(current_E);  % the sum of the total energy of all paths in the bin. - Zhichao
end
n_vdd=0;
n_vtn=0;
n_vtp=0;

%%%%%%Calculate Gradient vector; 
%%%%%Variable orders:
%%%%%vdd1,vdd2,...,vddn,vtn1,vtn2,...,vtnn,vtp1,vtp2,...,vtpn,w11,w12,..w1n,w21,..,wnn
if(flag_vdd)
    for i=1:length(vdd_allocation(:,1))
        n_v=n_v+1;
        G_D(n_v,1)=0;
        G_P(n_v,1)=0;
        G_A(n_v,1)=0;
        id_=find(num==(vdd_allocation(i,1)*2));
        Xi(n_v,1)=Xc(1,id_);
        ub(n_v,1)=ub_i(1,id_);
        lb(n_v,1)=lb_i(1,id_);
        trust_radius(n_v,1)=trust_int_v;
        for j=1:length(find(vdd_allocation(i,:)>0))
            id=vdd_allocation(i,j);
            G_D(n_v,1) = G_D(n_v,1) + delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D(1,1);
            G_P(n_v,1) = G_P(n_v,1) + optimizing_stage_weight(id)*10^(E(id))*Derivative(id).d_E(1,1)/power_base;
        end
    end
    % n_vdd=n_v;
    % if 2-vdd case, now n_v = 2
end
n_vdd=n_v;

if(flag_vtn)
    for i=1:length(vtn_allocation(:,1))
        n_v=n_v+1;
        G_D(n_v,1)=0;
        G_P(n_v,1)=0;
        G_A(n_v,1)=0;
        id_=find(num==(vtn_allocation(i,1)*2));
        Xi(n_v,1)=Xc(2,id_);
        ub(n_v,1)=ub_i(2,id_);
        lb(n_v,1)=lb_i(2,id_);
        trust_radius(n_v,1)=trust_int_v;
        for j=1:length(find(vtn_allocation(i,:)>0))
            id=vtn_allocation(i,j);
            G_D(n_v,1) = G_D(n_v,1) + delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D(2,1);
            G_P(n_v,1) = G_P(n_v,1) + optimizing_stage_weight(id)*10^(E(id))*Derivative(id).d_E(2,1)/power_base;
        end
    end
    % n_vtn=n_v;
end
% if 2-vtn case, now n_v = 4
n_vtn=n_v;

if(flag_vtp)
    for i=1:length(vtp_allocation(:,1))
        n_v=n_v+1;
        G_D(n_v,1)=0;
        G_P(n_v,1)=0;
        G_A(n_v,1)=0;
        id_=find(num==(vtp_allocation(i,1)*2));
        Xi(n_v,1)=Xc(3,id_);
        ub(n_v,1)=ub_i(3,id_);
        lb(n_v,1)=lb_i(3,id_);
        trust_radius(n_v,1)=trust_int_v;
        for j=1:length(find(vtp_allocation(i,:)>0))
            id=vtp_allocation(i,j);
            G_D(n_v,1) = G_D(n_v,1) + delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D(3,1);
            G_P(n_v,1) = G_P(n_v,1) + optimizing_stage_weight(id)*10^(E(id))*Derivative(id).d_E(3,1)/power_base;
        end
    end
    % n_vtp=n_v;
end
% if 2-vtp case, now n_v = 6
n_vtp=n_v;

%%%
for i_stage=1:n_bins
    for i=1:num(i_stage)    % i = 1:2 or 1:4 or 1:6 or 1:8 or 1:10
        id=num(i_stage)/2;      % num = 2*optimizing_stage.     Hence, id = 1 or 2 or 3 or 4 or 5, 
        n_v=n_v+1;
        Xi(n_v,1)=Xc(i+3,i_stage);
        ub(n_v,1)=ub_i(i+3,i_stage);
        lb(n_v,1)=lb_i(i+3,i_stage);
        nf1 = max(2, ceil((1-dwli)*Xc(i+3,i_stage)/width_eff_mult));   % Here different from the file 'write_derivation'. I remove the flag_w.
        nf2 = ceil((1+dwui)*Xc(i+3,i_stage)/width_eff_mult);   % Or when flag_w=1, then nf1=nf2, the derivation computation will have an error of dividing by zero
        trust_radius(n_v,1) = nf2 - nf1;  
        % trust_radius(n_v,1)=trust_int_w;        % Here trust_int_w = dwui + dwli
        G_D(n_v,1)=delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D(3+i,1);
        G_P(n_v,1)=optimizing_stage_weight(id)*10^(E(id))*Derivative(id).d_E(3+i,1)/power_base;
        % fprintf('n_v: %d\n',n_v);
        % G_P
        G_A(n_v,1)= Derivative(id).d_A(3+i,1);
        % G_A
    end
end
% fprintf('i_stage: %i\n',num(i_stage));
% G_D
% G_P
% G_A
% fprintf('UP CASE n_vdd: %i - n_vtn: %i - n_vtp: %i\n',n_vdd,n_vtn,n_vtp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Calculate Hessian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flag_vdd)
    for i=1:length(vdd_allocation(:,1))
        H_D(i,i)=0;
        H_P(i,i)=0;
        for j=1:length(find(vdd_allocation(i,:)>0))%%%%%d2D/dvdd2
            id=vdd_allocation(i,j);
            for j2=1:length(find(vdd_allocation(i,:)>0))
                id2=vdd_allocation(i,j2);       % id and id2: 1,1  1,2  1,3  2,1  2,2  2,3 ...
                H_D(i,i)=H_D(i,i)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                *Derivative(id).d_D(1,1)*Derivative(id2).d_D(1,1);
                H_P(i,i)=H_P(i,i)-log(10)*W(id2)*10^E(id2)*W(id)*10^E(id)*Derivative(id2).d_E(1,1)...
                    *Derivative(id).d_E(1,1)/power_base^2;
            end
            H_D(i,i)=H_D(i,i)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-2)*Derivative(id).d_D(1,1)^2 ...
                   +delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D2(1,1);
            H_P(i,i)=H_P(i,i)+log(10)*W(id)*10^E(id)*(Derivative(id).d_E(1,1))^2/power_base+W(id)*10^E(id)*Derivative(id).d_E2(1,1)/power_base;
        end

        % fprintf('before vdd alloc\n');
        % H_D
        for j=i+1:length(vdd_allocation(:,1)) %%%%d2D/(dvdd1 dvdd2)
            H_D(i,j)=0;
            H_P(i,j)=0;
            for j1=1:length(find(vdd_allocation(i,:)>0))
                id1=vdd_allocation(i,j1);   % for vdd1
                for j2=1:length(find(vdd_allocation(j,:)>0))
                    id2=vdd_allocation(j,j2);   % for vdd2
                    H_D(i,j)=H_D(i,j)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(1,1)*Derivative(id2).d_D(1,1);
                    H_P(i,j)=H_P(i,j)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(1,1)...
                    *Derivative(id1).d_E(1,1)/power_base^2;
                end
            end
        end
        % fprintf('before vtn alloc\n');
        % H_D
        for j=1:length(vtn_allocation(:,1))%%%%d2D/dvdd1 dvtn1
            H_D(i,n_vdd+j)=0;
            H_P(i,n_vdd+j)=0;
            for j1=1:length(find(vdd_allocation(i,:)>0))    % j1 is for Vdd1/Vdd2
                for j2=1:length(find(vtn_allocation(j,:)>0))    % j2 is for Vtn1/Vtn2
                    id1=vdd_allocation(i,j1);
                    id2=vtn_allocation(j,j2);
                    H_D(i,n_vdd+j)=H_D(i,n_vdd+j)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(1,1)*Derivative(id2).d_D(2,1);
                    H_P(i,n_vdd+j)=H_P(i,n_vdd+j)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(2,1)...
                    *Derivative(id1).d_E(1,1)/power_base^2;
                    if(id1==id2)
                        H_D(i,n_vdd+j)=H_D(i,n_vdd+j)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(1,1)...
                            *Derivative(id1).d_D(2,1)+delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(1,2);
                        H_P(i,n_vdd+j)=H_P(i,n_vdd+j)+log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(1,1)*Derivative(id1).d_E(2,1)...
                            /power_base+W(id1)*10^E(id1)*Derivative(id1).d_E2(1,2)/power_base;
                    end
                end
            end
        end
        % fprintf('before vtp alloc\n');
        % H_D
        for j=1:length(vtp_allocation(:,1))%%%%d2D/dvdd1 dvtp1
            H_D(i,n_vtn+j)=0;
            H_P(i,n_vtn+j)=0;
            for j1=1:length(find(vdd_allocation(i,:)>0))
                for j2=1:length(find(vtp_allocation(j,:)>0))
                    id1=vdd_allocation(i,j1);
                    id2=vtp_allocation(j,j2);
                    H_D(i,n_vtn+j)=H_D(i,n_vtn+j)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(1,1)*Derivative(id2).d_D(3,1);
                    H_P(i,n_vtn+j)=H_P(i,n_vtn+j)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(3,1)...
                    *Derivative(id1).d_E(1,1)/power_base^2;
                    if(id1==id2)
                        H_D(i,n_vtn+j)=H_D(i,n_vtn+j)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(1,1)...
                            *Derivative(id1).d_D(3,1)+delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(1,3);
                        H_P(i,n_vtn+j)=H_P(i,n_vtn+j)+log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(1,1)*Derivative(id1).d_E(3,1)...
                            /power_base+W(id1)*10^E(id1)*Derivative(id1).d_E2(1,3)/power_base;
                    end
                end
            end
        end

        % fprintf('before nVTP - d2D/dvdd1 dw1\n');
        % H_D
        c_i=n_vtp;
        for i_stage=1:n_bins %%% d2D/dvdd1 dw1
            id2=num(i_stage)/2;
            for j=1:num(i_stage)
                c_i=c_i+1;
                % fprintf('c_i up: %i\n',c_i);
                H_D(i,c_i)=0;
                H_P(i,c_i)=0;
                for j2=1:length(find(vdd_allocation(i,:)>0))
                    id1=vdd_allocation(i,j2);
                    H_D(i,c_i)=H_D(i,c_i)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(1,1)*Derivative(id2).d_D(3+j,1);                    % d_D(3+j,1) is d_D(4,1) -> d_D(13,1)
                    H_P(i,c_i)=H_P(i,c_i)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(3+j,1)...
                    *Derivative(id1).d_E(1,1)/power_base^2;
                    if(id1==id2)
                        H_D(i,c_i)=H_D(i,c_i)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(1,1)...
                            *Derivative(id1).d_D(3+j,1)+delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(1,3+j);                  
                        H_P(i,c_i)=H_P(i,c_i)+log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(1,1)*Derivative(id1).d_E(3+j,1) /power_base+W(id1)*10^E(id1)*Derivative(id1).d_E2(1,3+j)/power_base;
                    end
                end
            end
        end
        % fprintf('after nVTP - d2D/dvdd1 dw1\n');
        % H_D
    end
end

% fprintf('After flag_vdd - n_vtn val: %i - n_vtp val: %i\n',n_vtn,n_vtp);
% H_D
% pause(1000000);

if(flag_vtn)
    for i=1:length(vtn_allocation(:,1))
        H_D(n_vdd+i,n_vdd+i)=0;
        H_P(n_vdd+i,n_vdd+i)=0;
        for j=1:length(find(vtn_allocation(i,:)>0))%%%%%d2D/dvtn^2
            id=vtn_allocation(i,j);
            for j2=1:length(find(vtn_allocation(i,:)>0))
                id2=vtn_allocation(i,j2);
                H_D(n_vdd+i,n_vdd+i)=H_D(n_vdd+i,n_vdd+i)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id))^(n_norm-1)*(plus_delay...
                    +D(id2))^(n_norm-1)*Derivative(id).d_D(2,1)*Derivative(id2).d_D(2,1);
                H_P(n_vdd+i,n_vdd+i)=H_P(n_vdd+i,n_vdd+i)-log(10)*W(id2)*10^E(id2)*W(id)*10^E(id)*Derivative(id2).d_E(2,1)...
                    *Derivative(id).d_E(2,1)/power_base^2;
            end
            H_D(n_vdd+i,n_vdd+i)=H_D(n_vdd+i,n_vdd+i)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-2)*Derivative(id).d_D(2,1)^2 ...
                +delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D2(2,2);
            H_P(n_vdd+i,n_vdd+i)=H_P(n_vdd+i,n_vdd+i)+log(10)*W(id)*10^E(id)*(Derivative(id).d_E(2,1))^2/power_base+W(id)*10^E(id)*Derivative(id).d_E2(2,2)/power_base;
        end
        
        % fprintf('(VTN) before vtn alloc\n');
        % H_D
        for j=i+1:length(vtn_allocation(:,1)) %%%%d2D/(dvtn1 dvtn2)
            H_D(n_vdd+i,n_vdd+j)=0;
            H_P(n_vdd+i,n_vdd+j)=0;
            for j1=1:length(find(vtn_allocation(i,:)>0))
                id1=vtn_allocation(i,j1);
                for j2=1:length(find(vtn_allocation(j,:)>0))
                    id2=vtn_allocation(j,j2);
                    H_D(n_vdd+i,n_vdd+j)=H_D(n_vdd+i,n_vdd+j)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(2,1)*Derivative(id2).d_D(2,1);
                    H_P(n_vdd+i,n_vdd+j)=H_P(n_vdd+i,n_vdd+j)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(2,1)...
                    *Derivative(id1).d_E(2,1)/power_base^2;
                end
            end
        end
        
        % new!!!! Added by Zhichao
        % H_D
        for j=1:length(vtp_allocation(:,1)) %%%%d2D/(dvtn1 dvtp1)
            H_D(n_vdd+i,n_vtn+j)=0;
            H_P(n_vdd+i,n_vtn+j)=0;
            for j1=1:length(find(vtn_allocation(i,:)>0))        % for vtn1/vtn2
                for j2=1:length(find(vtp_allocation(j,:)>0))    % for vtp1/vtp2
                    id1=vtn_allocation(i,j1);
                    id2=vtp_allocation(j,j2);
                    H_D(n_vdd+i,n_vtn+j)=H_D(n_vdd+i,n_vtn+j)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(2,1)*Derivative(id2).d_D(3,1);
                    H_P(n_vdd+i,n_vtn+j)=H_P(n_vdd+i,n_vtn+j)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(3,1)...
                    *Derivative(id1).d_E(2,1)/power_base^2;
                    if(id1==id2)
                        H_D(n_vdd+i,n_vtn+j)=H_D(n_vdd+i,n_vtn+j)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(2,1)...
                            *Derivative(id1).d_D(3,1)+delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(2,3);
                        H_P(n_vdd+i,n_vtn+j)=H_P(n_vdd+i,n_vtn+j)+log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(2,1)*Derivative(id1).d_E(3,1)...
                            /power_base+W(id1)*10^E(id1)*Derivative(id1).d_E2(2,3)/power_base;
                    end
                end
            end
        end

        c_i=n_vtp;
        % fprintf('(VTN) after VTN alloc - c_i: %i\n',c_i);
        % H_D
        for i_stage=1:n_bins %%% d2D/dvtn1 dw1
            id2=num(i_stage)/2;
            for j=1:num(i_stage)        % for w1, w2, w3, ..., w10
                c_i=c_i+1;
                H_D(n_vdd+i,c_i)=0;
                H_P(n_vdd+i,c_i)=0;
                for j2=1:length(find(vtn_allocation(i,:)>0))        % for vtn1/vtn2
                    id1=vtn_allocation(i,j2);
                    H_D(n_vdd+i,c_i)=H_D(n_vdd+i,c_i)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(2,1)*Derivative(id2).d_D(3+j,1);
                    H_P(n_vdd+i,c_i)=H_P(n_vdd+i,c_i)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(3+j,1)...
                    *Derivative(id1).d_E(2,1)/power_base^2;
                    if(id1==id2)
                        H_D(n_vdd+i,c_i)=H_D(n_vdd+i,c_i)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(2,1)...
                            *Derivative(id1).d_D(3+j,1)+delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(2,3+j);                  
                        H_P(n_vdd+i,c_i)=H_P(n_vdd+i,c_i)+log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(2,1)*Derivative(id1).d_E(3+j,1)/power_base+W(id1)*10^E(id1)*Derivative(id1).d_E2(2,3+j)/power_base;

                    end
                end
            end
        end
        % fprintf('(VTN) after d2D/dvdd1 dw1 END - c_i: %i\n',c_i);
        % H_D
    end
end 

% H_D
% fprintf('pause at grad descent');
% pause(1000000);

if(flag_vtp)
    for i=1:length(vtp_allocation(:,1))
        H_D(n_vtn+i,n_vtn+i)=0;
        H_P(n_vtn+i,n_vtn+i)=0;
        for j=1:length(find(vtp_allocation(i,:)>0))%%%%%d^2D/dvtp^2
            id=vtp_allocation(i,j);
            for j2=1:length(find(vtp_allocation(i,:)>0))
                id2=vtp_allocation(i,j2);
                H_D(n_vtn+i,n_vtn+i)=H_D(n_vtn+i,n_vtn+i)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id))^(n_norm-1)*(plus_delay...
                    +D(id2))^(n_norm-1)*Derivative(id).d_D(3,1)*Derivative(id2).d_D(3,1);
                H_P(n_vtn+i,n_vtn+i)=H_P(n_vtn+i,n_vtn+i)-log(10)*W(id2)*10^E(id2)*W(id)*10^E(id)*Derivative(id2).d_E(3,1)...
                    *Derivative(id).d_E(3,1)/power_base^2;
            end
            H_D(n_vtn+i,n_vtn+i)=H_D(n_vtn+i,n_vtn+i)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-2)*Derivative(id).d_D(3,1)^2 ...
                +delay_base^(1/n_norm-1)*(plus_delay+D(id))^(n_norm-1)*Derivative(id).d_D2(3,3);
            H_P(n_vtn+i,n_vtn+i)=H_P(n_vtn+i,n_vtn+i)+log(10)*W(id)*10^E(id)*(Derivative(id).d_E(3,1))^2/power_base+W(id)*10^E(id)*Derivative(id).d_E2(3,3)/power_base;
        end
        % fprintf('(VTP) before vtp alloc\n');
        % H_D
        for j=i+1:length(vtp_allocation(:,1)) %%%%d2D/(dvtp1 dvtp2)
            H_D(n_vtn+i,n_vtn+j)=0;
            H_P(n_vtn+i,n_vtn+j)=0;
            for j1=1:length(find(vtp_allocation(i,:)>0))
                id1=vtp_allocation(i,j1);
                for j2=1:length(find(vtp_allocation(j,:)>0))
                    id2=vtp_allocation(j,j2);
                    H_D(n_vtn+i,n_vtn+j)=H_D(n_vtn+i,n_vtn+j)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(3,1)*Derivative(id2).d_D(3,1);
                    H_P(n_vtn+i,n_vtn+j)=H_P(n_vtn+i,n_vtn+j)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(3,1)...
                    *Derivative(id1).d_E(3,1)/power_base^2;
                end
            end
        end

        c_i=n_vtp;
        % fprintf('(VTP) after VTP alloc - c_i: %i\n',c_i);
        % H_D
        for i_stage=1:n_bins %%% d2D/dvdd1 dw1
            id2=num(i_stage)/2;
            for j=1:num(i_stage)
                c_i=c_i+1;
                H_D(n_vtn+i,c_i)=0;
                H_P(n_vtn+i,c_i)=0;
                for j2=1:length(find(vtp_allocation(i,:)>0))
                    id1=vtp_allocation(i,j2);
                    H_D(n_vtn+i,c_i)=H_D(n_vtn+i,c_i)+(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)...
                        *Derivative(id1).d_D(3,1)*Derivative(id2).d_D(3+j,1);
                    H_P(n_vtn+i,c_i)=H_P(n_vtn+i,c_i)-log(10)*W(id2)*10^E(id2)*W(id1)*10^E(id1)*Derivative(id2).d_E(3+j,1)...
                    *Derivative(id1).d_E(3,1)/power_base^2;
                    if(id1==id2)
                        H_D(n_vtn+i,c_i)=H_D(n_vtn+i,c_i)+(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(3,1)...
                            *Derivative(id1).d_D(3+j,1)+delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(3,3+j);                  
                        H_P(n_vtn+i,c_i)=H_P(n_vtn+i,c_i)+log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(3,1)*Derivative(id1).d_E(3+j,1)/power_base+W(id1)*10^E(id1)*Derivative(id1).d_E2(3,3+j)/power_base;

                    end
                end
            end
        end
        % fprintf('(VTP) after d2D/dvdd1 dw1 END - c_i: %i\n',c_i);
        % H_D
    end
end 

% fprintf('After flag_vtp - n_vtp val: %i\n',n_vtp);
% H_D

c_i=n_vtp;
for i_stage=1:n_bins
    id1=num(i_stage)/2;
    for j1=1:num(i_stage)
        c_i=c_i+1;
        for j2=j1:num(i_stage)
            H_D(c_i,c_i+j2-j1)=(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(2*(n_norm-1))*Derivative(id1).d_D(3+j1,1)*Derivative(id1).d_D(3+j2,1)...
                +(n_norm-1)*delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-2)*Derivative(id1).d_D(3+j1,1)*Derivative(id1).d_D(3+j2,1)...
                +delay_base^(1/n_norm-1)*(plus_delay+D(id1))^(n_norm-1)*Derivative(id1).d_D2(3+j1,3+j2);
            H_P(c_i,c_i+j2-j1)=-log(10)*(W(id1)*10^E(id1))^2*Derivative(id1).d_E(3+j1)*Derivative(id1).d_E(3+j2,1)/power_base^2 ...
                +log(10)*W(id1)*10^E(id1)*Derivative(id1).d_E(3+j1,1)*Derivative(id1).d_E(3+j2,1)/power_base...
                +W(id1)*10^E(id1)*Derivative(id1).d_E2(3+j1,3+j2)/power_base;
        end
        c_j=c_i+num(i_stage)-j1;
        for i2=i_stage+1:n_bins
            id2=num(i2)/2;
            for j2=1:num(i2)
                c_j=c_j+1;
                H_D(c_i,c_j)=(1-n_norm)*delay_base^(1/n_norm-2)*(plus_delay+D(id1))^(n_norm-1)*(plus_delay+D(id2))^(n_norm-1)*Derivative(id1).d_D(3+j1,1)*Derivative(id2).d_D(3+j2,1);
                H_P(c_i,c_j)=-log(10)*W(id1)*10^E(id1)*W(id2)*10^E(id2)*Derivative(id1).d_E(3+j1,1)*Derivative(id2).d_E(3+j2,1)/power_base^2;
            end
        end
    end
end

% fprintf('After BIN WITH VTP THINGY (before transf) - c_i: %i\n',c_i);
% H_D
% pause(235983468)

for i=2:length(Xi)      % Modified by CZC, 'for i=2:length(Xi)' instead of 'for i=3:length(Xi)'
    for j=1:i-1
        H_D(i,j)=H_D(j,i);
        H_P(i,j)=H_P(j,i);
    end
end

% fprintf('After transf\n');
% H_D
% pause(12590235890);
% fprintf('printing Xi\n');
% Xi
switch mode
    case 1
        f = w_1*(delay_base^(1/n_norm)-plus_delay)+w_2*log10(power_base);
    case 2
        f = w_1*(delay_base^(1/n_norm)-plus_delay)+w_2*Derivative(1).area;
    case 3
        f = w_1*Derivative(1).area + w_2 * log10(power_base);
end
switch mode
    case 1
        G=G_D*w_1+G_P*w_2;
        H=H_D*w_1+H_P*w_2;
    case 2
        G=G_D*w_1+G_A*w_2;
        H=H_D*w_1;
    case 3
        G=G_A*w_1+G_P*w_2;
        H=H_P*w_2;
end
X1=Xi;
X2=Xi;
f1=f;
f2=f;
trust_ub=zeros(size(Xi));
trust_lb=zeros(size(Xi));

for i_rt=1:num_points
    trust_ub(1:n_vtn,1)=Xi(1:n_vtn,1)+trust_radius(1:n_vtn,1)*Rt/num_points*i_rt;   % trust_radius = dvi
    trust_lb(1:n_vtn,1)=Xi(1:n_vtn,1)-trust_radius(1:n_vtn,1)*Rt/num_points*i_rt;   
    % fprintf('trust 1:\n');
    % trust_ub
    % trust_lb
    trust_ub(n_vtn+1:n_vtp,1)=Xi(n_vtn+1:n_vtp,1)+trust_radius(n_vtn+1:n_vtp,1)*Rt/num_points*i_rt;
    trust_lb(n_vtn+1:n_vtp,1)=Xi(n_vtn+1:n_vtp,1)-trust_radius(n_vtn+1:n_vtp,1)*Rt/num_points*i_rt;
    % fprintf('trust 2\n');
    % trust_ub
    % trust_lb
    trust_ub(n_vtp+1:length(Xi),1)=Xi(n_vtp+1:length(Xi),1).*(1+trust_radius(n_vtp+1:length(Xi),1)*Rt/num_points*i_rt);
    % trust_lb(n_vtp+1:length(Xi),1)=Xi(n_vtp+1:length(Xi),1).*(1-trust_radius(n_vtp+1:length(Xi),1)/4*Rt/num_points*i_rt);
    % Modified by czc
    trust_lb(n_vtp+1:length(Xi),1)=Xi(n_vtp+1:length(Xi),1).*(1-trust_radius(n_vtp+1:length(Xi),1)*Rt/num_points*i_rt);

    % trust_lb(n_vtp+1:length(Xi),1)=max(trust_lb(n_vtp+1:length(Xi),1),2*width_eff_mult);

    % fprintf('trust 3:\n');
    % trust_ub
    % trust_lb


    %%%%%Newton method
    G1=G+H*(X1-Xi);
    % fprintf('start\n');
    % G
    % H
    % X1
    % pause(10000);

    switch mode
        case 1
            [X1, f1]=Logarithmic_barrier_N(X1,f1,G1,H,min(ub,trust_ub),max(lb,trust_lb),G_A,zeros(size(H_D)),max_limit_metric - Derivative(1).area );
        case 2
            [X1, f1]=Logarithmic_barrier_N(X1,f1,G1,H,min(ub,trust_ub),max(lb,trust_lb),G_P, H_P , log10(max_limit_metric * metric_scale) - log10(power_base));
        case 3
            [X1, f1]=Logarithmic_barrier_N(X1,f1,G1,H,min(ub,trust_ub),max(lb,trust_lb),G_D, H_D , log10(max_limit_metric * metric_scale) - (delay_base^(1/n_norm) -plus_delay) );
    end
    % fprintf('------ after switch case ------\n');
    % X1
    % f1
    % pause(10000);
    % fprintf('Derivative 1 area: %d - ',Derivative(1).area);
    % fprintf('end log1\n');

    if(flag_vdd)
        for i_vdd=1:length(vdd_allocation(:,1))
            for i_id=1:length(find(vdd_allocation(i_vdd,:)>0))
                id=find(num==(vdd_allocation(i_vdd,i_id)*2));
                Xf(i_rt*2-1).X(1,id)=X1(i_vdd,1);
            end
        end
    else
        n_vdd=length(vdd_allocation(:,1));
        for i_vdd=1:length(vdd_allocation(:,1))
            for i_id=1:length(find(vdd_allocation(i_vdd,:)>0))
                id=find(num==(vdd_allocation(i_vdd,i_id)*2));
                Xf(i_rt*2-1).X(1,id)=Xc(1,id);
            end
        end
        n_vdd=0;
    end
    if(flag_vtn)
        for i_vtn=1:length(vtn_allocation(:,1))
            for i_id=1:length(find(vtn_allocation(i_vtn,:)>0))
                id=find(num==(vtn_allocation(i_vtn,i_id)*2));
                if(abs(X1(n_vdd+i_vtn,1))<0.001)
                    Xf(i_rt*2-1).X(2,id)=0;
                else
                    Xf(i_rt*2-1).X(2,id)=X1(n_vdd+i_vtn,1);
                end
            end
        end
    else
        for i_vtn=1:length(vtn_allocation(:,1))
            for i_id=1:length(find(vtn_allocation(i_vtn,:)>0))
                id=find(num==(vtn_allocation(i_vtn,i_id)*2));
                Xf(i_rt*2-1).X(2,id)=Xc(2,id);
            end
        end
        n_vtn=n_vdd;
    end
    % fprintf('>> printing Xf_Vtn_TOP\n');
    % fprintf('i_rt*2-1: %i - Xf1:%d\n',(i_rt*2-1),Xf(i_rt*2-1).X(2,id));

    if(flag_vtp)
        for i_vtp=1:length(vtp_allocation(:,1))
            for i_id=1:length(find(vtp_allocation(i_vtp,:)>0))
                id=find(num==(vtp_allocation(i_vtp,i_id)*2));
                if(abs(X1(n_vtn+i_vtp,1))<0.001)
                    Xf(i_rt*2-1).X(3,id)=0;
                else
                    Xf(i_rt*2-1).X(3,id)=X1(n_vtn+i_vtp,1);
                end
            end
        end
    else
        for i_vtp=1:length(vtp_allocation(:,1))
            for i_id=1:length(find(vtp_allocation(i_vtp,:)>0))
                id=find(num==(vtp_allocation(i_vtp,i_id)*2));
                Xf(i_rt*2-1).X(3,id)=Xc(3,id);
            end
        end
        n_vtp=n_vtn;
    end
    % fprintf('printing Xf_Vtp_TOP\n');
    % fprintf('i_rt*2-1: %i - Xf1:%d\n',(i_rt*2-1), Xf(i_rt*2-1).X(3,id));

    c_i=n_vtp;
    for i_stage=1:n_bins
        for i_w=1:num(i_stage)
            c_i=c_i+1;
            if(flag_w)
                % if you only optimize the width of stage 1 and copy it to other stages, use this:
                % if(i_w<=2)
                %     Xf(i_rt*2-1).X(3+i_w,i_stage)=X1(c_i,1);            % Here update the width
                %     if(i_w==1)
                %         c_i_wid_ind_nand = c_i;
                %     else    % i_w == 2
                %         c_i_wid_ind_inv = c_i;
                %     end
                % else
                %     if(mod(i_w,2)==1)
                %         Xf(i_rt*2-1).X(3+i_w,i_stage)=X1(c_i_wid_ind_nand,1);
                %     else
                %         Xf(i_rt*2-1).X(3+i_w,i_stage)=X1(c_i_wid_ind_inv,1);
                %     end
                %     % Xf(i_rt*2).X(3+i_w,i_stage)=X1(c_i_wid_ind,1);            % Here update the width
                % end
                % if you want all width optimize independently, use this:
                Xf(i_rt*2-1).X(3+i_w,i_stage)=X1(c_i,1);            % Here update the width
            else
                Xf(i_rt*2-1).X(3+i_w,i_stage)=Xc(3+i_w,i_stage);    % Here don't update the width, Modified by CZC
            end
        end
    end
    % c_i=n_vtp;
    % for i_stage=1:n_bins
    %     for i_w=1:num(i_stage)
    %         c_i=c_i+1;
    %         Xf(i_rt*2-1).X(num(i_stage)+2+i_w,i_stage)=X1(c_i,1);
    %     end
    % end
    Jf(i_rt*2-1,1)=delay_base^(1/n_norm)-plus_delay+(X1-Xi)'*G_D+1/2*(X1-Xi)'*H_D*(X1-Xi);
    Jf(i_rt*2-1,2)=log10(power_base)+(X1-Xi)'*G_P+1/2*(X1-Xi)'*H_P*(X1-Xi)-Jf(i_rt*2-1,1);
    Jf(i_rt*2-1,5)=Derivative(1).area + (X1-Xi)'*G_A;

    %%%%%%%%%%%%Gradient method
    G2=G+H*(X2-Xi);
    switch mode
        % metric is used to scale max power or delay up in mode 2's logarithmic barrier method to avoid model error
        case 1
            [X2, f2]=Logarithmic_barrier_G(X2,f2,G2,H,min(ub,trust_ub),max(lb,trust_lb), G_A,zeros(size(H_D)),max_limit_metric - Derivative(1).area);
        case 2
            [X2, f2]=Logarithmic_barrier_G(X2,f2,G2,H,min(ub,trust_ub),max(lb,trust_lb), G_P, H_P , log10(max_limit_metric * metric_scale) - log10(power_base));
        case 3
            [X2, f2]=Logarithmic_barrier_N(X2,f2,G2,H,min(ub,trust_ub),max(lb,trust_lb), G_D, H_D , log10(max_limit_metric * metric_scale) - (delay_base^(1/n_norm) -plus_delay) );
    end
    
    if(flag_vdd)
        for i_vdd=1:length(vdd_allocation(:,1))
            for i_id=1:length(find(vdd_allocation(i_vdd,:)>0))
                id=find(num==(vdd_allocation(i_vdd,i_id)*2));
                Xf(i_rt*2).X(1,id)=X2(i_vdd,1);
            end
        end
    else
        for i_vdd=1:length(vdd_allocation(:,1))
            for i_id=1:length(find(vdd_allocation(i_vdd,:)>0))
                id=find(num==(vdd_allocation(i_vdd,i_id)*2));
                Xf(i_rt*2).X(1,id)=Xc(1,id);
            end
        end
        n_vdd=0;
    end
    % fprintf('1. n_vdd: %i, n_vtn: %i, n_vtp: %i - lentgh Xi: %i\n',n_vdd,n_vtn,n_vtp,length(Xi));
    if(flag_vtn)
        for i_vtn=1:length(vtn_allocation(:,1))
            for i_id=1:length(find(vtn_allocation(i_vtn,:)>0))
                id=find(num==(vtn_allocation(i_vtn,i_id)*2));
                if(abs(X2(n_vdd+i_vtn,1))<0.001)
                    Xf(i_rt*2).X(2,id)=0;
                else
                    Xf(i_rt*2).X(2,id)=X2(n_vdd+i_vtn,1);
                end
            end
        end
    else
        for i_vtn=1:length(vtn_allocation(:,1))
            for i_id=1:length(find(vtn_allocation(i_vtn,:)>0))
                id=find(num==(vtn_allocation(i_vtn,i_id)*2));
                Xf(i_rt*2).X(2,id)=Xc(2,id);
                end
        end
        n_vtn=n_vdd;
    end
    % fprintf('2. n_vdd: %i, n_vtn: %i, n_vtp: %i - lentgh Xi: %i\n',n_vdd,n_vtn,n_vtp,length(Xi));
    % fprintf('>> printing Xf_Vtn_BOT\n');
    % fprintf('i_rt*2-1:%i - Xf2:%d\n',(i_rt*2-1),Xf(i_rt*2-1).X(2,id));
    % fprintf('i_rt*2:%i - Xf2:%d\n',(i_rt*2),Xf(i_rt*2).X(2,id));

    if(flag_vtp)
        for i_vtp=1:length(vtp_allocation(:,1))
            for i_id=1:length(find(vtp_allocation(i_vtp,:)>0))
                id=find(num==(vtp_allocation(i_vtp,i_id)*2));
                if(abs(X2(n_vtn+i_vtp,1))<0.001)
                    Xf(i_rt*2).X(3,id)=0;
                else
                    Xf(i_rt*2).X(3,id)=X2(n_vtn+i_vtp,1);
                end
            end
        end
    else
        for i_vtp=1:length(vtp_allocation(:,1))
            for i_id=1:length(find(vtp_allocation(i_vtp,:)>0))
                id=find(num==(vtp_allocation(i_vtp,i_id)*2));
                Xf(i_rt*2).X(3,id)=Xc(3,id);
            end
        end
        n_vtp=n_vtn;
    end
    % fprintf('3.n_vdd: %i, n_vtn: %i, n_vtp: %i - lentgh Xi: %i\n',n_vdd,n_vtn,n_vtp,length(Xi));
    % fprintf('printing Xf_Vtp_BOT\n');
    % fprintf('i_rt*2-1:%i - Xf2:%d\n',(i_rt*2-1),Xf(i_rt*2-1).X(3,id));
    % fprintf('i_rt*2:%i - Xf2:%d\n',(i_rt*2),Xf(i_rt*2).X(3,id));

    c_i=n_vtp;
    for i_stage=1:n_bins
        for i_w=1:num(i_stage)
            c_i=c_i+1;
            if(flag_w)
                % if you only optimize the width of stage 1 and copy it to other stages, use this:
                % if(i_w<=2)
                %     Xf(i_rt*2).X(3+i_w,i_stage)=X2(c_i,1);            % Here update the width
                %     if(i_w==1)
                %         c_i_wid_ind_nand = c_i;
                %     else    % i_w == 2
                %         c_i_wid_ind_inv = c_i;
                %     end
                % else
                %     if(mod(i_w,2)==1)
                %         Xf(i_rt*2).X(3+i_w,i_stage)=X2(c_i_wid_ind_nand,1);
                %     else
                %         Xf(i_rt*2).X(3+i_w,i_stage)=X2(c_i_wid_ind_inv,1);
                %     end
                %     % Xf(i_rt*2).X(3+i_w,i_stage)=X2(c_i_wid_ind,1);            % Here update the width
                % end
                % if you want all width optimize independently, use this:
                Xf(i_rt*2).X(3+i_w,i_stage)=X2(c_i,1);            % Here update the width
            else
                Xf(i_rt*2).X(3+i_w,i_stage)=Xc(3+i_w,i_stage);    % Here don't update the width, Modified by CZC
            end
            % fprintf('printing XfW_1\n');
            % Xf(i_rt*2).X(2+i_w,i_stage)
        end
    end
    % c_i=n_vtp;
    % for i_stage=1:n_bins
    %     for i_w=1:num(i_stage)
    %         c_i=c_i+1;
    %         Xf(i_rt*2).X(num(i_stage)+2+i_w,i_stage)=X2(c_i,1);
    %         % fprintf('printing XfW_2\n');
    %         % Xf(i_rt*2).X(num(i_stage)+2+i_w,i_stage)
    %     end
    % end
    Jf(i_rt*2,1)=delay_base^(1/n_norm)-plus_delay+(X2-Xi)'*G_D+1/2*(X2-Xi)'*H_D*(X2-Xi);
    Jf(i_rt*2,2)=log10(power_base)+(X2-Xi)'*G_P+1/2*(X2-Xi)'*H_P*(X2-Xi)-Jf(i_rt*2,1);
    Jf(i_rt*2,5)=Derivative(1).area + (X1-Xi)'*G_A;
end