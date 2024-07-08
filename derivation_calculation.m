n_non0_partial=n_partial;
% n_partial=1+n+flag_vdd+flag_vt+(n+flag_vdd+flag_vt+1)*(n+flag_vdd+flag_vt)/2;
% n_partial=1+n+2+(n+2+1)*(n+3)/2;
% n_partial=1+n+flag_vdd+flag_vtn+flag_vtp+(n+flag_vdd+flag_vtn+flag_vtp)*(n+flag_vdd+flag_vtn+flag_vtp+1)/2;
n_partial=1+n+3+(n+3+1)*(n+3)/2;

% fprintf('deriv calc NEW - n_non0_partial: %i - n_partial:%i\n',n_non0_partial,n_partial);
% D
% fprintf('ULTRA UP\n');


if(flag_vdd==0) % Here we put everything on their own place
    D(3:n_non0_partial+1)=D(2:n_non0_partial);
    DE(3:n_non0_partial+1)=DE(2:n_non0_partial);
    LE(3:n_non0_partial+1)=LE(2:n_non0_partial);
    D(2)=D(1);  % if flag_vdd == 1, D(2) should be the delay of X = (vdd+dv/2, vtn-dv/2, vtp-dv/2, w1, w2, ...)
    DE(2)=DE(1);
    LE(2)=LE(1);
    n_non0_partial=n_non0_partial+1;
     D((n+2+flag_vtn+flag_vtp)+(n+3)+1:n_non0_partial+(n+3))= D((n+2+flag_vtn+flag_vtp)+1:n_non0_partial);        % Here 'n_non0_partial+(n+3)' is the last element of D, and will be the 'n_partial' where everything is optimized if flag_vtn, flag_vtp are 1.
    DE((n+2+flag_vtn+flag_vtp)+(n+3)+1:n_non0_partial+(n+3))=DE((n+2+flag_vtn+flag_vtp)+1:n_non0_partial);
    LE((n+2+flag_vtn+flag_vtp)+(n+3)+1:n_non0_partial+(n+3))=LE((n+2+flag_vtn+flag_vtp)+1:n_non0_partial);
     D((n+2+flag_vtn+flag_vtp)+1:(n+2+flag_vtn+flag_vtp)+(n+3))= D(1);
    DE((n+2+flag_vtn+flag_vtp)+1:(n+2+flag_vtn+flag_vtp)+(n+3))=DE(1);
    LE((n+2+flag_vtn+flag_vtp)+1:(n+2+flag_vtn+flag_vtp)+(n+3))=LE(1);
    n_non0_partial=n_non0_partial+(n+3); 
end
% fprintf('chk1 n_non0_partial: %i - n_partial:%i\n',n_non0_partial,n_partial);

if(flag_vtn==0)
    % fprintf('D firstN\n');
    % D
    D(4:n_non0_partial+1)=D(3:n_non0_partial);      % 若flag_vdd=0，则在上面的条件中已经把vdd的位置给补上，所以此处没有问题
    DE(4:n_non0_partial+1)=DE(3:n_non0_partial);
    LE(4:n_non0_partial+1)=LE(3:n_non0_partial);
    D(3)=D(1);
    LE(3)=LE(1);
    DE(3)=DE(1);
    n_non0_partial=n_non0_partial+1;
    % fprintf('D mid 0N\n');
    % D
    % fprintf('FLAGVT0 - BEG1: %i - END1: %i\n',n+3+flag_vtp+3,n_non0_partial);
    if(flag_vdd==1)
         D((n+3+flag_vtp)+2+1:n_non0_partial+1)= D((n+3+flag_vtp)+2:n_non0_partial);
        DE((n+3+flag_vtp)+2+1:n_non0_partial+1)=DE((n+3+flag_vtp)+2:n_non0_partial);
        LE((n+3+flag_vtp)+2+1:n_non0_partial+1)=LE((n+3+flag_vtp)+2:n_non0_partial);
        n_non0_partial=n_non0_partial+1;
    end
    % D
     D((n+3+flag_vtp)+2)= D(2);
    DE((n+3+flag_vtp)+2)=DE(2);
    LE((n+3+flag_vtp)+2)=LE(2);
    % fprintf('D mid 2N\n');
    % D
     D((n+3+flag_vtp)+(n+2+flag_vtp)+(n+2)+1:n_non0_partial+(n+2))= D((n+3+flag_vtp)+(n+2+flag_vtp)+1:n_non0_partial);      % modified by czc
    DE((n+3+flag_vtp)+(n+2+flag_vtp)+(n+2)+1:n_non0_partial+(n+2))=DE((n+3+flag_vtp)+(n+2+flag_vtp)+1:n_non0_partial);
    LE((n+3+flag_vtp)+(n+2+flag_vtp)+(n+2)+1:n_non0_partial+(n+2))=LE((n+3+flag_vtp)+(n+2+flag_vtp)+1:n_non0_partial);
    % fprintf('D mid 3N\n');
    % D
     D((n+3+flag_vtp)+(n+2+flag_vtp)+1:(n+3+flag_vtp)+(n+2+flag_vtp)+(n+2))= D(1);      % modified by czc
    DE((n+3+flag_vtp)+(n+2+flag_vtp)+1:(n+3+flag_vtp)+(n+2+flag_vtp)+(n+2))=DE(1);
    LE((n+3+flag_vtp)+(n+2+flag_vtp)+1:(n+3+flag_vtp)+(n+2+flag_vtp)+(n+2))=LE(1);
    n_non0_partial=n_non0_partial+(n+2);
    % fprintf('D lastN\n');
    % D 
end

if(flag_vtp==0)
    % fprintf('D firstP\n');
    % D
     D(5:n_non0_partial+1)= D(4:n_non0_partial);        %  补上vtp一阶导数的空位
    DE(5:n_non0_partial+1)=DE(4:n_non0_partial);
    LE(5:n_non0_partial+1)=LE(4:n_non0_partial);
     D(4)= D(1);    
    LE(4)=LE(1);
    DE(4)=DE(1);
    n_non0_partial=n_non0_partial+1;
    % fprintf('D mid 0P\n');
    % D
    % fprintf('FLAGVTP0 - BEG1: %i - END1: %i\n',n+3+flag_vtp+3,n_non0_partial);
    if(flag_vdd==1)
        % this belongs to flag_vdd = 1
         D((n+3+1)+3+1:n_non0_partial+1)= D((n+3+1)+3:n_non0_partial);        % modified by czc
        DE((n+3+1)+3+1:n_non0_partial+1)=DE((n+3+1)+3:n_non0_partial);
        LE((n+3+1)+3+1:n_non0_partial+1)=LE((n+3+1)+3:n_non0_partial);
        n_non0_partial=n_non0_partial+1;
    end
     D((n+3+1)+3)= D(2);        % 由于flag_vtp=0，所以partial(Vdd)·partial(Vtp)=partial(Vdd)
    DE((n+3+1)+3)=DE(2);
    LE((n+3+1)+3)=LE(2);
    % fprintf('D mid 1P-VDD\n');
    % D
    if(flag_vtn==1)
        % this belongs to flag_vtn = 1
         D((n+3+1)+(n+3)+2+1:n_non0_partial+1)= D(n+3+1+(n+3)+2:n_non0_partial);        % modified by czc
        DE((n+3+1)+(n+3)+2+1:n_non0_partial+1)=DE(n+3+1+(n+3)+2:n_non0_partial);
        LE((n+3+1)+(n+3)+2+1:n_non0_partial+1)=LE(n+3+1+(n+3)+2:n_non0_partial);
        n_non0_partial=n_non0_partial+1;
    end
    % fprintf('D mid 1P-VTN\n');
    % D
     D((n+3+1)+(n+3)+2)= D(3);        % modified by czc
    DE((n+3+1)+(n+3)+2)=DE(3);
    LE((n+3+1)+(n+3)+2)=LE(3);
    % fprintf('D mid 2P\n');
    % D

     D((n+3+1)+(n+3)+(n+2)+(n+1)+1:n_non0_partial+(n+1))= D((n+3+1)+(n+3)+(n+2)+1:n_non0_partial);        % modified by czc
    DE((n+3+1)+(n+3)+(n+2)+(n+1)+1:n_non0_partial+(n+1))=DE((n+3+1)+(n+3)+(n+2)+1:n_non0_partial);
    LE((n+3+1)+(n+3)+(n+2)+(n+1)+1:n_non0_partial+(n+1))=LE((n+3+1)+(n+3)+(n+2)+1:n_non0_partial);
    % fprintf('D mid 3P\n');
    % D
     D((n+3+1)+(n+3)+(n+2)+(n+1)+1:(n+3+1)+(n+3)+(n+2)+1)= D(1);        % modified by czc
    DE((n+3+1)+(n+3)+(n+2)+(n+1)+1:(n+3+1)+(n+3)+(n+2)+1)=DE(1);
    LE((n+3+1)+(n+3)+(n+2)+(n+1)+1:(n+3+1)+(n+3)+(n+2)+1)=LE(1);
    n_non0_partial=n_non0_partial+(n+1);
    % fprintf('D lastP\n');
    % D 
end
% AFTER what above does, D, DE, LE's length become n_partial
% fprintf('chk3 n_non0_partial: %i - n_partial:%i\n',n_non0_partial,n_partial);


% pause(10000);

% Dynamic Energy
DE(1:n_partial)=((DE(1:n_partial)-LE(1:n_partial))*10*Delay_i(i_stage)/2)/calib_dyn;
% Leakage Energy
LE(1:n_partial)=(LE(1:n_partial)*max(Delay_i)*delay_weight)/calib_leak;

%if(min(DE(1:n_partial))<0)
 %   disp('Dynamic power<0');
  %  fde0=1;
    %break;
%     continue;
%end
% D(1:n_partial)=log10(D(1:n_partial));
D(1:n_partial)=log10(D(1:n_partial)/calib_delay);
E(1:n_partial)=log10(LE(1:n_partial)+activity*DE(1:n_partial));

for i=1:3
    d_D(i,1)=(D(i+1)-D(1))/dv;
    d_E(i,1)=(E(i+1)-E(1))/dv;
end

% if (enforce_vtn_vtp_shift_eq)       % modified by czc
%     d_D(3,1) = d_D(2,1);    % vtn_shift_change = vtp_shift_change
%     d_E(3,1) = d_E(2,1);    % vtn_shift_change = vtp_shift_change
% end
% fprintf('1\n');
% d_D

for i=4:n+3
    nf1 = max(2, ceil((1-dwl)*Xc(i,i_stage)/width_eff_mult));   % Here different from the file 'write_derivation'. I remove the flag_w.
    nf2 = ceil((1+dwu)*Xc(i,i_stage)/width_eff_mult);   % Or when flag_w=1, then nf1=nf2, the derivation computation will have an error of dividing by zero
    delta_wid = (nf2 - nf1) * width_eff_mult;
    d_D(i,1)=(D(i+1)-D(1)) / delta_wid * flag_w; 
    d_E(i,1)=(E(i+1)-E(1)) / delta_wid * flag_w;  
end



for j=1:n+3
    for i=j:n+3
        nf1 = max(2, ceil((1-dwl)*Xc(j,i_stage)/width_eff_mult));   % Here different from the file 'write_derivation'. I remove the flag_w.
        nf2 = ceil((1+dwu)*Xc(j,i_stage)/width_eff_mult);   % Or when flag_w=1, then nf1=nf2, the derivation computation will have an error of dividing by zero
        delta_wid_j = (nf2 - nf1) * width_eff_mult;

        nf1 = max(2, ceil((1-dwl)*Xc(i,i_stage)/width_eff_mult));   % Here different from the file 'write_derivation'. I remove the flag_w.
        nf2 = ceil((1+dwu)*Xc(i,i_stage)/width_eff_mult);   % Or when flag_w=1, then nf1=nf2, the derivation computation will have an error of dividing by zero
        delta_wid_i = (nf2 - nf1) * width_eff_mult;

        temp_access_deriv_calc = ((n+3+1)+(n+3+1+1-j))*j/2+i+1-j;     % derivation set中的index（即D,E,DE,LE中的index）
        % fprintf('j: %i - i:%i - whats here: %i\n',j,i,temp_access_deriv_calc);
        if(i<=3)
            % fprintf('case a\n');
            if(j<=3)
                d_D2(i,j)=((D(temp_access_deriv_calc)-D(j+1))/dv - d_D(i)) / dv;
                d_E2(i,j)=((E(temp_access_deriv_calc)-E(j+1))/dv - d_E(i)) / dv;
            else        % 应该不会进这个loop
                d_D2(i,j)=((D(temp_access_deriv_calc)-D(j+1))/dv - d_D(i)) / delta_wid_j * flag_w;
                d_E2(i,j)=((E(temp_access_deriv_calc)-E(j+1))/dv - d_E(i)) / delta_wid_j * flag_w;
            end
        else
            % fprintf('case b\n');
            if(j<=3)
                d_D2(i,j)=((D(temp_access_deriv_calc)-D(j+1))/delta_wid_i - d_D(i)) / dv * flag_w; 
                d_E2(i,j)=((E(temp_access_deriv_calc)-E(j+1))/delta_wid_i - d_E(i)) / dv * flag_w; 
            else
                d_D2(i,j)=((D(temp_access_deriv_calc)-D(j+1))/delta_wid_i - d_D(i)) / delta_wid_j * flag_w;
                d_E2(i,j)=((E(temp_access_deriv_calc)-E(j+1))/delta_wid_i - d_E(i)) / delta_wid_j * flag_w;
            end
        end
    end
end


if(flag_vdd==0)
    d_D2(:,1)=0;
    d_E2(:,1)=0;
end
if(flag_vtn==0)
    d_D2(:,2)=0;
    d_E2(:,2)=0;
end
if(flag_vtp==0)
    d_D2(:,3)=0;
    d_E2(:,3)=0;
end


for i=1:n+3
    for j=1:i-1
        d_D2(j,i)=d_D2(i,j);
        d_E2(j,i)=d_E2(i,j);
    end
end

