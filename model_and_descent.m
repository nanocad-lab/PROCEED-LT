fprintf('model and descent NEW\n');
if(l~=1&&l~=nof(k-1)&&nof(k-1)>2)   % The initial point is not the most left point or not the most right point
    switch mode
        case 1      % One is the right point and another is the left point
            wp1=-J(l,2)+J(l-1,2);       % POWER DIFFERENCE with the left point
            wp2=J(l,1)-J(l-1,1);        % DELAY DIFFERENCE with the left point
            wq1=-J(l+1,2)+J(l,2);       % POWER DIFFERENCE with the right point
            wq2=J(l+1,1)-J(l,1);        % DELAY DIFFERENCE with the right point
        case 2
            wp1=-J(l,5)+J(l-1,5);
            wp2=J(l,1)-J(l-1,1);
            wq1=-J(l+1,5)+J(l,5);
            wq2=J(l+1,1)-J(l,1);
        case 3
            wp1=-J(l,2)+J(l-1,2);
            wp2=J(l,5)-J(l-1,5);
            wq1=-J(l+1,2)+J(l,2);
            wq2=J(l+1,5)-J(l,5);
    end
    cp=wp1+wp2;
    cq=wq1+wq2;
    
    
    [Xtemp,ntemp,Jtemp]=Gradient_descent(1, 0, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, ...
                                            2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp, ...
                                            vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
    Xkt=Xtemp;
    nk=ntemp;   % nk = ntemp = num_points_return = 2 * number_points = 10
    Jk=Jtemp;

    [Xtemp,ntemp,Jtemp]=Gradient_descent(0, 1, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
    Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);  % Another group of points after gradient descent
    Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
    nk=nk+ntemp;

    [Xtemp,ntemp,Jtemp]=Gradient_descent(wp1/cp, wp2/cp, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
    Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);
    Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
    nk=nk+ntemp;

    % modified by czc, wq1/cp -> wq1/cq, wq2/cp -> wq2/cq
    [Xtemp,ntemp,Jtemp]=Gradient_descent(wq1/cq, wq2/cq, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale );
    Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);
    Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
    nk=nk+ntemp;

    for i_nk=1:nk
        Xkt(i_nk).X(1,n_bins+1)=Rt;
        Xkt(i_nk).X(2,n_bins+1)=0;
        Xkt(i_nk).X(3,n_bins+1)=0;
    end

else
    if(nof(k-1)==2)
        switch mode
            case 1
                wp1=-J(2,2)+J(1,2);
                wp2=J(2,1)-J(1,1);
            case 2
                wp1=-J(2,5)+J(1,5);
                wp2=J(2,1)-J(1,1);
            case 3
                wp1=-J(2,2)+J(1,2);
                wp2=J(2,5)-J(1,5);
        end
        cp=wp1+wp2;
        
        [Xtemp,ntemp,Jtemp]=Gradient_descent(1, 0, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
        Xkt=Xtemp;
        Jk=Jtemp;
        nk=ntemp;

        [Xtemp,ntemp,Jtemp]=Gradient_descent(0, 1, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
        Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);
        Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
        nk=nk+ntemp;

        [Xtemp,ntemp,Jtemp]=Gradient_descent(wp1/cp, wp2/cp, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
        Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);
        Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
        nk=nk+ntemp;

        for i_nk=1:nk
            Xkt(i_nk).X(1,n_bins+1)=Rt;
            Xkt(i_nk).X(2,n_bins+1)=0;
            Xkt(i_nk).X(3,n_bins+1)=0;
        end
    else
        wp1=0.5;
        wp2=0.5;
        cp=wp1+wp2;
        
        [Xtemp,ntemp,Jtemp]=Gradient_descent(1, 0, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
        Xkt=Xtemp;
        Jk=Jtemp;
        nk=ntemp;

        [Xtemp,ntemp,Jtemp]=Gradient_descent(0, 1, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
        Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);
        Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
        nk=nk+ntemp;

        [Xtemp,ntemp,Jtemp]=Gradient_descent(wp1/cp, wp2/cp, Xc(:,1:n_bins), Rt/2, dvi, dwui, dwli, Derivative, 2*optimizing_stage, optimizing_stage_weight, flag_w, flag_vdd,flag_vtn,flag_vtp,vdd_allocation,vtn_allocation,vtp_allocation,ub_i,lb_i, mode, max_limit_metric,metric_scale);
        Xkt((nk+1):(nk+ntemp))=Xtemp(1:ntemp);
        Jk((nk+1):(nk+ntemp),1:5)=Jtemp(:,1:5);
        nk=nk+ntemp;

        for i_nk=1:nk
            Xkt(i_nk).X(1,n_bins+1)=Rt;
            Xkt(i_nk).X(2,n_bins+1)=0;
            Xkt(i_nk).X(3,n_bins+1)=0;
        end
    end
end
Jk(:,2)=Jk(:,2)-log10(delay_weight);        % WHY????
switch mode %%select colomn to rank data
    case 1
        rankc1 = 1; %delay
        rankc2 = 2; %power
    case 2
        rankc1 = 5;
        rankc2 = 1; %area
    case 3
        rankc1 = 5;
        rankc2 = 2; %power
end
if(Rt<Rt_th)
    for i_nk=1:nk
        Jk(i_nk,3)=i_nk;
    end
    J_temp=[];
    J_temp=sortrows(Jk,rankc2); %rank by power in ascending order
    i_temp=1;
    while(i_temp<nk)
        if(isnan(J_temp(i_temp+1,rankc2)))
            nk=i_temp;
            break;
        end
        if(J_temp(i_temp,rankc1)<J_temp(i_temp+1,rankc1))
            J_temp(i_temp+1,:)=[];
            nk=nk-1;
        else
            i_temp=i_temp+1;
        end
    end
    X_temp=[];
    % nk
    for i_nk=1:nk
        X_temp(i_nk).X=Xkt(J_temp(i_nk,3)).X;
        % X_temp(i_nk).X
    end
    Xkt=[];
    Xkt=X_temp;
    Jk=[];
    Jk=J_temp;
end