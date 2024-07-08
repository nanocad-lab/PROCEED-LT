function [X, f]=Logarithmic_barrier_G(Xi,f,G,H,ub,lb, G_l,H_l, max_limit_metric)
%%%%%%%%%%%%%%%Gradient descent
f_s=0;
term=0;
alpha=0.5;
beta=0.5;
m=length(Xi);
t=m*2/0.02;
X=Xi;
L_G_maxlimit = - (G_l + H*X) / ( X'*G_l + 0.5*X'*H_l*X - max_limit_metric) ;%Max metric limit
L_G = -1./(X-ub)-1./(X-lb) + t*G + L_G_maxlimit; %%upper and lower bound 
L_H_maxlimit = - H_l / ( X'*G_l + 0.5*X'*H_l*X - max_limit_metric) + (G_l + H*X) * (G_l+H*X)' / ( X'*G_l + 0.5*X'*H_l*X - max_limit_metric)^2;
L_H=diag(1./(X-lb).^2+1./(X-ub).^2)+t*H + L_H_maxlimit;
Xnt=-L_G;
Decrease=norm(Xnt,2);
n_iter=0;
while(Decrease>t*0.01 && term==0)
    t_=1;
    while(L_G'*(t_*Xnt)+1/2*(t_*Xnt)'*L_H*(t_*Xnt)>alpha*L_G'*t_*Xnt || min(X+t_*Xnt<ub-0.0001)==0 ||min(X+t_*Xnt>lb+0.0001)==0 || min( (X+t_*Xnt)'*G_l + 0.5*(X + t_*Xnt)'*H*(X+ t_*Xnt)<max_limit_metric ) ==0 )
        t_=t_*beta;
        n_iter=n_iter+1;
        if(n_iter>10000)
            if(f_s==0)            
                f_s=1;
                n_iter=0;
                t_=0;
                break;
            else
                disp('error: number of iteration > 10000 in logic barrier Gradient');
                term=1;
                break;
            end
        end
    end
    if(term~=1)
        X=X+t_*Xnt;
    end
    L_G=t*G+t*H*(X-Xi)-1./(X-ub)-1./(X-lb) + L_G_maxlimit;
    L_H=t*H+diag(1./(X-lb).^2+1./(X-ub).^2) + L_H_maxlimit;
    Xnt=-L_G;
    Decrease=norm(Xnt,2);
    if(f_s)
        argmax=find(abs(L_G)==max(abs(L_G)));
        Xnt=zeros(size(L_G));
        Xnt(argmax)=-L_G(argmax);
        Decrease=norm(Xnt,inf);
    end
    if(term)
        break;
    end
end
roundn = @(x,n) 10.^n .* round(x/10.^n);
X=roundn(X,-4);

f=f+(X-Xi)'*G+1/2*(X-Xi)'*H*(X-Xi);
