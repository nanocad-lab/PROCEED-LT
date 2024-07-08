function [X, f]=Logarithmic_barrier_N(Xi,f,G,H,ub,lb, G_l,H_l, max_limit_metric)
%%%%%%Newton step
term=0;
[V D]=eig(H);       % Calculate the eigenvalues and eigenvectors of H
D=diag(D);  % Extract the eigenvalues from the diagonal matrix
D(find(D<0))=0;
H_=V*diag(D)*V';%%%%%%%modify H to be positive semi-definite
alpha=0.5;
beta=0.75;
m=length(Xi);
t=m*2/0.02;
X=Xi;
L_G_maxlimit = - (G_l + H*X) / ( X'*G_l + 0.5*X'*H_l*X - max_limit_metric) ;%Max metric limit
L_G=-1./(X-ub)-1./(X-lb) + t*G  +L_G_maxlimit;
L_H_maxlimit = - H_l / ( X'*G_l + 0.5*X'*H_l*X - max_limit_metric) + (G_l + H*X) * (G_l+H*X)' / ( X'*G_l + 0.5*X'*H_l*X - max_limit_metric)^2;
L_H=diag(1./(X-lb).^2+1./(X-ub).^2)+t*H_ + L_H_maxlimit;
Xnt=-L_H\L_G;
Decrease=-L_G'*Xnt;
n_iter=0;
while(Decrease>t*0.001)
    t_=1;
    while(L_G'*(t_*Xnt)+1/2*(t_*Xnt)'*L_H*(t_*Xnt)>alpha*L_G'*t_*Xnt || min(X+t_*Xnt<ub-0.0001)==0 ||min(X+t_*Xnt>lb+0.0001)==0 ||min(X+t_*Xnt>lb+0.0001)==0 || min( (X+t_*Xnt)'*G_l + 0.5*(X + t_*Xnt)'*H*(X+ t_*Xnt)<max_limit_metric ) ==0 )
        t_=t_*beta;
        n_iter=n_iter+1;
        if(n_iter>10000)
            term=1;
            disp('error: number of iteration > 10000 in logic barrier');
            break;
        end
    end
    if(term~=1)
        X=X+t_*Xnt;
    end
    L_G=t*G+t*H_*(X-Xi)-1./(X-ub)-1./(X-lb) + L_G_maxlimit;
    L_H=t*H_+diag(1./(X-lb).^2+1./(X-ub).^2) + L_H_maxlimit;
    Xnt=-L_H\L_G;
    Decrease=-L_G'*Xnt;
    if(term)
        break;
    end
end
roundn = @(x,n) 10.^n .* round(x/10.^n);
X=roundn(X,-4);

f=f+(X-Xi)'*G+1/2*(X-Xi)'*H*(X-Xi);

