function [X, f]=Logic_barrier_N(Xi,f,G,H,ub,lb)
%%%%%%Newton step
[V D]=eig(H);
D=diag(D);
D(find(D<0))=0;
H_=V*diag(D)*V';%%%%%%%modify H to be positive semidefinite
alpha=0.5;
beta=0.75;
m=length(Xi);
t=m*2/0.02;
X=Xi;
L_G=-1./(X-ub)-1./(X-lb) + t*G;
L_H=diag(1./(X-lb).^2+1./(X-ub).^2)+t*H_;
Xnt=-L_H\L_G;
Decrease=-L_G'*Xnt;
n_iter=0;
while(Decrease>t*0.001)
    t_=1;
    while(L_G'*(t_*Xnt)+1/2*(t_*Xnt)'*L_H*(t_*Xnt)>alpha*L_G'*t_*Xnt || min(X+t_*Xnt<ub)==0 ||min(X+t_*Xnt>lb)==0)
        t_=t_*beta;
        n_iter=n_iter+1;
        if(n_iter>100000)
            disp('error: number of iteration > 100000 in logic barrier');
            break;
        end
    end
    X=X+t_*Xnt;
    L_G=t*G+t*H_*(X-Xi)-1./(X-ub)-1./(X-lb);
    L_H=t*H_+diag(1./(X-lb).^2+1./(X-ub).^2);
    Xnt=-L_H\L_G;
    Decrease=-L_G'*Xnt;
end
f=f+(X-Xi)'*G+1/2*(X-Xi)'*H*(X-Xi);

