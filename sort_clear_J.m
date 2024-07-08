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
