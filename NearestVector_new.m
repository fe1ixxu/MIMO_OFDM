function [w,S,m] = NearestVector_new(A,iteration,alpha,lo)

N=size(A,2);
w=1*ones(N,1);
w2=[w;1];

d=N;
c=N;
de=sqrt(c*lo/d);
ss=(A*w).*conj(A*w);
m=10*log10(max(ss)/mean(ss));
tstart=tic;
for i=1:iteration

    A1=A'*A;
    phi=angle(A*w);
    b=sqrt(alpha)*A'*exp(1j*phi);
    Ahat=[A1 b;b' alpha];
    if lo ==1
        w2=exp(j*angle(Ahat*[w;1]));
        w=w2(1:end-1);
    else

    Z=Ahat*[w;1];

% %%% nearest vector in LOW PAR
    Z=Z/norm(Z);
    
    for k=0:d-1

    [Zm, M]=sort(Z);
    Zm=Zm(1:d-k);
    Zm_abs=abs(Zm);
    M=M(1:d-k);
    if size(Zm_abs,1)~=size(unique(Zm_abs),1)
        continue;
        
    end
    gamma=sqrt((c-k*de^2)/(sum(abs(Zm).^2)));
    if(gamma*Zm<=de)
        w2=de*exp(1j*angle(Z));
        w2(M)=gamma*Zm;
        break;
    else
        continue; 
    end
    end

    w=w2(1:end-1);
    
    end
    ss=(A*w).*conj(A*w);
    m=[m,10*log10(max(ss)/mean(ss))];
    
end
tend=toc(tstart);
fprintf("Elapsed time is %f seconds with setting k as %d \n",tend,sqrt(lo))
S=A*w;
