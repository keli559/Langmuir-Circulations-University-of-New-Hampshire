%helmholtz equation solving as a function in Matlab
function u=helmholtz(N, f, lambda)

x=cos((0:N)'.*pi/N);
n=(3:N+1)';
%f_hat has N+3 terms from 1, 2,... N+3, where N, N+1, N+2, N+3 are 0
f_hat=chft(f); 
f_hat(N:N+1)=0;
f_hat=[f_hat',0,0]';

an=1./(4.*(n-1).*(n-2));%an has 2 less terms than N+1*N+1 matrix
bn=1./(2*((n-1).^2-1));%bn has N-1 terms, same as an
dn=1./(4*(n-1).*n);%N-1 terms

beta=[1 1 ones(1,N-3) zeros(1,4)]';
cn=[2 ones(1,N-2) zeros(1,4)]';

known_vec=zeros(N+1,1);
for ii=(3:N+1)
    known_vec(ii)=cn(ii-2)*an(ii-2)*f_hat(ii-2)-bn(ii-2)...
        *beta(ii)*f_hat(ii)+dn(ii-2)*beta(ii+2)*f_hat(ii+2);
end

A=zeros(N+1);

A(1,:)=ones;
A(2,:)=cos((0:N)*pi);

clear ii

for ii=(3:N+1)
    A(ii,ii-2)=-lambda*cn(ii-2)*an(ii-2);
    A(ii,ii)=1+lambda*bn(ii-2)*beta(ii);
end

clear ii

for ii=(3:N-1)
    A(ii,ii+2)=-lambda*dn(ii-2)*beta(ii+2);
end

A_spar=sparse(A);

%u_hat=dlusolve(N, A, known_vec);
u_hat=(A_spar\known_vec);


u=chifft(u_hat);
u(1)=0;
u(N+1)=0;
