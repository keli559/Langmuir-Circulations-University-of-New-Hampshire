%CHFT Chebyshev fourier transform

function U=chft(v)
N= length(v)-1;
if N==0, w=0;
    return, 
end


ii=0:N-1;
v=v(:);
V=[v;flipud(v(2:N))]; % transform x -> theta
ck=ones(N+1,1);
ck(1)=ck(1)/2;
ck(N+1)=ck(N+1)/2;
U_total=real(fft(V))/N;
U=U_total(1:N+1).*ck;
