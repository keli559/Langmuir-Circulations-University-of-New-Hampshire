function vprime=ffdfft(v)
% v has N dots
%assume v is in row.

N=length(v);
vhat=fft(v);
k=[0:N/2-1 0 -N/2+1:-1];
vprimehat=i*k.*vhat;
vprime=ifft(vprimehat);