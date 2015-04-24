%chebshev inverse transform by using ifft function

function w=chifft(f_hat)
N=length(f_hat)-1;
f_hat(1)=2*f_hat(1);
f_hat(N+1)=2*f_hat(N+1);
f_hat_fft=[f_hat;flipud(f_hat(2:N))];
w_fft=real(ifft(f_hat_fft))*N;
w(1:N+1)=w_fft(1:N+1);
w=w';

