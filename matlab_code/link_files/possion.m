function phi=possion(M,N,g,xi,zeta)
z=cos((0:N)'.*pi/N);
y=(0:M-1)'.*(2.0*pi)/M;
clear ii jj
%set up g(z,y) and g_yhat(z,k)
%for ii=1:N+1
%    for jj=1:M
%        g(ii,jj)=-((pi^2.0)/4.0+1)*cos((pi/2.0)*z(ii))*cos(y(jj));
%    end
%end
for ii=1:N+1
    for jj=1:M
        g(ii,jj)=g(ii,jj)*(-1/xi^2);
    end
end
for ii=1:N+1
    g_yhat(ii,:)=fft(g(ii,:));%g_yhat(z,k)
end
%clear ii jj
lambda=[0:M/2-1 0 -M/2+1:-1]'.*[0:M/2-1 0 -M/2+1:-1]'*(zeta/xi)^2;
for jj=1:M
%for every k, there is a bunch of Helmholtz
%     equation both in real and imaginary parts
    phi_yhat_real(:,jj)=helmholtz(N, ...
        real(g_yhat(:,jj)), lambda(jj));
    phi_yhat_imag(:,jj)=helmholtz(N, ...
        imag(g_yhat(:,jj)), lambda(jj));
    phi_yhat(:,jj)=phi_yhat_real(:,jj)+i*phi_yhat_imag(:,jj);
end
for ii=1:N+1
    phi(ii,1:M)=real(ifft(phi_yhat(ii,1:M)));
end