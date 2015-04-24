%Gu calculation, 
%On Z, it is chebyshev, so points on Z cover both ends of the domain. N+1
%On Y, it is FFT, so last point in Y domain doesn't count. M
%Where M and N are the number of the division of Y and Z domain
%respectively.

%This function is following the given sample parameter regime without 
%rotation factor 

function G=Gu(Psi, u, vs, epsilon, zeta, xi, La_t, M, N)
% In this coding ii (1:N+1) counts on Z, jj (1:M) counts on Y. Ex, Psi(ii, jj)
% derivative of the non linear terms Psi_Z, u_Y Psi_Y, u_Z
% vs adnnd us issue:
% treat them as N+1*M matrix. This is more general in the form.
for jj=1:M
    Psi_Z(:,jj)=chedfft(Psi(:,jj));
    u_Z(:,jj)=chedfft(u(:,jj));
end
for ii=1:N+1
    u_Y(ii,:)=ffdfft(u(ii,:));
    Psi_Y(ii,:)=ffdfft(Psi(ii,:));
end
% advec represents the advection term :
%       advec=Psi_Z*u_Y-Psi_Y*u_Z
for ii=1:N+1
    for jj=1:M
        advec(ii,jj)=Psi_Z(ii,jj)*u_Y(ii,jj)-Psi_Y(ii,jj)*u_Z(ii,jj);
    end
end
for ii=1:N+1
    for jj=1:M
        G(ii,jj)=-epsilon*zeta*xi*advec(ii,jj)-(epsilon*zeta/La_t^2)...
            *vs(ii,jj)*u_Y(ii,jj)+vs(ii,jj)/La_t^2;
    end
end



