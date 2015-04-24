function G_Omega=GOmega(Psi, Omega, T, u, vs, epsilon, zeta,...
    xi, La_t, M, N)

for jj=1:M
    Psi_Z(:,jj)=chedfft(Psi(:,jj));
    Omega_Z(:,jj)=chedfft(Omega(:,jj));
    us_Z(:,jj)=1/xi*ones(N+1,1);%For given sample testing case us_z=1 all the time, therefore us_Z=1/xi
    u_Z(:,jj)=chedfft(u(:,jj));
end
for jj=1:M
    Psi_ZZ(:,jj)=chedfft(Psi_Z(:,jj));
end

for ii=1:N+1
    Omega_Y(ii,:)=ffdfft(Omega(ii,:));
    Psi_Y(ii,:)=ffdfft(Psi(ii,:));
    T_Y(ii,:)=ffdfft(T(ii,:));
    u_Y(ii,:)=ffdfft(u(ii,:));
    us_Y(ii,:)=zeros(1,M);% us doesn't depend on y according to given sample case
end
for ii=1:N+1
    Psi_YY(ii,:)=ffdfft(Psi_Y(ii,:));
end
for ii=1:N+1
    for jj=1:M
        vsPsi_YY(ii,jj)=vs(ii,jj)*Psi_YY(ii,jj);
        vsPsi_ZZ(ii,jj)=vs(ii,jj)*Psi_ZZ(ii,jj);
    end
end

for ii=1:N+1
    vsPsi_YYY(ii,:)=ffdfft(vsPsi_YY(ii,:));
    vsPsi_YZZ(ii,:)=ffdfft(vsPsi_ZZ(ii,:));
end

% advec represents the advection term :
%       advec=Psi_Z*T_Y-Psi_Y*T_Z
for ii=1:N+1
    for jj=1:M
        advPsiOmega(ii,jj)=Psi_Z(ii,jj)*Omega_Y(ii,jj)...
            -Psi_Y(ii,jj)*Omega_Z(ii,jj);
        advuus(ii,jj)=u_Z(ii,jj)*us_Y(ii,jj)-u_Y(ii,jj)*us_Z(ii,jj);
    end
end
for ii=1:N+1
    for jj=1:M
        G_Omega(ii,jj)=-epsilon*zeta*xi*advPsiOmega(ii,jj)+...
            (epsilon*zeta/La_t^2)*(xi*advuus(ii,jj)+zeta^2*vsPsi_YYY(ii,jj)...
            +xi^2*vsPsi_YZZ(ii,jj));
    end
end
