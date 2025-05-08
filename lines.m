function [tout,b_out,c_out] = lines(bac,c,K,Db,Kchi,Gamma,c_crit,kt,Nr,vr,t0,tf,tol,dr)

options = odeset('RelTol',tol,'AbsTol',tol);

[t15s,P15s] = ode15s(@f_dF,[t0 tf],[bac c],options);

tout  = t15s(end);
b_out = P15s(end,1:Nr);
c_out = P15s(end,Nr+1:2*Nr);

function dF = f_dF(t,vin)

    % A*dF = b, where dF = [du0/dr,...,duN/dr]

    
    A  = eye(2*Nr);
    b  = zeros(2*Nr,1);
    vb = vin(1:Nr);
    vc = vin(Nr+1:2*Nr);

    A(Nr+1:2*Nr,Nr+1:2*Nr) = diag(vr).*A(Nr+1:2*Nr,Nr+1:2*Nr);

    A(1,1:3)   =    [3 -4 1];
    A(Nr,Nr-2:Nr) = [1 -4 3];
    A(Nr+1,Nr+1:Nr+3) = [3 -4 1];
    A(2*Nr,2*Nr-2:2*Nr) = [1 -4 3];

    % Terry Hwa's PNAS sensing function. cm = 0 True K-S model
    chemder1 = Kchi./(Kchi+vc).^2;
    chemder2 = -2*Kchi./(Kchi+vc).^3;

    funct = 1/2*(1+tanh((vc(2:Nr-1)-c_crit)/0.8e-6));

    b(2:Nr-1)  = funct.*(Db*vr(2:Nr-1)'.*(vb(3:Nr)-2*vb(2:Nr-1)+vb(1:Nr-2))./(dr^2)+...
                Db*(vb(3:Nr)-vb(1:Nr-2))/(2*dr)-...
                vr(2:Nr-1)'.*(vb(3:Nr)-vb(1:Nr-2))/(2*dr).*chemder1(2:Nr-1).*(vc(3:Nr)-vc(1:Nr-2))/(2*dr) -...
                vr(2:Nr-1)'.*vb(2:Nr-1).*(vc(3:Nr)-2*vc(2:Nr-1)+vc(1:Nr-2))/(dr^2).*chemder1(2:Nr-1)-...
                vr(2:Nr-1)'.*vb(2:Nr-1).*((vc(3:Nr)-vc(1:Nr-2))/(2*dr)).^2.*chemder2(2:Nr-1)-...
                vb(2:Nr-1).*((vc(3:Nr)-vc(1:Nr-2))/(2*dr)).*chemder1(2:Nr-1));
    
    b(Nr+2:2*Nr-1) = Gamma^(-1)*(vr(2:Nr-1)'.*(vc(3:Nr)-2*vc(2:Nr-1)+vc(1:Nr-2))./(dr^2) +...
                    (vc(3:Nr)-vc(1:Nr-2))/(2*dr)-funct.*vr(2:Nr-1)'.*vb(2:Nr-1).*vc(2:Nr-1)./(vc(2:Nr-1)+K));

    b(2*Nr) = -kt*(vc(Nr)-1)*dr;
    
    dF   = A\b;

end

end