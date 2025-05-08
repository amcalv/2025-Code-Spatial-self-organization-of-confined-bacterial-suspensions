clear all
close all 
clc

% Save data
datadir = './data/';
qarch   = 0;

% Dimensionless parameters
K         = 0.001/0.274;         % Half-velocity
Rinf      = 2.2;         % Dimensionless domain size
Db        = 219.95/2.4e3;
Kchi      = 0.0007/0.274;
c_crit    = 5e-6/0.274;
Gamma     = 0.92308;
kt        = 1;

% Numerical parameters
Nr   = 150;                  % Number of nodes in r
vr   = linspace(1e-5,Rinf,Nr);     % r \in [0,Rinf]
dr   = Rinf/(Nr-1);             % \Delta r
tend = 48;                  % Final time
ntf  = 600;                  % Number of time steps
vt   = linspace(0,tend,ntf); % Time vector
tol  = 1e-4;                 % Tolerance

% Initial bacterial and nutrient profile 
%b0  = ones(Nr,1).*heaviside(Rinit-vr)';
b0  = ones(Nr,1);
c0  = ones(Nr,1);

% Initial plots
figure(1)
semilogy(vr,c0,'r','LineWidth',2)
hold on
semilogy(vr,b0,'k','LineWidth',2)
axis([0 Rinf 1e-3 1e3])
hold off
drawnow

    
for nt = 1:ntf-1

   t0 = vt(nt);
   tf = vt(nt+1);

   [t,b,c] = lines(b0,c0,K,Db,Kchi,Gamma,c_crit,kt,Nr,vr,t0,tf,tol,dr);
   
   b0  = b;
   c0  = c;
    
   
   %%%%% Time/dependent plots %%%%%
   figure(1)
   semilogy(vr,b,'k','LineWidth',2)
   hold on
   semilogy(vr,c,'r','LineWidth',2)
   axis([0 Rinf 1e-6 1e2])
   xlabel('Coordinate r')
   ylabel('Bacterial and nutrient concentration b(r,t), c(r,t)')
   hold off
   drawnow
   warning off

   disp([' t = ' num2str(t,'%.9e')])

end



% Save data

if qarch    
    if ~isdir(datadir)
        mkdir(datadir)           
    end

    namearchtran = [datadir 'Example_' num2str(Xi) '.dat'];
    fid = fopen(namearchtran,'wt');
    fprintf(fid,'%.12e %.12e\n',[vt;b0]);
    fclose(fid);
end