clc,clear
close all




source = [-7]; % INPUT ELEMENT
slope = 0.8 % slope parameter
Ps = 4.16; %Power per source
ZspanM = 750; %Propagation distance
start_slop = 600;
steepness=2.038;

order_slope=2 % NOT USED
abrupt=1; abrupt_start = start_slop + 280; % NOT USED

incline=1;	max_incline = 1.3; % NOT USED
incline_slope= (max_incline-slope)/slope * ZspanM/start_slop; % NOT USED

random_z=0; % NOT USED
max_ran = 1;	ran_mag = 0.05; % NOT USED

loss=0; % NOT USED

main_slope=1 ;

 
 
 
 
M  = 30;  % number of supermodes
H=1*zeros(M);
H=H+diag(ones(1,M-1),1);
H=H+1*diag(ones(1,M-1),-1);

dslpoe = 0.0425*slope

for io=1:M/2
	H(io,io)= dslpoe*(io)^main_slope  / ((M/2)^(main_slope-1) *main_slope ) ;
end
for io=M/2+1:M
	H(io,io)= dslpoe*((M-io))^main_slope / ((M/2)^(main_slope-1) *main_slope );
end

if (loss)
	H=H + 0.0007i*diag(ones(size(H,1),1));
end



M=size(H,1);
[V,D] = eig(H);
DD = diag(real(D));
[DD,IX]=sort((DD));
D=diag(DD);
V=V(:,IX);
scatter(1:size(DD,1),DD);
plot(V(:,end));


E_in=zeros(size(H,1),1);
g=1*ones(size(H,1),1);

for si=1:length(source)
	E_in(M/2+source(si)) = sqrt(Ps);
end










% RK45 constans 
c21 = 1/4;
c31 = 3/32;      c32 = 9/32;
c41 = 1932/2197; c42 = -7200/2197; c43 = 7296/2197;
c51 = 439/216;   c52 = -8;         c53 = 3680/513;    c54 = -845/4104;
c61 = -8/27;     c62 = 2;          c63 = -3544/2565;  c64 = 1859/4104; c65 = -11/40;
cw1 = 25/216;    cw3 = 1408/2565;  cw4 = 2197/4104;   cw5 = -1/5;
cz1 = 16/135;    cz3 = 6656/12825; cz4 = 28561/56430; cz5 = -9/50;     cz6 = 2/55;
k1 = E_in*0; k2 = k1; k3 = k1; k4 = k1; k5 = k1; k6 = k1;
E_inC=E_in*0;HX_L=0;Ent=0;HX_NL=0;PX1=0;PX2=0;HXtot2=0;HX_L=0;HX_NL2=0;


Zspan = 0:1:ZspanM; % - saving points
dz = 0.01;
t_t=1;  
t=0;
zz = 1;
aa = 1;
step = Zspan(3)-Zspan(2);
tt = 0;
epsilon = 10^-9;
DDia=diag(diag(H));
scale=30;
h=H;	
		
while t<=Zspan(end)

	if ismember(aa,zz)==1  		
		H=h;

		if(incline) 
			H=H + incline_slope.* t / (ZspanM) .*DDia;
		end 
		
		if (t>start_slop)
			if(~abrupt) 
				H(M/2,M/2)= steepness*(t-start_slop)^order_slope / 100^order_slope + H(M/2,M/2);
			else
				H(M/2,M/2)= steepness*(min(t,abrupt_start)-start_slop)^order_slope / 100^order_slope + H(M/2,M/2);
			end
		end		
		
		if(random_z) 
			for ir=0:max_ran
				H(M/2+ir,M/2+ir) = H(M/2+ir,M/2+ir) + ran_mag*rand();
				if(ir>0) 	H(M/2-ir,M/2-ir) = H(M/2-ir,M/2-ir) + ran_mag*rand();	end
			end
		end


	
		E_inC(:,t_t)=E_in;
		CX_z(:,t_t) = gather(abs(V'*E_in).^2);
		Pot(:,t_t) = diag(H);

		HX_L(:,t_t) = (-(sum(D*CX_z(:,t_t),1)));
		HX_NL(:,t_t) = ((sum((1/2)*abs(E_in).^4,1))); 

		time(t_t) = t;
		t_t = t_t+1;
		disp([num2str(100*t/Zspan(end))  '%'])
		aa=-1;  
	end

	f = @(E_in,H,dz) dz*1i*( H*E_in + (abs(E_in).^2).*E_in );

	k1 = f(E_in ,H,dz);
	k2 = f(E_in + c21*k1,H,dz);
	k3 = f(E_in + c31*k1+c32*k2,H,dz);
	k4 = f(E_in + c41*k1+c42*k2+c43*k3,H,dz);
	k5 = f(E_in + c51*k1+c52*k2+c53*k3+c54*k4,H,dz);  
	k6 = f(E_in + c61*k1+c62*k2+c63*k3+c64*k4+c65*k5, H, dz); 

	W5 =  E_in + cw1*k1 + cw3*k3 + cw4*k4 + cw5*k5; 
	W6 =  E_in + cz1*k1 + cz3*k3 + cz4*k4 + cz5*k5 + cz6*k6;

	kappa = norm(W5-W6);            


	if kappa<=epsilon
		t = t+dz;
		E_in = W6;    
		dz=dz*1.1;   
        if tt-t<0
			aa=1;
			tt=tt+step;
		end
	else
		dz = dz/1.1;
	end
	
end







surf(abs(E_inC).^2); shading interp

xlim([0 ZspanM])
ylim([7 23])
view(113,32)
caxis([0 3])
ax = gca;
ax.FontSize = 14; 

eff = (abs(E_inC(M/2,end)).^2) / sum(abs(E_inC(:,end)).^2) * 100;
fprintf('Funneling efficiency: %.2f%%\n', eff);


figure
surf(real(Pot)); shading interp


(sum(abs(E_inC(:,1)).^2))/(sum(abs(E_inC(:,end)).^2))





