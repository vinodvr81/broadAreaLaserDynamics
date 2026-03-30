


clear
clc
clf

psatur =   10;
pb1 = 0.05;
pb2 = 0.005;

piasa =  0.5;
pdim = 100;
dim = 50000;
dt = 4e-2;
noise=0.0000;

tic

pd(1)=1e-3;
DcarG(1)=0.;
dcarA(1)=0.;    


for k= 1:pdim/2    
    piagr(k) = k*0.4e-1;
    
% pd(1)=1e-3;
% DcarG(1)=0.;
% dcarA(1)=0.;    
 

for i = 1:dim
    dpd(i)= dt*(2*(DcarG(i) + dcarA(i) - 1)*pd(i)+ noise*rand); 
    dDcarG(i) = -dt*pb1*((DcarG(i)*(1 + pd(i)))- piagr(k));   
    ddcarA(i) = -dt*pb2*((dcarA(i)*(1 + psatur*pd(i)))+ piasa);  
    
    pd(i+1)= pd(i)+dpd(i);
    DcarG(i+1)= DcarG(i)+dDcarG(i);
    dcarA(i+1)= dcarA(i)+ddcarA(i);
     
end

power(k)=sum(pd)(dim+1);
pd(1)= pd(dim+1)+0.001;
DcarG(1)= DcarG(dim+1)+0.001;
dcarA(1)= dcarA(dim+1)+0.001;
end

% 
for k = (pdim/2+1):pdim
    piagr(k) = piagr(k-1)- 0.4e-1;
% pd(1)=1e-3;
% DcarG(1)=0.;
% dcarA(1)=0.;    


    for i = 1:dim
    dpd(i)= dt*(2*(DcarG(i) + dcarA(i) - 1)*pd(i)+ noise*rand); 
    dDcarG(i) = -dt*pb1*((DcarG(i)*(1 + pd(i)))- piagr(k));   
    ddcarA(i) = -dt*pb2*((dcarA(i)*(1 + psatur*pd(i)))+ piasa);  
    
    pd(i+1)= pd(i)+dpd(i);
    DcarG(i+1)= DcarG(i)+dDcarG(i);
    dcarA(i+1)= dcarA(i)+ddcarA(i);
        
    end
power(k)=pd(dim+1);
pd(1)= pd(dim+1)+0.001;
DcarG(1)= DcarG(dim+1)+0.001;
dcarA(1)= dcarA(dim+1)+0.001;
end

 %plot(piagr,power,'-')
 s=10;
givennumber = 6000;
gamm=0.5;

for i=1:givennumber 
 mu(i)=0.; 
 D0(i)=0.;
 d0(i)=0.;
 intensity(i)=i*(0.2*1e-3);
end
 
 
for i=1:givennumber   
  mu(i) = (1 + gamm/(1+s*intensity(i)))*(1 + intensity(i));
  D0(i) = mu(i)/(1+intensity(i));
  d0(i) = -gamm/(1+s*intensity(i));
end


plot(mu,intensity,piagr,power)
 
    
