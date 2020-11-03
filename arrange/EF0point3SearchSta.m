%step 1: clears memory and graphics window 
clear
clc
clf

%step 2:initialization of parameters
psatur =   10;
pb1 = 0.05;
pb2 = 0.005;

diffusioncoeff= 1.5e-4;
dt = 1e-3;
palpha = 2;
pbeta = 0;

trans = 0.1;
propagationlength = 1;

numberofturns = 1; 
dimx = 16; 
dimz = 16;

chipsize = 20;
dx = chipsize/(dimx-1) ;
dz = propagationlength/(dimz-1);
piagr(1:dimx,1:dimz) = 1.44;
piasa(1:dimx,1:dimz) = 0.5;
% piagr(1:dimx,1:dimz) = 0.;
% piasa(1:dimx,1:dimz) = 0.;
% 
% 
% for j = 1:dimx
% for i = 1:3*dimz/4    
% piagr(j,i) = 1.47;
% end
% for i = 3*dimz/4:dimz
%  piagr(j,i) = 0.;
% end
% end
% 
% 
% for j = 1:dimx
% for i = 1:3*dimz/4    
%    piasa(j,i) = 0.;
% end
% for i = 3*dimz/4:dimz
%     piasa(j,i) = 0.5;
% end
% end

        


pDcarG(1:dimx,1:dimz) =0.;
pdcarA(1:dimx,1:dimz) = 0.;
pef(1:dimx,1:dimz) =(0+0*1i);
pnoise(1:dimx,1:dimz) = 0.;

pefieldFFT(1:dimx) = (0+0*1i);
pefieldFFTx(1:dimx) = (0+0*1i);
pefieldFFTz(1:dimz) = (0+0*1i);
pzefieldFFT(1:dimz) = (0+0*1i);

xefieldFFT  =  pefieldFFT';
efieldFFTx = pefieldFFTx';
efieldFFTz  = pefieldFFTz'; 

zefieldFFT = pzefieldFFT';  

for j = 1:dimx
for i = 1:(dimz/2-3)    
pef(j,i) = (0.05+0.05*1i);
pDcarG(j,i) =  1.15; 
pdcarA(j,i) = -0.57;
pnoise(j,i)=0.;
end
for i = (dimz/2-2):(dimz/2+2)
pef(j,i) = (0.3+0.3*1i);
pDcarG(j,i) =  1.43; 
pdcarA(j,i) = -0.23;
pnoise(j,i)=0.0001*(1+1*1i);
end
for i = (dimz/2+3):dimz
 pef(j,i) = (0.05+0.05*1i);
pDcarG(j,i) =  1.15; 
pdcarA(j,i) = -0.57;
pnoise(j,i)=0.;
end
end

ef  = pef';
DcarG = pDcarG';
dcarA =pdcarA';
noise = pnoise';

% ef  = pef;
% DcarG = pDcarG;
% dcarA =pdcarA;


pkx(1:dimx) = 0.;
pkz(1:dimz) = 0.;


dkx = 2*pi/(dimx*dx);
dkz = 2*pi/(dimz*dz);


for i=1:dimx/2.
         pkx(i) = i * dkx; 
end
for i=dimx/2:dimx
         pkx(i) = (dimx - i) * dkx;
end

  
for i=1:dimz/2.
         pkz(i) = i * dkz; 
end
for i=dimz/2:dimz
         pkz(i) = (dimz - i) * dkz;
end

kx = pkx'; 
kz = pkz'; 

loopnum = 0;
   

 
 while numberofturns > 0
 
    %step 1: calculating the time part of electric field
    def = dt*(((1-1i*palpha).*DcarG+  (1-1i*pbeta).*dcarA-1).*ef);
     
    
    dDcarG = -dt*pb1*((DcarG.*(1 + abs(ef).^2))- piagr);   
    
    ddcarA = -dt*pb2*((dcarA.*(1 + psatur*abs(ef).^2))+ piasa); 
    
    ef = def + ef;
    DcarG = DcarG + dDcarG;
    dcarA = dcarA + ddcarA;
    
    
    %step 2: calculating the diffraction of x part of the electric field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:dimz
     xefieldFFT  =  fft(ef(1:dimx,i));
 
     
     
   efieldFFTx = xefieldFFT.*exp(-1i*kx.*kx*dt);
     

%    efieldFFTx = (real(xefieldFFT).*cos(kx.*kx*dt) - imag(xefieldFFT) .* sin(kx.*kx*dt)) + ...
%                 1i*(imag(xefieldFFT).*cos(kx.*kx*dt) + real(xefieldFFT) .* sin(kx.*kx*dt));
%  

    ef(1:dimx,i)=ifft(efieldFFTx);    
    end
    
    %step 3: calculating the diffraction of z part of  electric field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:dimx
     zzefieldFFT  =  fft(ef(i,1:dimz))';
 


   efieldFFTz = zzefieldFFT.*exp(-diffusioncoeff*kz.*kz*dt);
 

    ef(i,1:dimz)=ifft(efieldFFTz)';    
    end
    
     %step 4: calculating the diffraction of z part of  electric field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:dimx
     zefieldFFT  =  fft(ef(i,1:dimz))';
     
        efieldFFTzz  =  zefieldFFT.*exp(-1i*(kz./trans)*dt);
     
%     efieldFFTzz  =  (real(zefieldFFT).*cos((kz./trans)*dt) + imag(zefieldFFT) .* sin((kz./trans)*dt) + ...
%                 1i*(imag(zefieldFFT).*cos((kz./trans)*dt) - real(zefieldFFT) .* sin((kz./trans)*dt)));
            
    ef(i,1:dimz)=ifft(efieldFFTzz)'; 
    end
   
%step last: plotting the graphs
  
   
   loopnum = loopnum + 1;
   
   
   if mod(loopnum, 5000) == 0
        subplot(2,2,1:2)
   surf(abs(ef).^2)
   title('Electric field') 
  
   
   subplot(2,2,3)
   surf(DcarG)
   title('Carrier profile in the gain region')
   
   subplot(2,2,4)
   surf(dcarA)
   title('Carrier profile in the saturable absorber region')
   pause (0.1)
    loopnum         
   end
 end   
   
 
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   