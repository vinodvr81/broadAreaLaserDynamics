
%step 1: clears memory and graphics window 
clear
clc
clf


load n1laset2dIma.mat
 fig=figure1;
 aviobj = avifile('exampleTrialColor3dIma.avi')


while loopnum < 5000*5000
 
    %step 1: calculating the time part of electric field
    def = dt*(((1-1i*palpha).*DcarG+  (1-1i*pbeta).*dcarA-1).*ef + noise*rand(dimx,dimz));
     
    
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
   imagesc (abs(ef).^2)  
   title('Electric field') 
      

   subplot(2,2,3)
   imagesc (DcarG)
   title('Carrier profile in the gain region')
   
   subplot(2,2,4)
   imagesc (dcarA)
   title('Carrier profile in the saturable absorber region')
   pause (0.1)
   
   F = getframe(fig);
    aviobj = addframe(aviobj,F); 
  
   end
 
 end      
   
 
     close(fig)
    aviobj = close(aviobj);