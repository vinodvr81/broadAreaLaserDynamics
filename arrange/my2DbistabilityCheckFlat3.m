

%step 1: clears memory and graphics window 
clear
clc
clf



%step 2:initialization of parameters
psatur =   10;
pb1 = 0.05;
pb2 = 0.005;

diffusioncoeff= 1.5e-4;
dt = 1e-2;
palpha = 2;
pbeta = 0;

trans = 0.1;
propagationlength = 1;


dimx = 16; 
dimz = 16;
dimp = 100;
chipsize = 20;
dx = chipsize/(dimx-1) ;
dz = propagationlength/(dimz-1);
% 
piagr(1:dimx,1:dimz) = 0.;
 piagr1(1:dimp,1:dimx,1:dimz) = 0.;
 piasa(1:dimx,1:dimz) = 0.;
%   piagr(1:dimx,1:dimz) = 1.44;
%  piasa(1:dimx,1:dimz) = 0.5;
% % %  
% piagr1(1:dimx,1:dimz) =1.44;



 piasa(1:dimx,1:dimz) = 0.5;
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
% 
for j = 1:dimx
for i = 1:(dimz/2-3)    
pef(j,i) = (0.05+0.05*1i);
pDcarG(j,i) =  1.15; 
pdcarA(j,i) = -0.57;
end
for i = (dimz/2-2):(dimz/2+2)
pef(j,i) = (0.3+0.3*1i);
pDcarG(j,i) =  1.43; 
pdcarA(j,i) = -0.23;
end
for i = (dimz/2+3):dimz
 pef(j,i) = (0.05+0.05*1i);
pDcarG(j,i) =  1.15; 
pdcarA(j,i) = -0.57;
end
end


% for j = 1:dimx
% for i = 1:(dimz/2-3)    
% pef(j,i) = (0.013+0.013*1i);
% pDcarG(j,i) =  0.; 
% pdcarA(j,i) = 0.;
% end
% for i = (dimz/2-2):(dimz/2+2)
% pef(j,i) = (0.013+0.013*1i);
% pDcarG(j,i) =  0.; 
% pdcarA(j,i) = 0.;
% end
% for i = (dimz/2+3):dimz
% pef(j,i) = (0.013+0.013*1i);
% pDcarG(j,i) =  0.; 
% pdcarA(j,i) = 0.;
% end
% end



% ef  = pef';
% DcarG = pDcarG';
% dcarA =pdcarA';
% noise = pnoise';

ef  = pef;
DcarG = pDcarG;
dcarA =pdcarA;

pkx(1:dimx) = 0.;
pkz(1:dimz) = 0.;


dkx = 2*pi/(dimx*dx);
dkz = 2*pi/(dimz*dz);


for i=1:(dimx/2)
         pkx(i) = i * dkx; 
end
for i=(dimx/2)+1:dimx
         pkx(i) = (dimx - i) * dkx;
end

  
for i=1:dimz/2
         pkz(i) = i * dkz; 
end
for i=dimz/2:dimz
         pkz(i) = (dimz - i) * dkz;
end

kx = pkx'; 
kz = pkz'; 



% figure1 = figure('Colormap',...
%     [0 0 0;0 0 0.0500000007450581;0 0 0.100000001490116;0 0 0.150000005960464;0 0 0.200000002980232;0 0 0.25;0 0 0.300000011920929;0 0 0.349999994039536;0 0 0.400000005960464;0 0 0.449999988079071;0 0 0.5;0 0 0.550000011920929;0 0 0.600000023841858;0 0 0.649999976158142;0 0 0.699999988079071;0 0 0.75;0 0 0.800000011920929;0 0 0.850000023841858;0 0 0.899999976158142;0 0 0.949999988079071;0 0 1;0 0.0588235296308994 1;0 0.117647059261799 1;0 0.176470592617989 1;0 0.235294118523598 1;0 0.294117659330368 1;0 0.352941185235977 1;0 0.411764711141586 1;0 0.470588237047195 1;0 0.529411792755127 1;0 0.588235318660736 1;0 0.647058844566345 1;0 0.705882370471954 1;0 0.764705896377563 1;0 0.823529422283173 1;0 0.882352948188782 1;0 0.941176474094391 1;0 1 1;0.0666666701436043 1 0.933333337306976;0.133333340287209 1 0.866666674613953;0.200000002980232 1 0.800000011920929;0.266666680574417 1 0.733333349227905;0.333333343267441 1 0.666666686534882;0.400000005960464 1 0.600000023841858;0.466666668653488 1 0.533333361148834;0.533333361148834 1 0.466666668653488;0.600000023841858 1 0.400000005960464;0.666666686534882 1 0.333333343267441;0.733333349227905 1 0.266666680574417;0.800000011920929 1 0.200000002980232;0.866666674613953 1 0.133333340287209;0.933333337306976 1 0.0666666701436043;1 1 0;1 0.888888895511627 0;1 0.777777791023254 0;1 0.666666686534882 0;1 0.555555582046509 0;1 0.444444447755814 0;1 0.333333343267441 0;1 0.222222223877907 0;1 0.111111111938953 0;1 0 0;0.75 0 0;0.5 0 0]);
% 
% % Create axes
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.583837209302326 0.775 0.341162790697674]);
% hold(axes1,'all');
% 
% subplot1 = subplot(2,2,3,'Parent',figure1);
% hold(subplot1,'all');
% 
% subplot2 = subplot(2,2,4,'Parent',figure1,'CLim',[-0.5535 -0.1525]);
% hold(subplot2,'all');
% 
%  fig=figure;
%  aviobj = avifile('bi2DIntiCond3.avi')

 noise = 0.0001;
power(1:dimp) = 0.;
 for k = 1:dimp/2

 piagr(1:dimx,1:dimz) = k*0.04; 
 piagr1(k,1:dimx,1:dimz) = piagr(1:dimx,1:dimz);
 loopnum = 0;
     
 while loopnum < 5000*10
 
    %step 1: calculating the time part of electric field
    def = dt*(((1-1i*palpha).*DcarG+  (1-1i*pbeta).*dcarA-1).*ef +noise*rand(dimx,dimz));
     
    
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
%     
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
   %intens(loopnum) = sum(sum(abs((ef).^2)));
%    intens1(loopnum) = sum(abs(ef(1:dimx,1).^2));
%    intens2(loopnum)= sum(abs(ef(1:dimx,16).^2));
%     
   
%    if mod(loopnum, 5000) == 0
% %     plot(intens)
% %      pause(0.1)
% %     loopnum
%    
% %        
% %         subplot(2,2,1:2)
% %    surf(abs(ef).^2)  
% %    title('Electric field') 
% %       
% % 
% %    subplot(2,2,3)
% %    surf (DcarG)
% %    title('Carrier profile in the gain region')
% %    
% %    subplot(2,2,4)
% %    surf (dcarA)
% %    title('Carrier profile in the saturable absorber region')
% %    pause(0.1)
% % %    pause (0.1)
% % %    
% %    F = getframe(fig);
% %     aviobj = addframe(aviobj,F); 
%   
%    end
 

 end 
 power(k) =  sum(sum(abs((ef).^2)));
 
 ef= ef;
DcarG= DcarG;
dcarA= dcarA ;
  end
   
 
     for k = (dimp/2)+1:dimp  
        
         piagr(1:dimx,1:dimz) =  piagr(1:dimx,1:dimz)- 0.04;
          piagr1(k,1:dimx,1:dimz) = piagr(1:dimx,1:dimz);
         loopnum = 0;
 
   while loopnum < 5000*10
 
    %step 1: calculating the time part of electric field
    def = dt*(((1-1i*palpha).*DcarG+  (1-1i*pbeta).*dcarA-1).*ef +noise*rand(dimx,dimz));
     
    
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
%     
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
  %intens(loopnum) = sum(sum(abs((ef).^2)));
%    intens1(loopnum) = sum(abs(ef(1:dimx,1).^2));
%    intens2(loopnum)= sum(abs(ef(1:dimx,16).^2));
%     
   
%    if mod(loopnum, 5000) == 0
% %     plot(intens)
% %      pause(0.1)
% %     loopnum
%    
% %        
% %         subplot(2,2,1:2)
% %    surf(abs(ef).^2)  
% %    title('Electric field') 
% %       
% % 
% %    subplot(2,2,3)
% %    surf (DcarG)
% %    title('Carrier profile in the gain region')
% %    
% %    subplot(2,2,4)
% %    surf (dcarA)
% %    title('Carrier profile in the saturable absorber region')
% %    pause(0.1)
% % %    pause (0.1)
% % %    
% %    F = getframe(fig);
% %     aviobj = addframe(aviobj,F); 
%   
    end
 
 power(k) = sum(sum(abs((ef).^2)));
 
 
 ef= ef;
DcarG= DcarG;
dcarA= dcarA ;
 end 
   
   
 
plot(piagr1(1:100,1,1),power,'-*')
title('power vs injection')

   
   
   
   
   
   
   
   
   
   
   
   