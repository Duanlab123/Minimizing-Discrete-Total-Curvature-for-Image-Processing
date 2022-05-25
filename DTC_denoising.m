
close all
clear all

addpath('img'); 
addpath('util');


%%%%%%%%%%%%%%%%%%%%%  cameraman image  %%%%%%%%%%%%%%
               

lambda = 0.6; %%           %%%%%  
p1= 30;                                          
a = 1;  %                                                        
b = 5;  %                                      
                                        


iter = 300;
threshold_res = 5e-5;

uclean = double(imread('Cameraman.png'));                      
        

sigma = 20;
u0 = uclean+sigma*randn(size(uclean));


[l1, l2] = size(u0);
Area = l1*l2;

lambda11 = u0*0; lambda12 = u0*0; 
x1 = u0*0; x2 = u0*0; 
u = u0;

residual = zeros(1,iter);
energy = zeros(1,iter);
error = zeros(1,iter);
%PSNR = zeros(1,iter); SSIM = zeros(1,iter);
relative_error = zeros(1,iter);
DtD = abs(psf2otf([1,-1],[l1, l2])).^2 + abs(psf2otf([1;-1],[l1, l2])).^2;
A = (p1)*DtD + lambda;

real_iter = iter;
record_flag = 0;
t1=clock;
for i=1:iter
    
    u_old = u;
    
    g = lambda*u0 - dxb(p1*x1 + lambda11) - dyb(p1*x2 + lambda12);
    g = fftn(g);
    u = real(ifftn(g./A));
    
    %%%%%%%%%%%%%%%%%%%%%  Update afa
   
    k = TotalCuNew(u);
    afa = a + b*(abs(k));
   % afa = a + b*(abs(k).^2);
   % afa = sqrt(a + b*(abs(k).^2));
  
    %%%%%%%%%%%%%%%%  For x
    
    xx = dxf(u) - lambda11/p1;
    xy = dyf(u) - lambda12/p1;
    x = sqrt(xx.^2 + xy.^2);
    x(x==0) = 1;
    x = max(x - afa/p1,0)./x;
    x1 = xx.*x;
    x2 = xy.*x;
    
    lambda11_old = lambda11;
    lambda12_old = lambda12;
    
    lambda11 = lambda11 + p1*(x1 - dxf(u));
    lambda12 = lambda12 + p1*(x2 - dyf(u));
    
    R11 = x1 - dxf(u);
    R12 = x2 - dyf(u);
    
    energy(i) = TCEuler_energy(u,u0,x1,x2,afa,lambda);
    
    residual(i) = sum( abs(R11(:))+abs(R12(:)) )/Area;
    error(i) = sum(sum(abs(lambda11-lambda11_old)+abs(lambda12-lambda12_old)))/Area;
    relative_error(i) =  sum(sum( abs(u-u_old) ))/sum(sum(u_old)); 
 
      if( relative_error(i) < threshold_res )
        if( record_flag==0 )
            real_iter = i;
            record_flag = 1;
        end
      end   
end
t2=clock;

t=etime(t2,t1);

fprintf(' The iteration number is: %10d\n', real_iter);
fprintf(' The iteration time is: %4.2fs', t);
fprintf(' The relative error is: %10.8f\n', relative_error(real_iter));

psnr_u = psnr(uint8(u),uint8(uclean))

iternum = 1:i;
figure;
plot(iternum,log(relative_error),'r','LineWidth',2); 
xlabel('Iteration')
ylabel('Relative error in u^k')
set(gca,'FontWeight','bold')

figure;
plot(iternum,log(residual),'r','LineWidth',2);
xlabel('Iteration')
ylabel('Relative residual')
set(gca,'FontWeight','bold')

figure;
plot(iternum,log(error),'r','LineWidth',2);  
xlabel('Iteration')
ylabel('Relative error in multiplier')
set(gca,'FontWeight','bold')

figure;
plot(iternum,log(energy),'b','LineWidth',2);   
xlabel('Iteration')
ylabel('Energy')
set(gca,'FontWeight','bold')

figure;
imshow(uint8(u0));

figure;
imshow(uint8(u));
