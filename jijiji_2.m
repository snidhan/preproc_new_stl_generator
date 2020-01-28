
file_out = '/home/jose/WORK/projects/high_re_wake/new_stl_generator/spd_ar6_a0_check.stl';

windows=0;
export=1;

addpath('/home/jose/WORK/projects/high_re_wake/new_stl_generator/distmesh')

%dx = 1/128;
%dx = 1/64;
%dx=0.03
dx = 1/10;
%dx=0.2;
%dx=0.1;
a = 3; b = 0.5; c = 0.5;                                                  
amax = a*1; bmax = b*1; cmax = c*1;                                 
fd = @(x) x(:,1).^2/a^2+x(:,2).^2/b^2+x(:,3).^2/c^2-1;                      

ah = 3; bh = 0.5; ch = 0.5;

%h = @(x) sqrt(1-x(:,1).^2/ah^2);                      

%h = @(x) sqrt(1-x(:,1).^2/ah^2)+0.1;                      
%h = @(x) 0.8*sqrt(1-x(:,1).^2/ah^2)+0.1;                      
h = @(x) 0.5*sqrt(1-x(:,1).^2/ah^2)+0.06;


%[p,t] = distmeshsurface(fd,@huniform,dx,[-amax,-bmax,-cmax ; amax,bmax,cmax]);
[p,t] = distmeshsurface(fd,h,dx,[-amax,-bmax,-cmax ; amax,bmax,cmax]);


idx = [3 2 1];
ps = p(:,idx);

%% Export 

if export

if windows
stlwrite_win(file_out,t,ps,'mode','ascii')
else
stlwrite(file_out,t,ps,'mode','ascii')
end

end
 
   
