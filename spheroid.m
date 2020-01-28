function [coords,dv,Ix,Iy,Vol]=spheroid(ratio,d,dx)
%
% ratio = a/b;
% d semiminor/semimajor axis (depndending on ratio)
% dx Eulerian grid spacing
%
% Ix, Iy=Iz are the moments of inertia divided by the mass
%
addpath distmesh
%
% prolate spheroid --> a > b=c for a prolate spheroid
% oblate spheroid  --> a < b=c for an oblate spheroid
% we always set the smallest axis as reference length
% setting the value to d corresponds to the radius of a sphere of diameter 2d
if ratio >1
    b = d; 
    a = b*ratio;
else
    a = d;
    b = a/ratio;
end
c=b;
%
Ix = 2/5*b^2;
Iy = 1/5*(a^2+b^2);
%
Vol = 4/3*pi*a*b^2;
%
if ratio >1
    e = sqrt(1-b^2/a^2);
elseif ratio <1
    e = sqrt(1-a^2/b^2);
end
%
% this equation for the area is valid only for a prolate spheroid
if ratio >1
    S = 2*pi*b^2*(1+a/(c*e)*asin(e));
elseif ratio <1
    S = 2*pi*b^2*(1+(1-e^2)/e*atanh(e));
elseif ratio ==1
    S = 4*pi*b^2; % this is a sphere
end
%
% grid spacing --> this will be a external input
%dx = 1/15;
% target lagrangian points
% !!!!
% REVISON 27/06/2106 ga: factor increased till 1.1 instead of 1.05
% for convergence of the oblate spheroid
% !!!!
factor = 1.1; % we add 5% (10%) more points than required
nlagr = round(factor*S/dx^2);
%
% we allow a tolerance of 5%. This way the number
% of Lagrangian points will be between the required value
% and 10% more than that. This should be OK.
ntol = round((factor-1)*nlagr);
%
% this defines the bounding box
amax = a*1.2;
bmax = b*1.2;
cmax = c*1.2;
%
% initial value for iteration
% if dx is large, you might try to decrease the multiplying factor below
h = dx*.90;
%
%
% this is the spheroid formula
fd=@(x) x(:,1).^2/a^2+x(:,2).^2/b^2+x(:,3).^2/c^2-1;
%
iter = 0;
flag=true;
while flag
    %
    iter = iter+1;
    %
    tic
    [p,t]=distmeshsurface(fd,@huniform,h,[-amax,-bmax,-cmax; amax,bmax,cmax]);
    toc
    %
    % mgv: I have tried to decrease the tolerance and the results did not improve
    % much
    %
    ntri   = size(t,1);
    nerr = ntri-nlagr;
    display([iter,ntri,nlagr])
    fac = 1 + (nerr/nlagr)/3;
    if abs(nerr) > ntol
        h = h*fac;
    else
        flag = false;
    end
end
%
area_tri = zeros(ntri,1);
coords   = zeros(ntri,3);
%
for it = 1:ntri
    ii = t(it,:);
    %
    v1 = p(ii(2),:)-p(ii(1),:);
    v2 = p(ii(3),:)-p(ii(1),:);
    v3  = cross(v1,v2);
    area_tri(it) = norm(v3)/2;
    coords(it,:) = (p(ii(1),:)+p(ii(2),:)+p(ii(3),:))/3;
end
%hold on
%plot3(coords(:,1),coords(:,2),coords(:,3),'k.')
%hold off
dv = area_tri*dx;

