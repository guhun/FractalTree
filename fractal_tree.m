function fractal_tree;
global cnt_stop lscale nran
close all

% cnt_stop = 6;
% lscale = .9;
% nran = 0;
% 
% x0 = [0;0];
% y0 = [0;1.5];
% l  = y0(2)*lscale;
% n = 3;
% theta = 80;


% Good settings:
% cnt_stop = 9;
% lscale = .9;
% nran = 8;
% stars=0;
% 
% x0 = [0;0];
% y0 = [0;1.5];
% l  = y0(2)*lscale;
% n = 10;
% theta = 80;

ABC=['abcdefghijklmnopqrstuvwxyz'];
name='tyrone';
for i=1:length(name); I(i)=find(name(i)==ABC); end
sprev = rng(sum(I),'multFibonacci');

% Good settings:
cnt_stop = 9;
lscale = .85;
nran = 13;
stars=0;

x0 = [0;0];
y0 = [0;1.5];
z0 = [0;0];
l  = y0(2)*lscale;
n = 15;
theta = 90;

cnt = 0;



[x1, y1, x2, y2, level] = fractal_tree_base(x0(2),y0(2),l,n,theta,90,cnt);

% [x1, y1, x2, y2, z1, z2, level] = fractal_tree_base3d(x0(2),y0(2),z0(2),l,n,theta,90,cnt);

maxthickness = max(level)+1;
hold off

% figure1 = figure('Color', [1 1 1]);
axes1 = axes('Visible','off','Color',[1 1 1]);
box(axes1,'on');
hold(axes1,'all');

l=6;
generate_tree_base(x0(1), y0(1), l/8, l*.05, l, 8)
if stars
    generate_random_stars( 500, [-1*l,1*l], [y0(1), 2*l], 20 )
    generate_milky_way( 20000, [-1*l,1*l], [y0(1), 2*l], 8, l/6, 5*l, [-l*3/4,-l/2],70  )
end
plot( x0, y0, 'k-','LineWidth',maxthickness+1);
hold on
% plot( x0, -y0, 'k-','LineWidth',maxthickness+1);
for i = 1:length(x1)
    plot([x1(i),x2(i)], [y1(i),y2(i)], 'k-','LineWidth',maxthickness-level(i)+1)
%     plot3([x1(i),x2(i)], [y1(i),y2(i)], [z1(i),z2(i)], 'k-','LineWidth',maxthickness-level(i)+1)
%     plot([x1(i),x2(i)], [y1(i),y2(i)], 'k-','LineWidth',maxthickness-level(i))
end
axis equal

return

% Tree 2
cnt_stop = 6;
lscale = 1.05;
nran = 5;

x0 = [4;4];
y0 = [0;1.2];
l  = y0(2)*lscale;
n = 7;
theta = 70;
cnt = 0;


[x1, y1, x2, y2, level] = fractal_tree_base(x0(2),y0(2),l,n,theta,90,cnt);

maxthickness = max(level)+1;
plot( x0, y0, 'k-','LineWidth',maxthickness+1);
hold on
for i = 1:length(x1)
    plot([x1(i),x2(i)], [y1(i),y2(i)], 'k-','LineWidth',maxthickness-level(i))
%     plot(-[x1(i),x2(i)], [y1(i),y2(i)], 'k-','LineWidth',maxthickness-level(i))
end
% axis equal

% Tree 3
cnt_stop = 6;
lscale = 0.9;
nran = 6;

x0 = [-4;-4];
y0 = [0;1.1];
l  = y0(2)*lscale;
n = 8;
theta = 120;
cnt = 0;


[x1, y1, x2, y2, level] = fractal_tree_base(x0(2),y0(2),l,n,theta,90,cnt);

maxthickness = max(level)+1;
plot( x0, y0, 'k-','LineWidth',maxthickness+1);
hold on
for i = 1:length(x1)
    plot([x1(i),x2(i)], [y1(i),y2(i)], 'k-','LineWidth',maxthickness-level(i))
%     plot(-[x1(i),x2(i)], [y1(i),y2(i)], 'k-','LineWidth',maxthickness-level(i))
end
% axis equal
axis([-10,10,0,10])



function [x1, y1, x2, y2, level] = fractal_tree_base(x0,y0,l,n,theta,theta0,cnt)
global cnt_stop lscale nran
% Generate a fractal tree
% x0 is a vector of size 2
% y0 is a vector of size 2
% l is the length of the pattern
% n is the number of splits
% theta is the angle of the fractal split


if (cnt>=cnt_stop)
    x1=[];y1=[];x2=[];y2=[];level=[];
    return
end



dtheta = theta/(n-1);
th_new = [-theta/2:dtheta:theta/2]';

% one-sided
% dtheta = theta/(n-1);
% th_new = [-theta:dtheta:0]';

% randomly select next branch
for j = 1:nran
irem = randi([1,length(th_new)],1,1);
th_new = th_new([1:irem-1,irem+1:length(th_new)]);
end

theta_new = th_new + theta0;
xnew = l*cosd( theta_new ) + x0;
ynew = l*sind( theta_new ) + y0;

x1 = repmat(x0,[length(xnew),1]);
x2 = xnew;
y1 = repmat(y0,[length(xnew),1]);
y2 = ynew;
level = repmat(cnt, [length(xnew),1]);
% x1 = x0;
% x2 = x0(2);
% y1 = y0;
% y2 = y0(2);
cnt = cnt + 1;
% lsc_ran = rand*.1
for i = 1:length(xnew)
    [x1n, y1n, x2n, y2n, L0] = fractal_tree_base(xnew(i), ynew(i), l*lscale, n, theta, theta_new(i), cnt);
    x1 = [x1; x1n];
    y1 = [y1; y1n];
    x2 = [x2; x2n];
    y2 = [y2; y2n];
    level = [level; L0];
end


function [x1, y1, z1, x2, y2, z2, level] = fractal_tree_base3d(x0,y0,z0,l,n,theta,theta0,cnt)
global cnt_stop lscale nran
% Generate a fractal tree
% x0 is a vector of size 2
% y0 is a vector of size 2
% l is the length of the pattern
% n is the number of splits
% theta is the angle of the fractal split


if (cnt>=cnt_stop)
    x1=[];y1=[];x2=[];y2=[];z1=[];z2=[];level=[];
    return
end



dtheta = theta/(n-1);
th_new = [-theta/2:dtheta:theta/2]';
phi_new = linspace(0,180,n+1);
phi_new = phi_new(1:end-1);
[th_new, phi_new] = meshgrid(th_new, phi_new);
% one-sided
% dtheta = theta/(n-1);
% th_new = [-theta:dtheta:0]';

% randomly select next branch
for j = 1:nran
irem = randi([1,length(th_new)],1,1);
th_new = th_new([1:irem-1,irem+1:length(th_new)]);
phi_new = phi_new([1:irem-1,irem+1:length(phi_new)]);
end

theta_new = th_new + theta0;

xnew = l*cosd( theta_new ).*cosd(phi_new) + x0;
ynew = l*sind( theta_new ).*cosd(phi_new) + y0;
znew = l*sind( phi_new ) + z0;



x1 = repmat(x0,[length(xnew),1]);
x2 = xnew;
y1 = repmat(y0,[length(xnew),1]);
y2 = ynew;
z1 = repmat(z0,[length(xnew),1]);
z2 = znew;

level = repmat(cnt, [length(xnew),1]);
% x1 = x0;
% x2 = x0(2);
% y1 = y0;
% y2 = y0(2);
cnt = cnt + 1;
% lsc_ran = rand*.1
for i = 1:length(xnew)
    [x1n, y1n, z1n, x2n, y2n, z2n, L0] = fractal_tree_base3d(xnew(i), ynew(i), znew(i), l*lscale, n, theta, theta_new(i), cnt);
    x1 = [x1; x1n];
    y1 = [y1; y1n];
    x2 = [x2; x2n];
    y2 = [y2; y2n];
    z1 = [z1; z1n];
    z2 = [z2; z2n];
    level = [level; L0];
end

function generate_tree_base(x0, y0, dx, dy, l, ps)

L = l;

x = x0-L/2;
y = y0;

even=1;
while L > 0

    plot([x,x+L],[y,y],'k-','LineWidth', ps);
    ps=max(ps-1,1);
        
    y=y-dy;
    L=L-dx;
    x=x0-L/2;
    dy=dy*(ps/(ps+1));
    
    
end

    
function generate_random_stars( nstars, xr, yr, sr )

for i = 1:nstars
    s = randi([1,sr],1,1);
    x = rand(1)*(xr(2)-xr(1))+xr(1);
    y = rand(1)*(yr(2)-yr(1))+yr(1);
    plot(x,y,'k.','MarkerSize',s);
end
    
function generate_milky_way( nstars, xr, yr, sr, dx, L, x0, th0  )


    for i = 1:nstars
      X = rand(1)*L;
      Y = dx*randn(1,1);
      r = sqrt(X*X+Y*Y);
      th = atand(Y/X)+th0;
      
      x = r*cosd(th)+x0(1);
      y = r*sind(th)+x0(2);
      
      s = randi([1,sr],1,1);
      
      if (x>=xr(1) && x<=xr(2) && y>=yr(1) && y<=yr(2))
        plot(x,y,'k.','MarkerSize',s);
      end

    end
    


 