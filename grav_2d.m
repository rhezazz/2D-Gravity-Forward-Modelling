%2D gravity modelling
%Mohammad Rheza Zamani
%Reference : Last, B.J. and Kubik, K. (1983) Compact Gravity Inversion. Geophysics, 48, 713-721.http://dx.doi.org/10.1190/1.1441501
clear all;
clc;
G = 6.67428*10^-11; %Gravity Constant (m^3.kg^-1.s^-2)
lengthx = 1000;
lengthz = 1000;
%Block dimenssion
dx = 20;
dz = 20;
%Block middle point 
dmx = dx/2;
dmz = dz/2;
%Number of vertical and horizontal block
nx = lengthx/dx;
nz = lengthz/dz;
%Total number of block
nb = nx*nz;
%Contrast density
V = zeros(nz,nx);
%V(24:25,1:12) = 1000;
%V(26:27,12:22) = 2000;
%V(28:29,22:32) = 3000;
%V(30:31,32:42) = 4000;
%V(32:33,42:50) = 5000;
V(1:12,24:25) = 1000;
V(12:22,26:27) = 2000;
%Make block model
for i = 1 : nx
    x(i) = dx*i - dmx;
end
xx = repmat(x,nz,1);
for j = 1 : nz
    z(j) = dz*j - dmz;
end
z1 = z';
zz = repmat(z1,1,nx);

%Kernel matrix
for i=1:nx
    for j = 1:nb
        r1 = sqrt((zz(j)-dz/2).^2 + (x(i)-xx(j)+dx/2).^2);
        r2 = sqrt((zz(j)+dz/2).^2 + (x(i)-xx(j)+dx/2).^2);
        r3 = sqrt((zz(j)-dz/2).^2 + (x(i)-xx(j)-dx/2).^2);
        r4 = sqrt((zz(j)+dz/2).^2 + (x(i)-xx(j)-dx/2).^2);
        theta1 = atan((x(i)-xx(j)+dx/2)/(zz(j)-dz/2));
        theta2 = atan((x(i)-xx(j)+dx/2)/(zz(j)+dz/2));
        theta3 = atan((x(i)-xx(j)-dx/2)/(zz(j)-dz/2));
        theta4 = atan((x(i)-xx(j)-dx/2)/(zz(j)+dz/2));    
        Kernell(i,j) = (2*G*((x(i)-xx(j)+dx/2)*log((r2*r3)/(r1*r4)) + dx*log(r4/r3) - (zz(j)+dz/2)*(theta4-theta2)+ (zz(j)-dz/2)*(theta3-theta1))).*10^5;
    end 
end

%Plot kernel matrix
figure(1)
imagesc(Kernell)
set(gcf, 'Position', get(0, 'Screensize'));
ylabel('Observation Points','FontWeight','bold','FontSize',10)
xlabel('The Blocks','FontWeight','bold','FontSize',10)
cb = colorbar;
cb.Label.String = 'Value';
cb.Location = 'southoutside';
colormap(jet)
%Calculated gravity response
V_rs = reshape(V,nb,1);
G = Kernell*V_rs;

figure(2)
subplot(2,1,1)
plot(x,G,'*-r')
xlabel('Distance (m)','FontWeight','bold','FontSize',10)
ylabel('Gravity Anomaly (mGal)','FontWeight','bold','FontSize',10)
title('Gravity Response','FontWeight','bold','FontSize',10)
grid on
subplot(2,1,2)
s = pcolor(xx,zz,V);
s.FaceColor = 'interp';
xlabel('Distance (m)','FontWeight','bold','FontSize',10)
ylabel('Depth(m)','FontWeight','bold','FontSize',10)
title('Subsurface Model','FontWeight','bold','FontSize',10)
cb = colorbar;
cb.Label.String = 'Contrast density (kg/m^3)';
cb.Location = 'southoutside';
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'Ydir','reverse')
colormap(jet)

data = [x' G];
saveas(figure(1),'Kernell gravity.png')
saveas(figure(2),'Model Gravity.png')
writematrix(data,'Data forward gravity.dat','Delimiter','tab')