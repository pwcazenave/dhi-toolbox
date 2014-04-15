clear
close all

% constants
inc=0.5;
H=0.8;
g=9.81;

% variables
T=1:inc:5;
D=[4:inc:60];
% D=[4,6,8,16,21,29,31];

for j=1:length(T);
    % 3 equations from Xu, D. et al. (1999)
    % Xu, D. et al (1999)
    lambda(j,2)=0.86*(g/(2*pi))*(T(j)^2);
    % common inverse method
    lambda(j,1)=(g/(2*pi))*(T(j)^2);
    % Neumann-spectrum derived relationship
    lambda(j,3)=(2/3)*(g/(2*pi))*(T(j)^2);
    % Ursell number from mike documentation, I think
    for i=1:length(D);
        Us(i,j)=(H*(lambda(j,2)^2))/D(i);
    end

    % from Soulsby (pages 70-72) - with fine increments of T, same as
    % common inverse method; coarsen the T interval, and they differ...
    omega(j)=(2*pi)/T(j);
    xi(j)=(omega(j)^2*D(j))/g;
    if xi(j)<=1
        nu(j)=xi(j)^0.5*(1+(0.2*xi(j)));
    else
        nu(j)=xi(j)*(1+(0.2*exp(2-(2*xi(j)))));
    end
    k(j)=nu(j)/D(j);
    lambda(j,4)=(2*pi)/k(j);

end

threshold=ones(size(D))*25;

figure(1)
scrsz=get(0,'ScreenSize');
set(gcf,'Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2])
plot(D,Us,'-x')
hold on
plot(D,threshold,'r:')
xlabel('Depth (m)')
ylabel('Ursell number (U_{s})')

items_period=cell(size(T));
for i=1:length(T)
    items_period{i}=['Period = ',num2str(T(i)),'s'];
end
% items={['Period = ',num2str(T)]};
% legend('Period 2s','Period 3s','Period 4s','Period 5s','Period 6s','Period 7s','Period 8s')
legend(items_period)

figure(2)
scrsz=get(0,'ScreenSize');
set(gcf,'Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2])
plot(T,lambda,'-x');
xlabel('Period (s)')
ylabel('Wavelength (m)')
legend('Common inverse relationshsip','Xu et al. (1999)','Neumann-spectrum derived relationship','Soulsby method')