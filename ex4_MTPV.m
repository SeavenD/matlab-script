
close all
clear
clc
linewidth=2;
fontsize=14;
pic = 0;
set(0,'DefaultLineLineWidth',linewidth)
set(0,'DefaultTextFontsize',fontsize)
set(0,'Defaultfigurecolor',[1 1 1])

phim = 0.07831;
Ld = 0.000136;
% Ld = 0.000263;
Lq = 0.000264;
dertL = Ld-Lq;
polar = 4;
Ismax = 800;
T_coff =3/2;

%% calculate Ismax circle
angle = pi/2.1:0.1:3.1*pi/2;
id_c = Ismax*cos(angle);
iq_c = Ismax*sin(angle);

%% calculate Flux&torque ellipse

id0 = 50:-1:-Ismax;
iq0 = -Ismax:Ismax;
for i=1:length(id0)
    idphase = (id0(i)*Ld + phim)^2;
    for j=1:length(iq0)
        iqphase = iq0(j)^2 * Lq^2;
        flux_ellipse (i,j) = sqrt(idphase+iqphase);   
        omega (i,j) = (340/1.732)./(2*pi*flux_ellipse (i,j));% Hz
        torque(i,j)= T_coff*polar*(dertL*id0(i) + phim).*iq0(j);% calculate torque
        if(abs(torque(i,j))>600)
            torque(i,j) = NaN;
        end
    end
end    
iq_few = 672:-1:254;
id_few = -sqrt(Ismax^2-iq_few.^2);
% for i=1:length(id_few)
    idphase_few = (id_few*Ld + phim).^2;
%     for j=1:length(iq_few)
        iqphase_few = iq_few.^2 * Lq^2;
        flux_ellipse_few  = sqrt(idphase_few +iqphase_few);   
        omega_few  = (340/1.732)./(flux_ellipse_few );
        torque_few = T_coff*polar*(dertL*id_few + phim).*iq_few;% calculate torque

%% calculate MTPV

OmegaMtpv = [2000:10:30000,50000,5000000];
FluxMtpv = (340/1.732)./OmegaMtpv;
sinthetaMtpv_num = -Lq*phim+sqrt((Lq*phim)^2 + 8*(dertL*FluxMtpv).^2);
sinthetaMtpv_dem = 4*dertL*FluxMtpv;
sinthetaMtpv = sinthetaMtpv_num./sinthetaMtpv_dem;
costhetaMtpv = sqrt(1-sinthetaMtpv.^2);
iqMtpv = FluxMtpv.*costhetaMtpv/Lq;
idMtpv = FluxMtpv.*sinthetaMtpv/Ld - phim/Ld;
for ii =1:length(iqMtpv)
    if((idMtpv(ii)^2 + iqMtpv(ii)^2)>Ismax^2)
        idMtpv(ii) = NaN;
        iqMtpv(ii) = NaN;
    end
end
torqueMtpv= T_coff*polar*(dertL.*idMtpv + phim).*iqMtpv;% calculate torque
idLMtpv =[idMtpv,idMtpv];
iqLMtpv =[iqMtpv,-iqMtpv];
%% calculate MTPA
Is =0:1:Ismax;
if(abs(dertL)<0.0001) % SPM
    cos_beta = 0.0;
else                  %IPM
    pp = sqrt(phim^2+8*dertL^2.*(Is.^2));
    cos_beta = (-phim+pp)./(4*(dertL.*Is));
end

id =Is.*cos_beta;
iq = sqrt(Is.^2 - id.^2);
id_L = [id,id];
iq_L = [iq,-iq];
Torque_mtpa = T_coff*polar*iq.*(dertL*id+phim);

%% plot
figure(1)
plot(id_c,iq_c,'r') % Ismax  circle
axis equal
xlabel('id(A)')
ylabel('iq(A)')
xlim([-Ismax-50,Ismax+50])
ylim([-Ismax-50,Ismax+50])

hold on
plot(id_L,iq_L,'b');% MTPA
axis equal

hold on
contour(id0,iq0,torque',30);
c = colorbar;
c.Label.String = 'Torque(Nm)';
xlim([-Ismax-50,Ismax+50])
ylim([-Ismax-50,Ismax+50])
hold on
% v = [0.068,0.0580,0.0485,0.0377,0.0272,0.01695,0.0055];
% contour(id0,iq0,flux_ellipse',v,'-.k','ShowText','on')%  same flux line
contour(id0,iq0,flux_ellipse','-.k','ShowText','on')%  same flux line
hold on
plot(idLMtpv,iqLMtpv,'g');
legend('Ismax','MTPA','Torque(Nm)','Flux','MTPV')
