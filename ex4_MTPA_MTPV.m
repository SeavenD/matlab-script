
close all
clear
clc
linewidth=2;
fontsize=14;
pic = 0;
set(0,'DefaultLineLineWidth',linewidth)
set(0,'DefaultTextFontsize',fontsize)
set(0,'Defaultfigurecolor',[1 1 1])
% 
phim = 0.07831;
Ld = 0.000136;
Lq = 0.000264;



% phim = 0.07228;
% Ld = 0.000134;
% Lq = 0.000369;

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
%         if(abs(torque(i,j))>600)
%             torque(i,j) = NaN;
%         end
%     end
% end 

%% calculate MTPV
% ud = -we*Lq*iq;
% uq = we*(phim+Ld*id)
% (ud^2-uq^2)/we^2 = (Lq*iq)^2 + (phim+Ld*id)^2
% r^2 = (Lq*iq)^2 + (phim+Ld*id)^2
OmegaMtpv = [2000:10:30000,50000,5000000];
% FluxMtpv =flux_ellipse;
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
% ttt = polar*(-0.5*dertL^2*iq.^2+iq*phim);

%% MTPA  lagrlangth 
% te_l = 1:227;
% k_cof1 = 8*polar*phim^2/dertL;
% k = (-k_cof1-sqrt(k_cof1^2 + 16.*te_l.^2))./(2.*te_l);
% id_k = (k.^2*phim/dertL)./(4-k.^2);
% iq_k = ((2*phim/dertL).*k)./(4-k.^2);
% 
%% calc Torque - angle
% ???=?[(?????)???????]???????
% ?es=????????sin
% ?e=?es+?e?
theta =0:0.1:pi;
Angle =theta/pi*180;
Ter = 3/2*polar*(dertL)*Ismax.*cos(theta)* Ismax.*sin(theta);
Tes = 3/2*polar* phim*Ismax.*sin(theta);
Te =Tes +Ter;
figure(1)
plot(Angle,Ter,Angle,Tes,Angle,Te)
h=legend('Ter','Tes','Te');set(h,'Location','Best');
title('Angle-Torque')
xlim([0,180])
xlabel('Angle(Degree)');
ylabel('Torque(Nm)');

%% plot
figure(2)
% plot(id_e,iq_e)
% % plot(id,iq)
% % hold on
% % plot(id,-iq)
% xlabel('id(A)')
% ylabel('iq(A)')
% contour(id0,iq0,adadad',35)
% colorbar
% plot(Torque_mtpa,iq,Torque_mtpa,id)
% hold on
plot(7.5*OmegaMtpv/pi,torqueMtpv,'g')
hold on
% plot([0]Torque_mtpa)
line([0,7.5*omega_few(1)/pi],[torque_few(1),torque_few(1)])
plot(7.5*omega_few/pi,torque_few,'r')
xlim([0,30000])
title('Speed - Torque')
xlabel('Speed(rpm)')
ylabel('Torque(Nm)');



idd = [0,id(2:end)];
iqq = [0,iq(2:end)];
Torque_mtpaa = [0,Torque_mtpa(2:end)];
[Ten_idn,ids] = polyfit(Torque_mtpaa,idd,4);
[Ten_iqn,iqs] = polyfit(Torque_mtpaa,iqq,4);

figure(3)
idpoly =(Ten_idn(1).*Torque_mtpaa.^4 + Ten_idn(2).*Torque_mtpaa.^3 +Ten_idn(3).*Torque_mtpaa.^2 +Ten_idn(4).*Torque_mtpaa + Ten_idn(5));
iqpoly =(Ten_iqn(1).*Torque_mtpaa.^4 + Ten_iqn(2).*Torque_mtpaa.^3 +Ten_iqn(3).*Torque_mtpaa.^2 +Ten_iqn(4).*Torque_mtpaa + Ten_iqn(5));
Torque_mtpa_Real = T_coff*polar*iqq.*(dertL*idd+phim);
Torque_mtpa_Poly = T_coff*polar*iqpoly.*(dertL*idpoly+phim);
subplot(2,2,1)
plot(Torque_mtpaa,id,'b--o',Torque_mtpaa,iq,'g--o',Torque_mtpaa,idpoly,'k',Torque_mtpaa,iqpoly,'r');
h=legend('Te-id','Te-iq','poly(Te-id)','poly(Te-iq)');set(h,'Location','Best');
title('Torque-idq')
xlabel('Torque(Nm)');
ylabel('Current(A)');
subplot(2,2,3)
idpoly(idpoly>=0)=0; %将idpoly中大于等于0的值全部替换为0
plot(Torque_mtpaa,id-idpoly,'b--o',Torque_mtpaa,iq-iqpoly,'g--o');
h=legend('id-idpoly','iq-iqpoly');set(h,'Location','Best');
title('idq-idqpoly')
xlabel('Torque(Nm)');
ylabel('Dert(A)');

subplot(2,2,[2,4])
plot(Torque_mtpaa,Torque_mtpa_Real - Torque_mtpa_Poly,'b--o');
title('Torque Devation between Set-Poly')
xlabel('Torque(Nm)');
ylabel('Dert(Nm)');




figure(4)
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
contour(id0,iq0,torque',20);% same torque line
c = colorbar;
c.Label.String = 'Torque(Nm)';
xlim([-Ismax-50,Ismax+50])
ylim([-Ismax-50,Ismax+50])
hold on
% v = [1000,600,400,300,233,200];
% contour(id0,iq0,omega',v,'-.k','ShowText','on')%  same flux line
contour(id0,iq0,flux_ellipse','-.k','ShowText','on')%  same flux line
hold on
plot(idLMtpv,iqLMtpv,'g');
legend('Ismax','MTPA','Torque(Nm)','Flux','MTPV')
