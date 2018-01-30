%% Script calcolo velocità auto by Mosconi 
% Dato un file con coordinate punti mappa, simula andamento auto e calcola
% speed e power.
%
% Idea: fisso un accelerazione costante (somma laterale e longitudinale) e
% costringo l'auto ad avere questa accelerazione.
% QuasiSteadyState model
%
% Errori: accelerazione positiva == accelerazione negativa, mai vero
%   non cosiderata efficienza meccanica (~60%)??
%% init
close all; clear;


%carica file coordinate circuito
%ogni step è di 1mm (LOL)
load percorsocurve.txt
figure; plot(percorsocurve(:,2),percorsocurve(:,1)); %non so orientamento corretto

g=9806;       %G-force 9.806 m/s2 -> 9806 mm/s2
atot= 1.5;    %Accelerazione limite imposta [m/s/s]
print=atot;   %???
massa = 320;  % [Kg]
l=0;
j=1;
dist=0;      %[mm]

%% Compute ragggio curvatura e maxspeed

for i=2:(size(percorsocurve,1)-1)
    l=l+rssq([percorsocurve(i,:)-percorsocurve(i-1,:)]);    %scorro "frammenti"
    if l>100    %100 = 10cm ???
        punti(j,1:2)=percorsocurve(i,1:2);
        punti(j,7)=l;
        dist=dist+l;
        punti(j,11)= dist;   %distanza percorsa
        
        %Calcolo raggio di curvatura (correct!) //possibile implementazione
        %anche con derivata, check wikipedia.
        x1=percorsocurve(i-1,1);
        y1=percorsocurve(i-1,2);
        x2=percorsocurve(i,1);
        y2=percorsocurve(i,2);
        x3=percorsocurve(i+1,1);
        y3=percorsocurve(i+1,2);
        x4 = (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
        y4 = (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
        punti(j,3)=rssq([percorsocurve(i,:)-[x4 y4]]); %raggio curvatura
        
        %saturate per segmenti rettilinei
        if punti(j,3)==inf
            punti(j,3)=1e6;
        end
        
        punti(j,4)=sqrt(atot*g*punti(j,3));     %velocità massima [mm/s]
        punti(j,5)=punti(j,4)*0.0036;           %velocità massima [kmh]
        
        if j>1
            punti(j,6)=(punti(j,4)^2-punti(j-1,4)^2)/(2*atot*g*l); %accelerazione long necessaria
        end
        
        j=j+1;
        l=0;
    end
end

%% controllo frenata
for i=2:size(punti,1)
    for j=i-1:-1:1
%         ay=punti(j+1,4)^2/(punti(j,3)*g);
%         ax=sqrt(atot^2-ay^2);
        vf=punti(j+1,4); %vmax ammessa
        s=punti(j+1,7); %step
        ri=punti(j,3); %rcurvatura
        vi=((g*ri*(2*s*(atot^2*g^2*ri^2 + 4*atot^2*g^2*s^2 - vf^4)^.5 + ri*vf^2))/(g*ri^2 + 4*g*s^2))^(1/2);
%         vi=sqrt(punti(j+1,4)^2+2*ax*g*punti(j,7)); %vi^2=vf^2-2axS
        if vi<punti(j,4) && isreal(vi)
            punti(j,4)=vi;
            punti(j,5)=punti(j,4)*0.0036;
            punti(j,8)=(punti(j+1,4)^2-punti(j,4)^2)/(2*g*punti(j+1,7)); %accelerazione x
            punti(j,9)=punti(j,4)^2/punti(j,3)/g; %accelerazione y
            punti(j,10)=sqrt(punti(j,8)^2+punti(j,9)^2);
        end
    end
end

%% Controllo accelerazione
punti=flip(punti,1);

for i=2:size(punti,1)
    for j=i-1:-1:1
%         ay=punti(j+1,4)^2/(punti(j,3)*g);
%         ax=sqrt(atot^2-ay^2);
        vf=punti(j+1,4);   %vmax ammessa
        s=punti(j+1,7);    %step
        ri=punti(j,3);     %rcurvatura
        vi=((g*ri*(2*s*(atot^2*g^2*ri^2 + 4*atot^2*g^2*s^2 - vf^4)^.5 + ri*vf^2))/(g*ri^2 + 4*g*s^2))^(1/2);
%         vi=sqrt(punti(j+1,4)^2+2*ax*g*punti(j,7)); %vi^2=vf^2-2axS
        if vi<punti(j,4) && isreal(vi)
            punti(j,4)=vi;
            punti(j,5)=punti(j,4)*0.0036;
            punti(j,8)=(punti(j+1,4)^2-punti(j,4)^2)/(2*g*punti(j+1,7)); %accelerazione x
            punti(j,9)=punti(j,4)^2/punti(j,3)/g; %accelerazione y [adimensional, G-normalized]
            punti(j,10)=sqrt(punti(j,8)^2+punti(j,9)^2);
        end
    end
end

punti=flip(punti,1);

%%
time=cumsum(punti(:,7)./punti(:,4));  %vettore tempo. Time = space / speed.

acclong=[1; diff(punti(:,4))]./[1; diff(time)];    %accelerazione longitudinale [mm/s2]

acclat=punti(:,9);                      %accelerazione laterale [G]
acctot=rssq([acclong./g acclat],2);        %accelerazione totale [G] (dovrebbe essere quella imposta da me)
    
percorsovel(:,1)= punti(:,11);  %space [mm]
percorsovel(:,2)= punti(:,4);   %speed [mm/s]

pow=acclong(acclong>0)/1000.*punti(acclong>0,4)/1000*massa;           %[W]
energ=sum(pow.*punti(acclong>0, 7)./punti(acclong>0, 4));             %J
energ_kWh=sum(pow.*punti(acclong>0, 7)./punti(acclong>0,4))/3.6e6;    %kWh


%% Plots
figure; 
subplot(3,1,1);
plot(time,acclat); title('Accelerazione Laterale'); ylabel('G');
subplot(3,1,2);
plot(time,acclong./g); title('Accelerazione Longitudinale'); ylabel('G');
subplot(3,1,3);
plot(time,acctot); title('Accelerazione Totale'); ylabel('G');
xlabel('s');

figure; 
plot(time,punti(:,4)./1000); xlabel('s'); ylabel('m/s'); title('Velocità'); grid on;

figure;
plot(time(acclong>0),pow./1000); ylabel('kW'); xlabel('s'); title(['Pwr [Tot Energy needed: ', num2str(energ_kWh), 'kWh]']);

fprintf(' Distanza percorsa: %f m \n', dist./1000);
fprintf(' Tempo impiegato: %.3f s \n', time(end));
fprintf(' Velocità media: %.3f m/s \n', mean(punti(:,4)./1000));
fprintf(' Accelerazione imposta: %.1f G \n' , atot);
fprintf(' Energy utilizzata: %.4f J || %.4f kWh \n', energ, energ_kWh);

