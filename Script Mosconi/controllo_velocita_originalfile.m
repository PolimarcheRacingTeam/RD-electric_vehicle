close all
clear
% p1 = [10 5];
% p2 = [6 7];
% x = 1:.1:10;
% m = (p1(1)-p2(1))/(p2(2)-p1(2));
% q = (p2(2)+p1(2))/2-(p2(1)+p1(1))*m/2;
% q = (p2(2)^2-p1(2)^2+p2(1)^2-p1(1)^2)/(2*(p2(2)-p1(2)));
% y = m*x+q;
% axis([-10 10 -10 10])
% hold on
% plot(x,y);
% plot(p1(1),p1(2),'o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','y','MarkerSize',5)
% plot(p2(1),p2(2),'o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','y','MarkerSize',5)
% hold off

% (y2-y1)y+(x2-x1)x+(x1^2+y1^2-(x2^2+y2)^2)/2=0,(y2-y3)y+(x2-x3)x+(x3^2+y3^2-(x2^2+y2)^2)/2=0

%syms x y x1 x2 y1 y2 x3 y3

%(y2-y1)y+(x2-x1)x+(x1^2+y1^2-(x2^2+y2)^2)/2=0,(y2-y3)y+(x2-x3)x+(x3^2+y3^2-(x2^2+y2)^2)/2=0
%[x4,y4] = solve([(y2-y1)*y+(x2-x1)*x+(x1^2+y1^2-(x2^2+y2)^2)/2==0,(y2-y3)*y+(x2-x3)*x+(x3^2+y3^2-(x2^2+y2)^2)/2==0],[x,y])
%s = solve([y == (x1-x2)/(y2-y1)*x+(y2^2-y1^2+x2^2-x1^2)/(2*(y2-y1)),y == (x3-x2)/(y2-y3)*x+(y2^2-y3^2+x2^2-x3^2)/(2*(y2-y3))],[x,y])

% x1 = 1;
% y1 = 4;
% x2 =5;
% y2 =6;
% x3 =6;
% y3 =10;
% x4 = (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
% y4 = (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
% axis([0 10 0 10])
% x = 1:.1:10;
% y = (x1-x2)/(y2-y1)*x+(y2^2-y1^2+x2^2-x1^2)/(2*(y2-y1));
% hold on
% plot(x,y)
% y = (x3-x2)/(y2-y3)*x+(y2^2-y3^2+x2^2-x3^2)/(2*(y2-y3));
% plot(x,y)
% plot(x1,y1,'o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','r','MarkerSize',5)
% plot(x2,y2,'o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','g','MarkerSize',5)
% plot(x3,y3,'o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','b','MarkerSize',5)
% plot(x4,y4,'o','MarkerEdgeColor','k','MarkerSize',5,'MarkerFaceColor','b','MarkerSize',5)
% hold off

load percorsocurve.txt
% plot(percorsocurve(:,2),percorsocurve(:,1))
g=9806;
atot=1.5;     %ACC LIMITE  
print=atot;
massa=300;
l=0;
j=1;
dist=0;
for i=2:(size(percorsocurve,1)-1)
    l=l+rssq([percorsocurve(i,:)-percorsocurve(i-1,:)]);
    if l>100
        punti(j,1:2)=percorsocurve(i,1:2);
        punti(j,7)=l;
        dist=dist+l;
        punti(j,11)=dist;
        x1=percorsocurve(i-1,1);
        y1=percorsocurve(i-1,2);
        x2=percorsocurve(i,1);
        y2=percorsocurve(i,2);
        x3=percorsocurve(i+1,1);
        y3=percorsocurve(i+1,2);
        x4 = (x1^2*y2 - x1^2*y3 - x2^2*y1 + x2^2*y3 + x3^2*y1 - x3^2*y2 + y1^2*y2 - y1^2*y3 - y1*y2^2 + y1*y3^2 + y2^2*y3 - y2*y3^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
        y4 = (- x1^2*x2 + x1^2*x3 + x1*x2^2 - x1*x3^2 + x1*y2^2 - x1*y3^2 - x2^2*x3 + x2*x3^2 - x2*y1^2 + x2*y3^2 + x3*y1^2 - x3*y2^2)/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
        punti(j,3)=rssq([percorsocurve(i,:)-[x4 y4]]); %raggio curvatura
        if punti(j,3)==inf
            punti(j,3)=1e6;
        end
        punti(j,4)=sqrt(atot*g*punti(j,3)); %velocità massima
        punti(j,5)=punti(j,4)*0.0036; %velocità massima kmh
        if j>1;
            punti(j,6)=(punti(j,4)^2-punti(j-1,4)^2)/(2*atot*g*l); %accelerazione long necessaria
        end
        j=j+1;
        l=0;
    end
end

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


punti2=flip(punti,1);

for i=2:size(punti2,1)
    for j=i-1:-1:1
%         ay=punti(j+1,4)^2/(punti(j,3)*g);
%         ax=sqrt(atot^2-ay^2);
        vf=punti2(j+1,4); %vmax ammessa
        s=punti2(j+1,7); %step
        ri=punti2(j,3); %rcurvatura
        vi=((g*ri*(2*s*(atot^2*g^2*ri^2 + 4*atot^2*g^2*s^2 - vf^4)^.5 + ri*vf^2))/(g*ri^2 + 4*g*s^2))^(1/2);
%         vi=sqrt(punti(j+1,4)^2+2*ax*g*punti(j,7)); %vi^2=vf^2-2axS
        if vi<punti2(j,4) && isreal(vi)
            punti2(j,4)=vi;
            punti2(j,5)=punti2(j,4)*0.0036;
            punti2(j,8)=(punti2(j+1,4)^2-punti2(j,4)^2)/(2*g*punti2(j+1,7)); %accelerazione x
            punti2(j,9)=punti2(j,4)^2/punti2(j,3)/g; %accelerazione y
            punti2(j,10)=sqrt(punti2(j,8)^2+punti2(j,9)^2);
        end
    end
end

punti2=flip(punti2,1);

% for i=1:(size(punti,1)-1)
%     for j=i:(size(punti,1)-1)
%         vi=punti(j,4);
%         s=punti(j,7);
%         rf=punti(j+1,3);
%         vf=((g*rf*(2*s*(atot^2*g^2*rf^2 + 4*atot^2*g^2*s^2 - vi^4)^(1/2) + rf*vi^2))/(g*rf^2 + 4*g*s^2))^(1/2);
%         if vf<punti(j+1,4) && isreal(vf)
%             punti(j+1,4)=vf;
%             punti(j+1,5)=punti(j+1,4)*0.0036;
%             punti(j+1,12)=(punti(j,4)^2-punti(j+1,4)^2)/(2*g*punti(j,7)); %accelerazione x
%             punti(j+1,13)=punti(j+1,4)^2/punti(j+1,3)/g; %accelerazione y
%             punti(j+1,14)=sqrt(punti(j+1,12)^2+punti(j+1,13)^2);
%         end
%     end
% end


time=cumsum(punti2(:,7)./punti2(:,4));

acclong=[1; diff(punti2(:,4))]./[1; diff(time)];

acclat=punti2(:,9);
acctot=rssq([acclong acclat],2);

percorsovel(:,1)=punti(:,11);
percorsovel(:,2)=punti(:,4);
% for i=1:size(percorsovel,1)
%     if percorsovel(i,2)==inf
%         percorsovel(i,2)=1e5;
%     end
% end

% accpos=acclong(acclong>0);
pow=acclong(acclong>0)/1000.*punti2(acclong>0,4)/1000*massa;
energ=sum(pow.*punti2(acclong>0,7)./punti2(acclong>0,4)); %J
energ_kWh=sum(pow.*punti2(acclong>0,7)./punti2(acclong>0,4))/3.6e6; %kwh

save(sprintf('percorsovel_%.1f_g.txt', atot),'percorsovel','-ASCII')
% percorsoacc(:,1)=punti(:,11);
% percorsoacc(:,2)=punti(:,
% syms vf vi ri g ay ax atot s
% 
% sol=solve([vi^2/ri==g*ay,ay^2+ax^2==atot^2,vi^2==vf^2+(2*s*g*ax)],[vi ax ay])