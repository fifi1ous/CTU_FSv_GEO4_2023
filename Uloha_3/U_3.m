clc; clear; format long g
RAD=pi/200;

B=[35.970234345864 74.0761702269699];
SS=[1,0,0;
    2,0,300;
    3,300,300
    4,300,0];
% pu=1/1000*RAD; pd=5/1000; Sxy=15/1000;
pu=1*RAD; pd=5; Sxy=15;
%% 
figure(1)
plot(-SS(:,2),-SS(:,3),'Marker','x','LineStyle','none','Color','k');
hold on
plot(-B(1),-B(2),'Marker','x','Color','r');
xlim([-320,20]); ylim([-320,20]);
text(-SS(:,2)+5,-SS(:,3)+5,num2str(SS(:,1)))
axis equal
saveas(gcf,'figure.fig');
hold off
close(gcf);
%% Rajón
% Výpočet
[PR1,PR2,EX1_kov,EX2_kov] = presnost_rajon(pu,pd,Sxy,B,SS);

MIN=find(PR2(:,5)==min(PR2(:,5))); MIN=MIN(1);                            %index minimální smerodatná odchylka
MAX=find(PR2(:,5)==max(PR2(:,5))); MAX=MAX(1);                            %index maximální směrodatná odchylka
PR1_min=PR1(MIN,:);          PR1_max=PR1(MAX,:);
PR2_min=PR2(MIN,:);          PR2_max=PR2(MAX,:);
PR_min=[PR1_min;PR2_min];   PR_max=[PR1_max;PR2_max];

R_MIN=PR1(MIN,1); S_MIN=PR1(MIN,2); R_MAX=PR1(MAX,1); S_MAX=PR1(MAX,2);
% Kovarianční matice pro body s minimální souřadnicovou směrodatnou oschylku 
KOV_mat_PR1_MIN=EX1_kov(2*R_MIN-1:2*R_MIN,2*S_MIN-1:2*S_MIN);
KOV_mat_PR2_MIN=EX2_kov(2*R_MIN-1:2*R_MIN,2*S_MIN-1:2*S_MIN);

% Kovarianční matice pro minimální směrodatnou oschylku 
KOV_mat_PR1_MAX=EX1_kov(2*R_MAX-1:2*R_MAX,2*S_MAX-1:2*S_MAX);
KOV_mat_PR2_MAX=EX2_kov(2*R_MAX-1:2*R_MAX,2*S_MAX-1:2*S_MAX);

%Grafické znázornění
%min
PR_min3=PR_min(:,end-1:end)*2.5;
PR_min2=PR_min;
PR_min2(:,end-1:end)=[];
PR_min2=[PR_min2,PR_min3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PR_min2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y(1,:),-X(1,:),'Color','#00E5EE');
plot(-Y(2,:),-X(2,:),'Color','#EE00EE');
plot(-SSA(:,1),-SSA(:,2),'Color','k','LineWidth',1);
plot(-SSC(:,1),-SSC(:,2),'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Elipsa chyb s vlivu podkladu [mm]','Směr na určovaný bod','Směr na orientaci','Location','southoutside')
title('Grafické znázornění přesnosti rajónu pro minimální Sxy (elipsy 2.5x zvětšeny)')
xlabel('Y[m]')
ylabel('X[m]')
hold off
%max
PR_max3=PR_max(:,end-1:end)*2.5;
PR_max2=PR_max;
PR_max2(:,end-1:end)=[];
PR_max2=[PR_max2,PR_max3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PR_max2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y(1,:),-X(1,:),'Color','#00E5EE');
plot(-Y(2,:),-X(2,:),'Color','#EE00EE');
plot(-SSA(:,1),-SSA(:,2),'Color','k','LineWidth',1);
plot(-SSC(:,1),-SSC(:,2),'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Elipsa chyb s vlivu podkladu [mm]','Směr na určovaný bod','Směr na orientaci','Location','southoutside')
title('Grafické znázornění přesnosti rajónu pro maximální Sxy (elipsy 2.5x zvětšeny)')
xlabel('Y[m]')
ylabel('X[m]')
hold off

%% Rajón zpět
[PRZ,EX1_kov] = presnost_rajon_zpet(pu,pd,B,SS);

MIN=find(PRZ(:,5)==min(PRZ(:,5))); MIN=MIN(1);                             %index minimální smerodatná odchylka
MAX=find(PRZ(:,5)==max(PRZ(:,5))); MAX=MAX(1);                             %index maximální směrodatná odchylka
PRZ_min=PRZ(MIN,:);          PRZ_max=PRZ(MAX,:);

R_MIN=PRZ(MIN,1); S_MIN=PRZ(MIN,2); R_MAX=PRZ(MAX,1); S_MAX=PRZ(MAX,2);

% Kovarianční matice pro body s minimální souřadnicovou směrodatnou oschylku 
KOV_mat_PRZ_MIN=EX1_kov(2*R_MIN-1:2*R_MIN,2*S_MIN-1:2*S_MIN);

% Kovarianční matice pro minimální směrodatnou oschylku 
KOV_mat_PRZ_MAX=EX1_kov(2*R_MAX-1:2*R_MAX,2*S_MAX-1:2*S_MAX);

% Grafické znázornění
%min
PRZ_min3=PRZ_min(:,end-1:end)*10;
PRZ_min2=PRZ_min;
PRZ_min2(:,end-1:end)=[];
PRZ_min2=[PRZ_min2,PRZ_min3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PRZ_min2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot(-SSA(:,1),-SSA(:,2),'Color','k','LineWidth',1);
plot(-SSB(:,1),-SSB(:,2),'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Směr i délka na orientaci','Směr na orientaci','Location','southoutside')
title('Grafické znázornění přesnosti rajónu zpět pro minimální Sxy (elipsa 10x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off
%max
PRZ_max3=PRZ_max(:,end-1:end)*2;
PRZ_max2=PRZ_max;
PRZ_max2(:,end-1:end)=[];
PRZ_max2=[PRZ_max2,PRZ_max3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PRZ_max2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot(-SSA(:,1),-SSA(:,2),'Color','k','LineWidth',1);
plot(-SSB(:,1),-SSB(:,2),'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Směr i délka na orientaci','Směr na orientaci','Location','southoutside')
title('Grafické znázornění přesnosti rajónu zpět pro maximální Sxy (elipsa 2x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off

%% Protínání z délek
[PPD,EX1_kov] = presnost_protina_z_delek(pd,B,SS);

MIN=find(PPD(:,5)==min(PPD(:,5))); MIN=MIN(1);                           %index minimální smerodatná odchylka
MAX=find(PPD(:,5)==max(PPD(:,5))); MAX=MAX(1);                           %index maximální směrodatná odchylka
PPD_min=PPD(MIN,:);          PPD_max=PPD(MAX,:);

R_MIN=PPD(MIN,1); S_MIN=PPD(MIN,2); R_MAX=PPD(MAX,1); S_MAX=PPD(MAX,2);

% Kovarianční matice pro body s minimální souřadnicovou směrodatnou oschylku 
KOV_mat_PPD_MIN=EX1_kov(2*R_MIN-1:2*R_MIN,2*S_MIN-1:2*S_MIN);

% Kovarianční matice pro minimální směrodatnou oschylku 
KOV_mat_PPD_MAX=EX1_kov(2*R_MAX-1:2*R_MAX,2*S_MAX-1:2*S_MAX);


% Grafické znázornění
%min
PPD_min3=PPD_min(:,end-1:end)*10;
PPD_min2=PPD_min;
PPD_min2(:,end-1:end)=[];
PPD_min2=[PPD_min2,PPD_min3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PPD_min2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot([-SSA(:,1);-SSB(:,1)],[-SSA(:,2);-SSB(:,2)],'Color','k','LineWidth',1);
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Měřená délka na určovaný bod','Location','southoutside')
title('Grafické znázornění přesnosti protínání z délek pro minimální Sxy (elipsa 10x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off
%max
PPD_max3=PPD_max(:,end-1:end)*5;
PPD_max2=PPD_max;
PPD_max2(:,end-1:end)=[];
PPD_max2=[PPD_max2,PPD_max3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PPD_max2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot([-SSA(:,1);-SSB(:,1)],[-SSA(:,2);-SSB(:,2)],'Color','k','LineWidth',1);
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Měřená délka','Location','southoutside')
title('Grafické znázornění přesnosti protínání z délek pro maximální Sxy (elipsa 5x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off

%% Protínání vpřed z úhlů
[PPUV,EX1_kov] = presnost_protina_z_uhlu_vpred(pu,B,SS);

MIN=find(PPUV(:,5)==min(PPUV(:,5))); MIN=MIN(1);                           %index minimální smerodatná odchylka
MAX=find(PPUV(:,5)==max(PPUV(:,5))); MAX=MAX(1);                           %index maximální směrodatná odchylka
PPUV_min=PPUV(MIN,:);          PPUV_max=PPUV(MAX,:);

R_MIN=PPUV(MIN,1); S_MIN=PPUV(MIN,2); R_MAX=PPUV(MAX,1); S_MAX=PPUV(MAX,2);

% Kovarianční matice pro body s minimální souřadnicovou směrodatnou oschylku 
KOV_mat_PPUV_MIN=EX1_kov(2*R_MIN-1:2*R_MIN,2*S_MIN-1:2*S_MIN);

% Kovarianční matice pro minimální směrodatnou oschylku
KOV_mat_PPUV_MAX=EX1_kov(2*R_MAX-1:2*R_MAX,2*S_MAX-1:2*S_MAX);



%Grafické znázornění
%min
PPUV_min3=PPUV_min(:,end-1:end)*10;
PPUV_min2=PPUV_min;
PPUV_min2(:,end-1:end)=[];
PPUV_min2=[PPUV_min2,PPUV_min3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PPUV_min2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot([-SSA(:,1);-SSB(:,1)],[-SSA(:,2);-SSB(:,2)],'Color','k','LineWidth',1);
plot(-SSC(:,1),-SSC(:,2),'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Měřené směry na bod','Směr na známý bod','Location','southoutside')
title('Grafické znázornění přesnosti protínání z úhlů vpřed pro minimální Sxy (elipsa 10x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off
%max
PPUV_max3=PPUV_max(:,end-1:end)*5;
PPUV_max2=PPUV_max;
PPUV_max2(:,end-1:end)=[];
PPUV_max2=[PPUV_max2,PPUV_max3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PPUV_max2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot([-SSA(:,1);-SSB(:,1)],[-SSA(:,2);-SSB(:,2)],'Color','k','LineWidth',1);
plot(-SSC(:,1),-SSC(:,2),'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Měřené směry na bod','Směr na známý bod','Location','southoutside')
title('Grafické znázornění přesnosti protínání z úhlů vpřed pro maximální Sxy (elipsa 5x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off

%% Protínání zpět 
[PZ,EX1_kov] = presnost_protinani_zpet(pu,B,SS);

MIN=find(PZ(:,6)==min(PZ(:,6))); MIN=MIN(1);                        %index minimální smerodatná odchylka
MAX=find(PZ(:,6)==max(PZ(:,6))); MAX=MAX(1);                        %index maximální směrodatná odchylka
PZ_min=PZ(MIN(1),:);          PZ_max=PZ(MAX,:);


% Kovarianční matice pro body s minimální souřadnicovou směrodatnou oschylku 
KOV_mat_PZ_MIN=EX1_kov(MIN*2-1:MIN*2,:);
% % Kovarianční matice pro minimální směrodatnou oschylku 
KOV_mat_PZ_MAX=EX1_kov(MAX*2-1:MAX*2,:);


%Grafické znázornění
%min
PZ_min3=PZ_min(:,end-1:end)*10;
PZ_min2=PZ_min;
PZ_min2(:,end-1:end)=[];
PZ_min2=[PZ_min2,PZ_min3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PZ_min2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot([-SSA(:,1);-SSB(:,1);-SSC(:,1)],[-SSA(:,2);-SSB(:,2);-SSC(:,2)],'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Měřené směry na orientace','Location','southoutside')
title('Grafické znázornění přesnosti protínání zpět pro minimální Sxy (elipsa 10x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off
%max
PZ_max3=PZ_max(:,end-1:end)*2;
PZ_max2=PZ_max;
PZ_max2(:,end-1:end)=[];
PZ_max2=[PZ_max2,PZ_max3];
[X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(PZ_max2,SS,B);
openfig('figure.fig') ;
hold on 
plot(-Y,-X,'Color','#00E5EE');
plot([-SSA(:,1);-SSB(:,1);-SSC(:,1)],[-SSA(:,2);-SSB(:,2);-SSC(:,2)],'Color','k','LineWidth',1,'LineStyle','--');
axis equal
legend('Body sítě','Určovaný bod','Elipsa chyb bez vlivu podkladu [mm]','Měřené směry na orientace','Location','southoutside')
title('Grafické znázornění přesnosti protínání zpět pro maximální Sxy (elipsa 2x zvětšena)')
xlabel('Y[m]')
ylabel('X[m]')
hold off