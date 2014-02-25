clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT PARAMETERS
%!!NOTA IMPORTANTE: tutte le tecniche di Shrinkage devono essere numerate
%successivamente a ShrinkageConstantCorrelation; in questa versione, se si testa un algoritmo di
%Shrinkage si deve necessariamente testare ShrinkageConstantCorrelation.
%
%nel file efficient_frontier decidi se commentare la riga
%                 'correlazioneA(find(correlazioneA<0))=0;'
%che pone a 0 le correlazioni negative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period='y'        %period {'y' 'm' 't' 's' 'b','n','2'}
shortselling='N'  %shortselling {'N','S'}
input='H'         %input {'N','S','H'}   N=NYSE100 (con indice=SP500), S=series_SI (con indice=media), H=series_HNFM (con indice=media)
window='m'        %window {'m','e'}   m=(giorniB=1mese), e=(giorniB=giorniA)
numalgoritmi=10   %numero di algoritmi da testare
algoritmo=cell(numalgoritmi,1);
algoritmo{1}='Markowitz';
algoritmo{2}='SIS_Pw';
algoritmo{3}='RMTCovarianceStanley';
algoritmo{4}='RMTCovarianceBouchaud';
algoritmo{5}='MSTAverage_correlation';
algoritmo{6}='MSTWeighted_correlation';
algoritmo{7}='MSTHausdorff_s';
algoritmo{9}='ShrinkageSIS_Pw';
algoritmo{10}='ShrinkageCommonCoVariance';
algoritmo{8}='ShrinkageConstantCorrelation';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if input=='N'
    load '../data/NYSE100'
    rend=logreturn;
    clear logreturn
    copyfile('../data/SP500index.mat', '../data/index.mat')
else if input=='S'
        load '../data/series_SI'
    else
        if input=='H'
           load '../data/series_HNFM'  
        end
    end
    rend=(series(1:end-1,:)-series(2:end,:))./series(2:end,:); %calcola i rendimenti dalle serie
    series_index=mean(series,2);%calcola l'indice
    logreturn=(series_index(1:end-1,:)-series_index(2:end,:))./series_index(2:end,:); %calcola i rendimenti 
    save(strcat('./','index.mat'), 'logreturn')
    clear series logreturn series_index
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET ADDITIONAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gg=[1,22,41,61,83,104,125,147,168,189,212,231,253,273,292,314,335,355,377,399,420,441,463,483,505,524,543,566,587,607,629,650,672,693,714,735,757,777,796,819,838,860,880,900,921,941,963,984,1004,1025,1044,1066,1086,1108,1129,1150,1173,1188,1211,1232,1252,1273,1291,1311,1333,1355,1375,1397,1419,1439,1462,1482,1503,1524,1543,1564,1585,1606,1627,1649,1670,1691,1714,1733,1755,1775,1794,1817,1838,1858,1879,1900,1922,1942,1964,1985,2007,2027,2046,2068,2089,2110,2132,2152,2175,2196,2217,2238,2259,2279,2298,2321,2340,2362,2384,2404,2427,2447,2469,2490,2510,2530,2549,2571,2591,2613,2634,2655,2678,2697,2720,2741,2761];
years={'1997','1998','1999','2000','2001','2002','2003','2004','2005','2006','2007'};

if period=='y'
    mesi={'year'};
    years=years(1:(end-1));
    NORM=1;
else if period=='m'
        mesi={'January','February','March','April','May','June','July','August','September','October','November','December'};
        NORM=12;
     else if period=='t'
            mesi={'T1','T2','T3','T4'};
            NORM=4;
          else if period=='s'
                mesi={'S1','S2'};
                NORM=2;
               else if period=='b'
                    mesi={'B1','B2','B3','B4','B5','B6'};
                    NORM=6;
                    else if period=='n' 
                        mesi={'nine_months'};
                        years={'N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13'};
                        NORM=12/9;
                        else if period=='2'
                            mesi={'2years'};
                            years={'97-98','98-99','99-00','00-01','01-02','02-03','03-04','04-05'};
                            NORM=0.5;
                            end
                        end
                   end
              end
         end
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(strcat('./',window,'/',input,'/',period,shortselling))
minimumeigenvaluetable=zeros(length(years)*length(mesi)-1,numalgoritmi);
rischiorealizzatotable=zeros(length(years)*length(mesi)-1,numalgoritmi);
rischiopredettotable=zeros(length(years)*length(mesi)-1,numalgoritmi);
[noth,azioni]=size(rend); 
oldminvarweights=zeros(azioni,numalgoritmi);


if (period=='2')
    giorni_iniziali=gg(1:(12):end);
else
giorni_iniziali=gg(1:(12/NORM):end);
end
for h=1:length(years)
    if (period=='y'|| period=='n'|| period=='2')
        else if h==length(years)
            mesi=mesi(1:(end-1));
        end
    end
    for hh=1:length(mesi) 
        if (period=='n' || period=='2')
            hhh=h;
        else
            hhh=(h-1)*NORM+hh;
        end
		giorno_iniziale=giorni_iniziali(hhh);
        if (period=='2')
            giorniA=giorni_iniziali(hhh+2)-giorni_iniziali(hhh);
            if (window=='e')
		    giorniB=giorni_iniziali(hhh+4)-giorni_iniziali(hhh+2);
            else if (window=='m')
            warning('giorniB=giorni_iniziali...........................................................')
            end
            end
        else
		giorniA=giorni_iniziali(hhh+1)-giorni_iniziali(hhh);
		if (window=='e')
            giorniB=giorni_iniziali(hhh+2)-giorni_iniziali(hhh+1); 
        end
        end
        if (window=='m')
        giorno=find(gg==giorni_iniziali(hhh+1))
	    giorniB=gg(giorno+1)-gg(giorno); %un mese
        end
		%OUTPUT
		disp(strcat(years(h),mesi(hh),shortselling));
		minvarweights=zeros(azioni,numalgoritmi);
		for i=1:numalgoritmi
            tic
			disp(strcat('processing_',algoritmo{i}))
			[rischiominvarpredetto,rischiominvarrealizzato,alpha,rischiopercentuale,pesiminvar,minimoautovalore]=efficient_frontier(algoritmo{i},rend,shortselling,giorno_iniziale,giorniA,giorniB);  
			minvarweights(:,i)=pesiminvar;
            minimumeigenvaluetable(hhh,i)=minimoautovalore;
			rischiorealizzatotable(hhh,i)=rischiominvarrealizzato*sqrt(250);%sqrt(giorni_iniziali(hhh+2)-giorni_iniziali(hhh+1));
            rischiopredettotable(hhh,i)=rischiominvarpredetto*sqrt(250);%sqrt(giorni_iniziali(hhh+2)-giorni_iniziali(hhh+1));
			%rischiopercentualetable(hhh,i)=rischiopercentuale;
			%alphatable(hhh,i)=alpha;
		    %Nefftable(hhh,i)=1./(sum(((pesiminvar).^2)')');
		    %turnover(hhh,i)=sum(abs(minvarweights(:,i)-oldminvarweights(:,i)));
            toc
        end
        oldminvarweights=minvarweights;
		clear i 
        name=strcat('pesi',years(h),mesi(hh),shortselling);
        filename=sprintf('%s',name{:});
		save(strcat('./',window,'/',input,'/',period,shortselling,'/',filename),'minvarweights')
		save(strcat('./',window,'/',input,'/',period,shortselling,'/',filename,'.txt'),'minvarweights','-ascii','-double')   
		
    end    
end
save(strcat('./',window,'/',input,'/',period,shortselling,'/minimumeigenvalue_',period,shortselling),'minimumeigenvaluetable') 
save(strcat('./',window,'/',input,'/',period,shortselling,'/minimumeigenvalue_',period,shortselling,'.txt'),'minimumeigenvaluetable','-ascii','-double')  
save(strcat('./',window,'/',input,'/',period,shortselling,'/rischiorealizzato_',period,shortselling),'rischiorealizzatotable') 
save(strcat('./',window,'/',input,'/',period,shortselling,'/rischiopredetto_',period,shortselling),'rischiopredettotable') 
save(strcat('./',window,'/',input,'/',period,shortselling,'/rischiorealizzato_',period,shortselling,'.txt'),'rischiorealizzatotable','-ascii','-double')
save(strcat('./',window,'/',input,'/',period,shortselling,'/rischiopredetto_',period,shortselling,'.txt'),'rischiopredettotable','-ascii','-double')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TABELLA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tabella=cell(numalgoritmi+1,5);
tabella{1,1}=strcat(period,shortselling);
tabella{1,2}='Predicted Risk';
tabella{1,4}='Realized Risk';
L=length(rischiopredettotable);
t2=mean(rischiopredettotable)*100;
t3=sqrt(var(rischiopredettotable)/L)*100;
t4=mean(rischiorealizzatotable)*100;
t5=sqrt(var(rischiorealizzatotable)/L)*100;
for i=1:numalgoritmi
	tabella{i+1,1}=algoritmo{i};
	tabella{i+1,2}=t2(i);
	tabella{i+1,3}=t3(i);
	tabella{i+1,4}=t4(i);
	tabella{i+1,5}=t5(i);
end
save(strcat('./',window,'/',input,'/',period,shortselling,'/tabella_',period,shortselling),'tabella')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CLEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear noth oldminvarweights  t2 t3 t4 t5 pesiminvar rischiominvarpredetto rischiominvarrealizzato name filename gg NORM giorni_iniziali  giorno_iniziale giorniA giorniB tabella

