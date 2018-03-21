%rsync -rzth --progress --delete  --include '*/' --include '*A1S.TAB' --exclude '*' frejon@squid.irfu.se:/data/LAP_ARCHIVE/RO-C*5-1*DERIV*/ /Users/frejon/Documents/RosettaArchive/PDS_Archives/DATASETS/SECOND_DELIVERY_VERSIONS/A1S/
% use bash script above first!
%
%Then run this to index, collect, sort and store information
%Plot as well if you like

kernelFile='/Users/frejon/Rosetta/Lapdog_GIT/metakernel_rosetta.txt';
%cspice_furnsh(kernelFile);



skip =1; %use this if you want to skip the indexing and file reading (long)
if ~skip
    
    %get all data files
    fileList = getAllFiles('/Users/frejon/Rosetta/temp/');
    clear dataraw     % make sure certain dataraw is cleared
    j=1; %counter
    for i=1:length(fileList) %loop all files
               
        i %output counter to prompt
        rfile=fileList{i,1}(1:end);

        %check if file has correct format and exists (to be safe)            
        if exist(rfile)==2 && strcmp(rfile(end-3:end),'.TAB') && str2double(rfile(end-5))==1;
            
            %rfile
            %import analysis file! lap_import.m is on git: https://github.com/frejon/lap_import/
            data_temp=lap_import(rfile);

            %%if you want UTC time, you need mice (google) and correct
            %%rosetta spice kernels
          
            %awkward way of getting UTC time via ephemeris time.
            formatin = 'YYYY-mm-ddTHH:MM:SS';         
            data_temp.t1 = datenum(cspice_et2utc(cspice_str2et(data_temp.START_TIME_UTC),'ISOC',0),formatin);

%          orbit.m does many things, in this case we want the radius from the comet           
%          orbit.m can be found here: https://github.com/frejon/Lapdog_GIT/blob/master/orbit.m

            if length(data_temp.t1) == 2 %orbit will get freaky results in this special case
                [altitude1,junk,junk] = orbit('Rosetta',cspice_str2et(data_temp.START_TIME_UTC(1)),'CHURYUMOV-GERASIMENKO','ECLIPJ2000','preloaded');
                altitude = [altitude1; altitude1];
            else
                
                [altitude,junk,junk] = orbit('Rosetta',cspice_str2et(data_temp.START_TIME_UTC),'CHURYUMOV-GERASIMENKO','ECLIPJ2000','preloaded');
            end
            
            data_temp.r = altitude;
            clear junk altitude
            try
                dataraw(j) = data_temp;  
                j=j+1;
                
            catch err
                
            end
            
        else %if .TABfile doesn't exist
            fprintf(1,'skipping %s...\n',rfile);
        end
       
        
    end
    



%------------------------------------------------------------------------------------------------
%select fields manually in this variable!!!!
%fld={'t1'  'r' 'Qualityfactor' 'Vph_knee' 'Iph0' 'Tph' 'asm_ni_aion' 'asm_ni_v_indep' 'ni_v_indep' 'ni_aion' 'asm_ni_v_dep' 'asm_ne_exp' 'ne_exp' 'asm_ne_linear' 'asm_ne_5eV' 'asm_Vsg', 'ni_v_dep'  'ne_linear' 'ne_5eV' 'Vsg' 'Te_linear' 'Te_exp'  'asm_Te_linear' 'asm_Te_exp' 'asm_v_ion' 'v_ion' 'v_aion' 'asm_Vsc_aion'  'Vsc_aion' 'asm_v_aion' 'asm_Vph_knee' 'Rsquared_linear' 'asm_Rsquared_linear' 'asm_Rsquared_exp' 'Rsquared_exp'};
%fld={'t1'  'r' 'Vbar' 'Vph_knee' 'Iph0' 'Tph' 'ni_v_indep' 'ni_aion' 'ne_exp' 'ni_v_dep'  'ne_linear' 'ne_5eV' 'Vsg' 'Te_linear' 'Te_exp' 'v_ion' 'v_aion' 'asm_Vsc_aion'  'Vsc_aion' 'Rsquared_linear' 'Rsquared_exp'};
%or all (automatically)
fld = fieldnames(dataraw(2));
fld = fld(5:end); %remove datestamps cellarrays, they'll only cause problems. we've got t1 instead.

%fld={'t1' 'r' 'asm_ni_aion','asm_ni_v_indep','ni_v_indep','ni_aion'};
%------------------------------------------------------------------------------------------------

len = length(fld);

A1S=[];
meand1=[];
mediand1=[];

%concenate data, but only the fields we want!
for j=1:length(dataraw)

    if j <2
        
        for k=1:len
            meand1.(sprintf('%s',fld{k,1})) =nanmean([dataraw(j).(sprintf('%s',fld{k,1}))]);
            A1S.(sprintf('%s',fld{k,1})) =[dataraw(j).(sprintf('%s',fld{k,1}))];
            mediand1.(sprintf('%s',fld{k,1})) =nanmedian([dataraw(j).(sprintf('%s',fld{k,1}))]);
 
        end
    else
        
        
        for k=1:len
            meand1.(sprintf('%s',fld{k,1})) = [[meand1.(sprintf('%s',fld{k,1}))];nanmean([dataraw(j).(sprintf('%s',fld{k,1}))])];
            mediand1.(sprintf('%s',fld{k,1})) =[[mediand1.(sprintf('%s',fld{k,1}))];nanmedian([dataraw(j).(sprintf('%s',fld{k,1}))])];
            A1S.(sprintf('%s',fld{k,1})) = [[A1S.(sprintf('%s',fld{k,1}))];[dataraw(j).(sprintf('%s',fld{k,1}))]];
        end
    end
    
    
end
'probe  1 done' %probe 2 not implemented, but it's easy.

  clear dataraw     % make sure certain dataraw is cleared
    j=1; %counter
    for i=1:length(fileList) %loop all files probe 2
               
        i %output counter to prompt
        rfile=fileList{i,1}(1:end)

        %check if file has correct format and exists (to be safe)            
        if exist(rfile)==2 && strcmp(rfile(end-3:end),'.TAB') && str2double(rfile(end-5))==2;
            
            
            %import analysis file! lap_import.m is on git: https://github.com/frejon/lap_import/
            data_temp=lap_import(rfile);

            %%if you want UTC time, you need mice (google) and correct
            %%rosetta spice kernels
          
            %awkward way of getting UTC time via ephemeris time.
            formatin = 'YYYY-mm-ddTHH:MM:SS';         
            data_temp.t1 = datenum(cspice_et2utc(cspice_str2et(data_temp.START_TIME_UTC),'ISOC',0),formatin);

%          orbit.m does many things, in this case we want the radius from the comet           
%          orbit.m can be found here: https://github.com/frejon/Lapdog_GIT/blob/master/orbit.m

            if length(data_temp.t1) == 2 %orbit will get freaky results in this special case
                [altitude1,junk,junk] = orbit('Rosetta',cspice_str2et(data_temp.START_TIME_UTC(1)),'CHURYUMOV-GERASIMENKO','ECLIPJ2000','preloaded');
                altitude = [altitude1; altitude1];
            else
                
                [altitude,junk,junk] = orbit('Rosetta',cspice_str2et(data_temp.START_TIME_UTC),'CHURYUMOV-GERASIMENKO','ECLIPJ2000','preloaded');
            end
            
            data_temp.r = altitude;
            clear junk altitude
            try
                dataraw(j) = data_temp;  
                j=j+1;
                
            catch err
                
            end
            
        else %if .TABfile doesn't exist
            fprintf(1,'skipping %s...\n',rfile);
        end
       
        
    end
    

    
len = length(fld);

A2S=[];
meand2=[];
mediand2=[];

%concenate data, but only the fields we want!
for j=1:length(dataraw)

    if j <2
        
        for k=1:len
            meand2.(sprintf('%s',fld{k,1})) =nanmean([dataraw(j).(sprintf('%s',fld{k,1}))]);
            A2S.(sprintf('%s',fld{k,1})) =[dataraw(j).(sprintf('%s',fld{k,1}))];
            mediand2.(sprintf('%s',fld{k,1})) =nanmedian([dataraw(j).(sprintf('%s',fld{k,1}))]);
 
        end
    else
        
        
        for k=1:len
            meand2.(sprintf('%s',fld{k,1})) = [[meand2.(sprintf('%s',fld{k,1}))];nanmean([dataraw(j).(sprintf('%s',fld{k,1}))])];
            mediand2.(sprintf('%s',fld{k,1})) =[[mediand2.(sprintf('%s',fld{k,1}))];nanmedian([dataraw(j).(sprintf('%s',fld{k,1}))])];
            A2S.(sprintf('%s',fld{k,1})) = [[A2S.(sprintf('%s',fld{k,1}))];[dataraw(j).(sprintf('%s',fld{k,1}))]];
        end
    end
    
    
end

'probe 2 done'

end
%rsquared selectioning!

Rsquared_threshold = 0.93; %taken from out of a hat.
ind = [];
ind.lin = A1S.Rsquared_linear > Rsquared_threshold;
ind.exp = A1S.Rsquared_exp > Rsquared_threshold;
% 
% ind.alin = A1S.asm_Rsquared_linear > Rsquared_threshold;
% ind.aexp = A1S.asm_Rsquared_exp > Rsquared_threshold;
% ind.q = ~A1S.Qualityfactor;


Vthr = -25;
plot(A1S.t1(A1S.Vbar>-Vthr & ind.lin),A1S.Vph_knee(A1S.Vbar>-Vthr  & ind.lin),'ro',A1S.t1((A1S.Vbar>-Vthr  & ind.lin)),A1S.Vbar(A1S.Vbar>-Vthr  & ind.lin),'bo')
ax= gca;


figure(1);
ind15=629;
ind16=1233;
nbins= 200;
h1=histogram(A1S.asm_ne_5eV(1:ind15-1),nbins);
hold on
h2 = histogram(A1S.ne_5eV(1:ind15-1),nbins);
hold off
ax = gca;
legend('asm\_ne_5eV 20140921','ne_5eV 20140921')
hold off
figure(2);
h3=histogram(A1S.asm_ne_5eV(ind15:ind16-1),nbins);
hold on
h4=histogram(A1S.ne_5eV(ind15:ind16-1),nbins);
hold off
ax = gca;
legend('asm\_ne_5eV 20150905','ne_5eV 20150905')
figure(3);
h5=histogram(A1S.asm_ne_5eV(ind16:end),nbins);
hold on
h6=histogram(A1S.ne_5eV(ind16:end),nbins);
hold off
ax = gca;
legend('asm\_ne_5eV 20160309','ne_5eV 20160309')
figure(21);
ind15=629;
ind16=1233;
nbins= 200;
h1=histogram(A1S.asm_ni_v_dep(1:ind15-1),nbins);
hold on
h2 = histogram(A1S.ni_v_dep(1:ind15-1),nbins);
ax = gca;
legend('asm\_ni\_v\_dep 20140921','ni_v_dep 20140921')
hold off
figure(22);
h3=histogram(A1S.asm_ni_v_dep(ind15:ind16-1),nbins);
hold on
h4=histogram(A1S.ni_v_dep(ind15:ind16-1),nbins);
hold off
ax = gca;
legend('asm\_ni\_v\_dep 20150905','ni_v_dep 20150905')
figure(23);
h5=histogram(A1S.asm_ni_v_dep(ind16:end),nbins);
hold on
h6=histogram(A1S.ni_v_dep(ind16:end),nbins);
hold off
ax = gca;
legend('asm\_ni\_v\_dep 20160309','ni_v_dep 20160309')
figure(31);
ind15=629;
ind16=1233;
nbins= 200;
h1=histogram(A1S.ASM_Iph0(1:ind15-1),nbins);
hold on
h2 = histogram(A1S.Iph0(1:ind15-1),nbins);
hold off
ax = gca;
legend('ASM\_Iph0 20140921','Iph0 20140921')
figure(32);
h3=histogram(A1S.ASM_Iph0(ind15:ind16-1),nbins);
hold on
h4=histogram(A1S.Iph0(ind15:ind16-1),nbins);
hold off
ax = gca;
legend('ASM\_Iph0 20150905','Iph0 20150905')
figure(33);
h5=histogram(A1S.ASM_Iph0(ind16:end),nbins);
hold on
h6=histogram(A1S.Iph0(ind16:end),nbins);
hold off
ax = gca;
legend('ASM\_Iph0 20160309','Iph0 20160309')
figure(41);
ind15=629;
ind16=1233;
nbins= 200;
h1=histogram(A1S.asm_Te_exp_belowVknee(1:ind15-1),nbins);
hold on
h2 = histogram(A1S.Te_exp_belowVknee(1:ind15-1),nbins);
hold off
ax = gca;
h1.BinLimits=[0 30];
h2.BinLimits=[0 30];
legend('asm\_Te\_exp\_belowVknee 20140921','Te_exp_belowVknee 20140921')
figure(42);
h3=histogram(A1S.asm_Te_exp_belowVknee(ind15:ind16-1),nbins);
hold on
h4=histogram(A1S.Te_exp_belowVknee(ind15:ind16-1),nbins);
hold off
h4.BinLimits=[0 30];
h3.BinLimits=[0 30];
ax = gca;
%ax.XAxis.Limits = [0 30];
legend('asm\_Te\_exp\_belowVknee 20150905','Te_exp_belowVknee 20150905')
figure(43);
h5=histogram(A1S.asm_Te_exp_belowVknee(ind16:end),nbins);
hold on
h6=histogram(A1S.Te_exp_belowVknee(ind16:end),nbins);
hold off
%h5.BinLimits=[0 30];
%h6.BinLimits=[0 30];
ax = gca;
legend('asm\_Te\_exp\_belowVknee 20160309','Te_exp_belowVknee 20160309')
figure(51);
ind15=629;
ind16=1233;
nbins= 200;
h1=histogram(A1S.asm_v_ion(1:ind15-1),nbins);
hold on
h2 = histogram(A1S.v_ion(1:ind15-1),nbins);
hold off
ax = gca;
legend('asm\_v\_ion 20140921','v_ion 20140921')
figure(52);
h3=histogram(A1S.asm_v_ion(ind15:ind16-1),nbins);
hold on
h4=histogram(A1S.v_ion(ind15:ind16-1),nbins);
hold off
ax = gca;
legend('asm\_v\_ion 20150905','v_ion 20150905')
figure(53);
h5=histogram(A1S.asm_v_ion(ind16:end),nbins);
hold on
h6=histogram(A1S.v_ion(ind16:end),nbins);
hold off
ax = gca;
legend('asm\_v\_ion 20160309','v_ion 20160309')





%
figure(300)

ind.lin = A1S.Rsquared_linear > 0.97;
plot((A1S.t1(ind.lin)),A1S.Iph0(ind.lin),'ro')

set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;


figure(301)
%ind.lin = A1S.Rsquared_linear > 0.97;
plot((A1S.t1(ind.lin)),A1S.Vsg(ind.lin),'ro',A1S.t1(ind.lin),A1S.Vph_knee(ind.lin),'bo',A1S.t1(ind.lin),A1S.Vbar(ind.lin),'go')
ax = gca;


set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;


legend('Vsg','Vph_knee')
%legend('asm\_ni\_aion','ni\_aion','asm\_ni\_v\_indep','ni\_v\_indep')

ax.Title.String=('Rsquared filtered Vsc vs time');
figure(30)
plot(A1S.t1(ind.lin),A1S.Vph_knee(ind.lin)-A1S.Vbar(ind.lin),'ro')
ax = gca;

%
figure(302)

plot(log10(A1S.r(ind.lin)),A1S.Iph0(ind.lin),'ro')
ax = gca;
%set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
%datetick('x',20,'keepticks')
grid on;


figure(303)
ind.exp = A1S.Rsquared_exp > 0.95;

plot(A1S.t1(ind.exp),log10(A1S.ne_exp(ind.exp))./(log10((A1S.r(ind.exp)))),'ro')
ax= gca;

%plot(log10(A1S.r(ind.exp)),A1S.Te_exp(ind.exp),'ro')
set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;



figure(304)

ind.lin = A1S.Rsquared_linear > 0.98;

plot((A1S.t1(ind.lin)),log10(A1S.ni_aion(ind.lin))./log10(A1S.r(ind.lin)),'ro')
ax = gca;

set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;


% % 
% % 
% % mind= ~meand1.Qualityfactor;
% % ind= ~A1S.Qualityfactor;
% % 
% % rmediand1=mediand1(~A1S.Qualityfactor);
% % rmeand1=meand1(~A1S.Qualityfactor);
% 
% 
% figure(440)
% plot(mediand1(mind).t1,log10(mediand1(mind).ni_v_indep),'ro',mediand1(mind).t1,log10(mediand1(mind).ni_aion),'r.',mediand1(mind).t1,log10(mediand1(mind).r),'bo')
% %plot(mediand1.t1,log10(mediand1.ni_v_indep),'ro',mediand1.t1,log10(mediand1.ni_aion),'r.',mediand1.t1,log10(mediand1.r),'bo')
% set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
% datetick('x',20,'keepticks')
% grid on;
% 
% 
% figure(410)
% plot(A1S.t1(ind),log10(A1S.ni_v_indep(ind)),'ro',A1S.t1(ind),log10(A1S.ni_aion(ind)),'r.',A1S.t1(ind),log10(A1S.r(ind)),'bo')
% %plot(A1S.t1,log10(A1S.ni_v_indep),'ro',A1S.t1,log10(A1S.ni_aion),'r.',A1S.t1,log10(A1S.r),'bo')
% set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
% datetick('x',20,'keepticks')
% grid on;




figure(441)


plot(log10(A1S.r),A1S.asm_ni_aion,'blacko',log10(mediand1.r),mediand1.asm_ni_aion,'ro:',log10(mediand1.r),500*1/sqrt(mediand1.r),'bo')
ax = gca;
% 
% plot(A1S.r,A1S.asm_ni_aion,'blacko',mediand1.r,mediand1.asm_ni_aion,'ro:')
% 


%plot(data1.t1,data1.asm_ni_v_indep,'ro',data2.t1,data2.asm_ni_v_indep,'blacko',data1.t1,data1.asm_ne_linear,'bo',data2.t1,data2.asm_ne_linear,'black+');
%axis([meand1.t1(1) meand1.t1(end) 0 10000])
% 
% set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
% datetick('x',20,'keepticks')
grid on;
%legend('probe 1 asm\_ni\_v\_dep','probe 2 asm\_ni\_v\_indep','probe 1 asm\_ne\_5eV','probe 1 asm\_ni\_aion','probe 1 asm\_ne\_exp');
legend('all','median')
%legend('asm\_ni\_aion','ni\_aion','asm\_ni\_v\_indep','ni\_v\_indep')

ax.Title.String=([sprintf('ni asm aion density estimations probe 1 vs log10radius all mission phases 19amu')]);



figure(444)


plot(A1S.t1,A1S.asm_ni_aion,'blacko',mediand1.t1,mediand1.asm_ni_aion,'ro:')
ax= gca;

%plot(data1.t1,data1.asm_ni_v_indep,'ro',data2.t1,data2.asm_ni_v_indep,'blacko',data1.t1,data1.asm_ne_linear,'bo',data2.t1,data2.asm_ne_linear,'black+');
%axis([meand1.t1(1) meand1.t1(end) 0 10000])

set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;
%legend('probe 1 asm\_ni\_v\_dep','probe 2 asm\_ni\_v\_indep','probe 1 asm\_ne\_5eV','probe 1 asm\_ni\_aion','probe 1 asm\_ne\_exp');
legend('all','median')
%legend('asm\_ni\_aion','ni\_aion','asm\_ni\_v\_indep','ni\_v\_indep')

ax.Title.String=([sprintf('ni aion density estimations probe 1 all mission phases 19amu')]);


figure(429)
plot(A1S.t1,log10(A1S.ni_aion),'blacko',A1S.t1,log10(A1S.asm_ni_aion),'black+',A1S.t1,log10(A1S.ni_v_indep),'ro',A1S.t1,log10(A1S.ne_linear),'bo',A1S.t1,log10(A1S.asm_ni_v_indep),'go',A1S.t1,log10(A1S.ne_exp),'yellowo',A1S.t1,log10(A1S.asm_ne_linear),'b+')
ax= gca;

set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;
legend('ni\_aion','asm\_ni\_aion','ni\_v\_indep','ne\_linear','asm\_ni\_v\_indep','ni\_v\_indep','ne\exp','asm\_ne\_linear')

ax.Title.String=([sprintf('median Density estimations probe 1 allmission phases 19amu')]);


figure(446)


%plot(A1S.t1,A1S.asm_ni_aion,'blacko',mediand1.t1,mediand1.asm_ni_aion,'ro:')
plot(mediand1.t1,mediand1.asm_ni_aion,'ro',mediand1.t1,mediand1.ni_aion,'bo',mediand1.t1,mediand1.asm_ni_v_indep,'go',mediand1.t1,mediand1.ni_v_indep,'blacko');
ax= gca;

%plot(data1.t1,data1.asm_ni_v_indep,'ro',data2.t1,data2.asm_ni_v_indep,'blacko',data1.t1,data1.asm_ne_linear,'bo',data2.t1,data2.asm_ne_linear,'black+');
%axis([meand1.t1(1) meand1.t1(end) 0 10000])

set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;
%legend('probe 1 asm\_ni\_v\_dep','probe 2 asm\_ni\_v\_indep','probe 1 asm\_ne\_5eV','probe 1 asm\_ni\_aion','probe 1 asm\_ne\_exp');
%legend('all','median')
legend('asm\_ni\_aion','ni\_aion','asm\_ni\_v\_indep','ni\_v\_indep')

ax.Title.String=([sprintf('median Density estimations probe 1 allmission phases 19amu')])



% 
 figure(445)
% 

plot(A1S.t1,A1S.asm_ni_aion,'blacko',A1S.t1,A1S.ni_aion,'ro')
%plot(data1.t1,data1.asm_ni_v_indep,'ro',data2.t1,data2.asm_ni_v_indep,'blacko',data1.t1,data1.asm_ne_linear,'bo',data2.t1,data2.asm_ne_linear,'black+');
%axis([meand1.t1(1) meand1.t1(end) 0 10000])
ax= gca;

set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
datetick('x',20,'keepticks')
grid on;
%legend('probe 1 asm\_ni\_v\_dep','probe 2 asm\_ni\_v\_indep','probe 1 asm\_ne\_5eV','probe 1 asm\_ni\_aion','probe 1 asm\_ne\_exp');
%legend('all','median')
legend('all','median','ni\_aion','asm\_ni\_v\_indep','ni\_v\_indep')

ax.Title.String=([sprintf('asm aion Density estimations probe 1 allmission phases 19amu')])

