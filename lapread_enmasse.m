%------------------------------------------------------------
% Author: Fredrik Johansson 09Jan2014 frejon@gmail.com
% The following routine needs:
% readsquidindex.m
% struct_string_replace.m (third party) (can be skipped if user is ignorant
% of (calibrated) index
% lap_import.m (see : https://github.com/frejon/lap_import )
% Description:  reads cached index from LAP archives
%   -tabindex contains resampled files from 'DERIVED' archive,
%   i.e. I*L,V*L,I*S,B*S,I*H,V*H.TAB. 
%   -index contains 'CALIB' archive files
%
%------------------------------------------------------------



%------------------------------------------------------------
%NEED THESE INPUTS strings

shortphase = '1501';
folderpath = '/data/LAP_ARCHIVE/';
%folderpath ='/Users/frejon/Documents/RosettaArchive/PDS_Archives/DATASETS/SECOND_DELIVERY_VERSIONS/';
indexpath = '/usr/local/src/gse/lapdog/index/index_1501_3_v2.mat';
tabindexpath = '/usr/local/src/gse/lapdog/tabindex/tabindex_1501_3.mat';
%tabindexpath = '/Users/frejon/Documents/RosettaArchive/Lapdog_GIT/tabindex/tabindex_1501_3.mat';



%what probe, which file type?
probe = 1;
typeID = sprintf('I%iS',probe); %I1S or I2S

%------------------------------------------------------------




%loop this to get index from all mission phases
[tabindex,index]= readsquidindex(shortphase,folderpath,indexpath,tabindexpath);

%get a list of the types of all files
tabtype = cellfun(@(x) x(end-6:end-4),tabindex(:,2),'un',0);

%tabindex does not contain A*S files, obtained by simple strrep
%(in loop)
A_ID = sprintf('A%iS',probe); %A1S or A2S
B_ID = sprintf('B%iS',probe); %A1S or A2S

%find index number of all files of type typeID.
ind_type= find(strcmp(typeID, tabtype));

clear dataraw
for i=1:length(ind_type);
    
    
    i
   
    rfile =tabindex{ind_type(i),1}; % get filepath
    %rfile = strrep(rfile,typeID,A_ID); %I*S-> A*S
    %if you need voltage bias steps, get it here
    bfile = strrep(rfile,typeID,B_ID); %I*S-> B*S

    if exist(rfile)==2
        
        formatin = 'YYYY-mm-ddTHH:MM:SS';
        
        temp=lap_import(rfile); %read file
        %OPTIONAL: use SPICE to get UTC string
        %temp.t1 = datenum(cspice_et2utc(cspice_str2et(temp.START_TIME_UTC),'ISOC',0),formatin);
        temp2=lap_import(bfile);
        temp.bias_potentials =temp2.bias_potentials;
        
        
        dataraw(i) = temp;
    else
                fprintf(1,'file missing: %s\n', rfile);
    end
    
    clear temp temp2
    
       
end

%everything stored in dataraw


