function [tabindex, index]= readsquidindex(shortphase,folderpath,indexpath,tabindexpath)
tabindex = [];
index = [];

% 
% tabindex
% %tabindex has format:
%{ ,1} filename
%{ ,2} shortfilename
%{ ,3} first index number
%{ ,4} end time(UTC)
%{ ,5} end time (S/C clock)
%{ ,6} number of columns
%{ ,7} number of rows


archivepath = '/Users/frejon/Documents/RosettaArchive/PDS_Archives/DATASETS/SECOND_DELIVERY_VERSIONS/RO-C-RPCLAP-3-1412-CALIB-V0.3';

%s_tabindexfile = sprintf('%s/tabindex/tabindex_%s.mat',indexpath,archiveid);
%fp = fopen(tabindexpath,'r');

%fp = -2;

if exist(tabindexpath)==2

 %   fclose(fp);
    load(tabindexpath);
    'lapdog: successfully loaded server tabindex'
    

substring = strrep(tabindex(1,1),tabindex(1,2),'');
%substring = 'archivepath/2014/MMM/DDD/'

lend = length(substring{1,1})-39-length(shortphase);

tabindexsubstring= substring{1,1}(1:lend); %works for all missionphases

%folderpath= '/Users/frejon/Documents/RosettaArchive/PDS_Archives/DATASETS/SECOND_DELIVERY_VERSIONS/';


tabindex(:,1) = cellfun(@(x) strrep(x,tabindexsubstring,folderpath),tabindex(:,1),'un',0);

    'lapdog: converted server to local tabindex'

%tabindexfile = sprintf('tabindex/tabindex_%s.mat',archiveid);
%save(tabindexfile,'tabindex');

else 
    'error, tabindexfile not found'
end





%s_indexfile = sprintf('server/index/index_%s.mat',archiveid);

if exist(indexpath)==2
    
        
        
    load(indexpath);
    'lapdog: succesfully loaded server index'

    
%substring = strrep(index(1,1),tabindex(1,2),'');
substring = '/data/LAP_ARCHIVE/cronworkfolder/';


%folderpath= '/Users/frejon/Documents/RosettaArchive/PDS_Archives/DATASETS/SECOND_DELIVERY_VERSIONS/';
index2 = index;
'replacing substrings...'
index = struct_string_replace(index2,substring,folderpath); %separate code
%indexfile = sprintf('index/index_%s.mat',archiveid);

%save(indexfile,'index');

else
    'error, indexfile not found'
end
end
