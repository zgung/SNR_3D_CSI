%#####################################################
%### converts a CSI Dicom file to FID text files   ###
%#####################################################
function [] = fidextract(directory)

% edited by minarikova.lenka@gmail.com
% %%%%%%%%%%%%%%% use:

% set the path to folder with 3D PRESS MRS SIEMENS dicoms
% the script will make a new folder inthere, named "Output"
% with the appropriate txt files with single spectra for jmrui

%clear all;
%directory = '~/Desktop/Current/Input/';

cd(directory);
disp(strcat('Processing:',directory));
%% <<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>
% read from spectroscopic file from the text part
% beginnig by "### ASCCONV BEGIN ###" and ending by "### ASCCONV END ###"
patient = dir(pwd);
for k = length(patient):-1:1 % find .IMA file
    fname = patient(k).name;
    if fname(1) == '.'
        patient(k) = [];
    elseif strcat(fname(end-2),fname(end-1),fname(end)) ~= strcat('IMA')
        patient(k) = [];
    end
end
fold_n = 1;

for aaa = 1:length(patient)

    voxel = read_ascconv_lenk(patient(aaa).name); % read parrameter from ascii part of the dicom
    nfo1 = dicominfo(patient(aaa).name); %read parameters from dicom header
    if isfield(nfo1.PatientName,'GivenName') == 0
        nfo1.PatientName.GivenName = strcat('noname');
    end
    fn = fullfile('Output',strcat(nfo1.PatientName.FamilyName,'_',nfo1.PatientName.GivenName));
    if exist(fn, 'dir')
        warning(.... creating a new folder);
        mkdir('Output',strcat(nfo1.PatientName.FamilyName,'_',nfo1.PatientName.GivenName,'_',fold_n));
        fold_n = fold_n + 1;
    else
        mkdir(fn)
    end

    outdir = strcat(directory,'Output/',strcat(nfo1.PatientName.FamilyName,'_',nfo1.PatientName.GivenName),'/');
    voxel.size_z = voxel.FoV_z / voxel.number_z; % the size of 1 voxel after zero filling in mm
    voxel.size_y = voxel.FoV_y / voxel.number_y;
    voxel.size_x = voxel.FoV_x / voxel.number_x;

    % fid = fopen('/tmp/Parameter.txt','r');
    % Para = fscanf(fid,'%i %i %i %i %f %f %f %f %f',12);
    % dicomdir = fscanf(fid,'%s',1);
    % outdir = fscanf(fid,'%s',1);
    % pat_name = fscanf(fid,'%s',1);

    % fclose(fid);
    dicomdir = patient(aaa).name;
    pat_name = strcat(nfo1.PatientName.FamilyName,'_',nfo1.PatientName.GivenName);
    COL= voxel.number_x;%Para(1);
    ROW= voxel.number_y;%Para(2);
    THK= voxel.number_z;%Para(3);
    vecSize = voxel.vecSize;%Para(4);
    samInt  = voxel.sampint;%Para(5);
    freq  = voxel.freq;%Para(6);
    COL_P= (voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x) - 1) * 2; %include only voxels inside the PRESS box
    ROW_P = (voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y) - 1) * 2;
    THK_P = (voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z) - 1) * 2;

    %%
    fid = fopen(dicomdir,'r');

    fseek(fid, -((vecSize*ROW*COL*THK*2*4)), 'eof');
    Imag = fread(fid,'float32');
    fclose(fid);


    %%
    a=Imag(1:2:end);
    b=Imag(2:2:end);
    c=complex(a,b);
    %plot(a(512*400:512*401))

    % create 4 d matrix of csi grid
    mag_data=zeros(COL,ROW,THK,vecSize);

    k=2;
    for z=1:THK
        for y=1:ROW
            for x=1:COL
                if (k==2)
                    mag_data(x,y,z,:)=c(1:vecSize);
                else
                    mag_data(x,y,z,1:vecSize) = c((k-2)*vecSize+1:(k-1)*vecSize);
                end
                k=k+1;
            end
        end
    end
    f=zeros(COL,ROW,THK,vecSize);
    %%
    h=1;
    disp('Processing ....')
    rstart = ROW/2-ROW_P/2+1;
    rende = ROW/2+ROW_P/2;
    cstart = COL/2-COL_P/2+1;
    cende = COL/2+COL_P/2;
    tstart = THK/2-THK_P/2+1;
    tende = THK/2+THK_P/2;
    disp('Percentage');
    %%

    for z=tende:-1:tstart % read the z direction in the oposite way
        for y=rstart:rende
            for x=cstart:cende

                    f(x,y,z,:) = mag_data(x,y,z,:);

                mat=squeeze(f(x,y,z,:));
                %t=[x;y;z];
                filename = sprintf('/%02d_%02d_%02d.txt',x,y,z);
                fi=strcat(outdir,filename);
                fid = fopen(fi,'w+');
                zeroPh=-125.;
                fprintf(fid,'SamplingInterval: %1.2f\n',samInt);
                fprintf(fid,'TransmitterFrequency: %d\n',freq);
                %fprintf(fid, 'ZeroOrderPhase: \n', zeroPh );
                fprintf(fid,'MagneticField: 3T\n');
                fprintf(fid,'TypeOfNucleus: 1H\n');
                fprintf(fid,'NameOfPatient: %s\n',pat_name);

                fu=[real(mat)'; imag(mat)'];
                fprintf(fid, '%d %d\n', fu);
                fclose(fid);
                %for gg=1:vecSize
                %    fprintf(fid,' %d ', fu(gg,1));
                %    fprintf(fid,' %d\n', fu(gg,2));
                %end    ;
                %fwrite(fid,imag(mat(:)) real(mat(:)),'float');
                h=h+1;
            end
            percent=h/(COL*THK*ROW)*100;

            per=sprintf('%d',round(percent));
        end
        disp(per);
    end
end
