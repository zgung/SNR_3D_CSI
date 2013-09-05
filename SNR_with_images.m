% minarikova.lenka@gmail.com

function [] = SNR_with_images(directory, cho_ppm, bdwtd, trnct, control,field)

% clear all;
% directory = '~/Desktop/Current/Ivica2/01/';
% cho_ppm = 3.21;
% bdwtd = 50;
% trnct = 0;
% control = 0;
% field = 3;

% !!!!!!!!!!!!!!!!!! readme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% for working you need my other function called read_ascconv_lenk.m 
%       for reading parameters from the dicom
% directory = '~/Patient_name' - where is a directory called "Spec" with 
%       a dicom 3D CSI file and a directory called "Dixon_1.0iso_PAT2_v2_W"
%       for water DIXON images
% cho_ppm = 3.2 - exact position of Choline peak if the spectrum is set
%       to begin at 8.76 ppm and ends at 0.64, with 1000 Hz bandwidth
% bdwtd = 50 - bandwidth of the peak
% trnct = number of points to truncate at the end of FID
% control = 1 - at first you need to control if the baseline correction 
%    is working properly and everything is set all right, if you do not
%    wish to continue controling and you need to stop the script running,
%    try pressing "ctrl + c" that should stop it

% the output is maximal, mean value of all SNRs of Cho and a table 
% with all SNRs in one row, all saved in txt files in Spec directory

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% search variables in spectroscopy file:
tic;

directory_spec = strcat(directory,'Spec/');
cd(directory_spec);
disp(strcat('Processing:',directory));
%% <<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>
snr.cho = cho_ppm;
press_big = 0; % choose if you want more voxel (even that no totally inside the press box)
snr.noise = 7.4; % Aus diesem Bereich wird das Rausch Siganl genommen +0.5 bis -0.5
% Anfangs und Endewert fuer die ppm Skala eintragen
anfang = 8.76; % Anfang der ppm Skala 8.76 normally
ende = 0.64; % Ende der ppm Skala 0.64 
% truncate, replace the last few points in fid with 0s:

% read the .txt header from spectroscopic file made of the text part
% beginnig by "### ASCCONV BEGIN ###" and ending by "### ASCCONV END ###"
slices = dir(pwd);
for k = length(slices):-1:1 % find .IMA file
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    elseif strcat(fname(end-2),fname(end-1),fname(end)) ~= strcat('IMA')
        slices(k) = [];
    end
end

disp('Searching for spec file...');
spect = dir(pwd);
for k = length(spect):-1:1
    fname = spect(k).name;
    if fname(1) == '.' %|| strcat(spect(k).name) == strcat('Header.txt') || spect(k).name == strcat('Signal.txt')
        spect(k) = [];
    end
    % check filetype
    strg = (spect(k).name); % name of the file
    points = strfind(strg,'.');
    last_point = max(points);
    filetype_end = (strg(last_point:end));
    end_IMA=strcmp(filetype_end,'.IMA');    
    if end_IMA == 0;
        spect(k) = [];
        continue
    end    
end
voxel = read_ascconv_lenk(slices(1).name); % read parrameter from spectroscopy dicom

%% determine CSI-in-press parameters:
vecSize = voxel.vecSize;
voxel.size_z = voxel.FoV_z / voxel.number_z;
voxel.size_y = voxel.FoV_y / voxel.number_y;
voxel.size_x = voxel.FoV_x / voxel.number_x;
% number of steps in each direction:
if press_big == 1
    voxel.step_x = voxel.number_x / 2 - ceil((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x); %include also voxels touching the pressbox fringes
    voxel.step_y = voxel.number_y / 2 - ceil((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y);
    voxel.step_z = voxel.number_z / 2 - ceil((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z);
elseif press_big == 0
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x) - 1; %include only voxels inside the PRESS box
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y) - 1;
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z) - 1;
else
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x); %include voxels with PRESSbox fringes inside
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y);
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z);
end

%% Slices v0.3
% by Matthias Riha, modified by minarikova.lenka@gmail.com

% Nur fuor Traversal Ebenen zu gebrauchen!!!!
% Das File bf.m muess im selben Order sein!!!!!!!

%% DEFINITIONS

ROW = voxel.number_x; % x
COL = voxel.number_y; % y
SLC = voxel.number_z; % z

bbox.col_start = COL / 2 - voxel.step_y + 1; %Sagital von
bbox.col_end = COL / 2 + voxel.step_y; %bis
bbox.row_start = ROW / 2 - voxel.step_x + 1; %Coronal von
bbox.row_end = ROW / 2 + voxel.step_x; %bis
bbox.slc_start = SLC / 2 - voxel.step_z + 1; %Diese beiden Werte muessen immer gleich sein
bbox.slc_end = SLC / 2 + voxel.step_z;                   

%% Einlesen der Daten
csi.file_in = strcat(spect(1,1).name); %Pfad immer an die Datei anpassen
% read CSI data
fid = fopen(csi.file_in,'r');
fseek(fid,-((vecSize*ROW*COL*SLC*2*4)),'eof');
csi.data = fread(fid,'float32');
fclose(fid);

csi.real = csi.data(1:2:end);
csi.imag = csi.data(2:2:end);
csi.complex = complex(csi.real,csi.imag);
%% create 4D matrix of the whole csi grid
csi.mat_complex = zeros(ROW,COL,SLC,vecSize);
csi.mat_real = zeros(ROW,COL,SLC,vecSize);

o = -1;
for z = 1:SLC
    for y = 1:COL
        for x = 1:ROW
            o = o + 1;
            csi.mat_cmplx(x,y,z,:) = csi.complex(o * vecSize + 1:(o + 1) * vecSize);	
            csi.mat_rl(x,y,z,:) = csi.real(o * vecSize + 1:(o + 1) * vecSize);
        end
    end
end
%%
i = 0;
SNR.mid = 0;
SNR.tab = 0;
vecSize = 2 * vecSize;
anfang = anfang * (-1);
c = 0;

% Define the ppm scale
div = (anfang) + ende;
int = div / vecSize;
e = anfang;
F = 0;
for i = 1:vecSize
    F(1,i) = e;
    e = e - int;
end
F = F * (-1);
Vecsize = (1:vecSize);
Vecsize = Vecsize.';
F_col = (F.');
tabulk = [Vecsize,F_col];

for xx = 1:vecSize % define the values for baseline correction
    if round(100 * real(tabulk(xx,2))) == round(snr.cho * 100)
        r_sd = real(tabulk(xx - round(bdwtd / 2),1)); % right minimum next to choline peak
        l_sd = real(tabulk(xx + round(bdwtd / 2),1)); % left minimum next to Cho peak
        
        break
    end
end
SNR.main = 0;
lala = 0;
%% PRESSbox Daten FFT und in Txt File schreiben
for z = bbox.slc_start:bbox.slc_end
    b = 0;
    c = c + 1;
    for y = bbox.col_start:bbox.col_end
        a = 0;
        b = b + 1;
        for x = bbox.row_start:bbox.row_end
            lala = lala + 1;
            i = i + 1;
            a = a + 1;
            % 1 Voxel aus der 4D Matrix auslesen
            csi.rshpd_cmplx = reshape(csi.mat_cmplx(x,y,z,:),[],1);
            % FFT eines Voxels
            %x1 = csi.rshpd_cmplx;
            x2 = csi.rshpd_cmplx(1:vecSize / 2 - trnct,1);
            rest1 = zeros(trnct,1);
            x1 = [x2;rest1];
            N = vecSize;
            X = fft(x1,N);
            X = real(fftshift(X));
            %Xnixfit = X;
            diff_ = var([X(r_sd - 28:r_sd);X(l_sd:l_sd + 28)]); % variance between two means of points at the edge of peak

            % Baseline FIT !!!!!!
            if control == 1
                disp(SNR.main);
                X = bf(X,[16,32,50,80,110,140,170,200,230,260,290,320,350,r_sd - 28,r_sd - 21,...
                    r_sd - 14,r_sd - 7,r_sd - 3,r_sd,l_sd,l_sd + 3,l_sd + 7,l_sd + 14,l_sd + 21,...
                    l_sd + 28,970],5,'spline','confirm');
                plot(F,X);
                set(gca,'XDir','reverse');
            else
                X = bf(X,[16,32,50,80,110,140,170,200,230,260,290,320,350,r_sd - 28,r_sd - 21,...
                    r_sd - 14,r_sd - 7,r_sd - 3,r_sd,l_sd,l_sd + 3,l_sd + 7,l_sd + 14,l_sd + 21,...
                    l_sd + 28,970],5,'spline');
            end
                 
            % ppm scale + values from X
            test_snr = [Vecsize,F_col,X];
            yy = 0;
            zz = 0;
            w = 0;
            snr.min = 0;
            snr.max = 0;
            snr.base = 0;
            for xx = 1:vecSize
                snr.schleife = test_snr(xx,:,:);
                ppm_wert = snr.schleife(:,2);
                if ppm_wert < snr.noise + 0.5
                    if ppm_wert > snr.noise - 0.5
                        yy = yy + 1;
                        snr.min(yy) = real(snr.schleife(:,3));
                    end
                end
                if ppm_wert < snr.cho + 0.15
                    if ppm_wert >snr.cho - 0.15
                        zz = zz + 1;
                        snr.max(zz) = real(snr.schleife(:,3));
                    end
                end
            end
            %SNR berechnung
            [~,snr.cho_max_ppm] = max(snr.max(2:length(snr.max) - 1));
            snr.cho_max_val = snr.max(snr.cho_max_ppm:snr.cho_max_ppm + 2);
            ima = mean(snr.cho_max_val);
            ims = std(snr.min(:));
            SNR.main = (ima) / (ims * 2);

            
            if diff_ > 10000
                SNR.main = 1;
            end
            if SNR.main < 2
                SNR.main = 1;
            end
            SNR.main;
            SNR.mid = SNR.mid + SNR.main;
            SNR.tab(b,a,c) = SNR.main;
            SNR.reshaped(lala) = SNR.main;
        end
    end
end
%% treshold for all values:
%reshape and add coordinates as are in siemens:
c.pix_width = voxel.step_y * 2; % number of voxels in the pressbox - x axis
c.pix_height = voxel.step_x * 2; % -||- - y axis
c.pix_depth = voxel.step_z * 2; % -||- - z axis
ooo = 0;
for k = 1:c.pix_depth
    for j = 1:c.pix_width
        for i = 1:c.pix_height
            ooo = ooo + 1;
            SNR.w_coor(ooo,1) = i + (voxel.number_x / 2 - voxel.step_x);
            SNR.w_coor(ooo,2) = j + (voxel.number_y / 2 - voxel.step_y);
            SNR.w_coor(ooo,3) = k + (voxel.number_z / 2 - voxel.step_z);
            SNR.w_coor(ooo,4) = SNR.reshaped(1,ooo);
        end
    end
end
%% saving important things:
% path = sprintf('Output_choSNR.txt');
% fid_write = fopen(path,'w');
% fprintf(fid_write,'%d\n',SNR.reshaped);
% fclose(fid_write);

pvc_mean = mean(SNR.reshaped);
path = sprintf('Output_cho_mean.txt');
fid_write = fopen(path,'w');
fprintf(fid_write,'%d\n',pvc_mean);
fclose(fid_write);
disp(pvc_mean);

pvc_max = max(SNR.reshaped);
path = sprintf('Output_cho_max.txt');
fid_write = fopen(path,'w');
fprintf(fid_write,'%d\n',pvc_max);
fclose(fid_write);
disp(pvc_max);
% save the bandwith and Cho ppm for the future
path = sprintf('Cho_ppm.txt');
fid_write = fopen(path,'w');
fprintf(fid_write,'%d\n',cho_ppm);
fclose(fid_write);
path = sprintf('Cho_bdwtd.txt');
fid_write = fopen(path,'w');
fprintf(fid_write,'%d\n',bdwtd);
fclose(fid_write);
dlmwrite('Output_choSNR_w_coor.txt', SNR.w_coor, 'delimiter', '\t', ...
         'precision', 6);

%% The second part:
% check coil elements:
if field == 3 % anyway we are not measuring dixons on 7 T
    if voxel.coilel3 ~= 0
        disp('!!Both coils ON!!');
    else
        if strcat(voxel.coilel1(4)) == strcat('R') && ...
                strcat(voxel.coilel2(4)) == strcat('R')
            disp('Right coil elements: ON');
        end
        if strcat(voxel.coilel1(4)) == strcat('L') && ...
                strcat(voxel.coilel2(4)) == strcat('L')
            disp('Left coil elements: ON');
        end
    end
else
    disp('I am not programmed to check coil elements on another field strength')
end
%% move to dixom folder: water first
disp('Processing images...');
cd(strcat(directory,'Dixon_1.0iso_PAT2_v2_W/'));
slices = dir(pwd);
for k = length(slices):-1:1
    fname = slices(k).name;
    if fname(1) == '.'
        slices(k) = [];
    end
end
% check if the slices are in a good order
Afields = fieldnames(slices);
Acell = struct2cell(slices);
sz = size(Acell); % Convert to a matrix
Acell = reshape(Acell, sz(1), []);      % Px(MxN)
% Make each field a column
Acell = Acell';                       % (MxN)xP
for k = 1:length(slices)
    [~,~,~,~,slc_n,~,~,~,~,~,~,~,~,~] = strread(slices(k,1).name,'%s %s %s %d %d %d %d %d %d %d %d %d %d %s','delimiter','.');
    Acell{k,4}= slc_n;
end
Acell = sortrows(Acell, 4); % Sort by first field "name"
Acell = reshape(Acell', sz); % Put back into original cell array format
slices = cell2struct(Acell, Afields, 1); % Convert to Struct
clear Acell; clear Afields;
%% determine CSI-in-press parameters:
voxel.size_z = voxel.FoV_z / voxel.number_z; % the size of 1 voxel after zero filling in mm
voxel.size_y = voxel.FoV_y / voxel.number_y;
voxel.size_x = voxel.FoV_x / voxel.number_x;
% number of steps in pressbox after zero filling in each direction / 2:
    voxel.step_x = voxel.number_x / 2 - fix((voxel.FoV_x - voxel.p_fov_x) / 2 / voxel.size_x) - 1; %include only voxels inside the PRESS box
    voxel.step_y = voxel.number_y / 2 - fix((voxel.FoV_y - voxel.p_fov_y) / 2 / voxel.size_y) - 1;
    voxel.step_z = voxel.number_z / 2 - fix((voxel.FoV_z - voxel.p_fov_z) / 2 / voxel.size_z) - 1;

nfo1 = dicominfo(slices(1).name);
% %%%%%%%%%%%%%%%%%%%%%%% Import the W images %%%%%%%%%%%%%%%%%%%%%%%%%
n = numel(slices);
imgs = cell(1,1);
for ii = 1:n
    imgs{1,1}(:,:,ii) = double(dicomread(slices(ii).name));
end
%% interpolate so 1 mm ~ 1 image voxel
[XI,YI,ZI] = meshgrid(1:(floor(10000 * 1 / nfo1.PixelSpacing(1,1)) / 10000):length(imgs{1,1}(1,:,1)),...
    1:(floor(10000 * 1 / nfo1.PixelSpacing(2,1)) / 10000):length(imgs{1,1}(:,1,1)),...
    1:(floor(10000 * 1 / nfo1.SliceThickness) / 10000):length(imgs{1,1}(1,1,:)));
imgs{2,1} = ba_interp3(imgs{1,1},XI,YI,ZI,'nearest'); % fast interpolation to 1x1x1 mm3 
[YI,XI,~] = size(imgs{2,1});
% make the center of coordinates the centre of each axial slice: horizontaly
if -nfo1.ImagePositionPatient(1,1) > XI -  ...
        -nfo1.ImagePositionPatient(1,1) == 1
        % add zeros to the end
        imgs{2,1}(:,XI + 1:-round(2 * ...
            nfo1.ImagePositionPatient(1,1)),:) = 0;
else
        % move image right and add zeros to the begging
        imgs{2,1}(:,XI - -round(2 * nfo1.ImagePositionPatient(1,1)) + ...
            1:2 * XI - -round(2 * nfo1.ImagePositionPatient(1,1)),:) = ...
            imgs{2,1}(:,1:XI,:);
        imgs{2,1}(:,1:XI - -round(2 * nfo1.ImagePositionPatient(1,1)),:) = 0;
end
% now, center image in vertical dimension
if -nfo1.ImagePositionPatient(2,1) > YI - -nfo1.ImagePositionPatient(2,1) == 1
        %add zeros to the end
        imgs{2,1}(YI + 1:-round(2 * ...
            nfo1.ImagePositionPatient(2,1)),:,:) = 0;
else
        % move image down and add zeros to the begging
        imgs{2,1}(YI - -round(2 * nfo1.ImagePositionPatient(2,1)) ...
            + 1:2 * YI - -round(2 * nfo1.ImagePositionPatient(2,1)),:,:) = ...
            imgs{2,1}(1:YI,:,:);
        imgs{2,1}(1:YI - -round(2 * nfo1.ImagePositionPatient(2,1)),:,:) = 0;
end
%% you need to cut the images in voxels dimensions (not in press box!)
c.cut_width = round(2 * voxel.step_y * voxel.size_y); % width in pixels of the cuted image
c.cut_height = round(2 * voxel.step_x * voxel.size_x); % height -||-
c.cut_depth = round(2 * voxel.step_z * voxel.size_z); % depth -||-
c.pix_width = voxel.step_y * 2; % number of voxels in the pressbox - x axis
c.pix_height = voxel.step_x * 2; % -||- - y axis
c.pix_depth = voxel.step_z * 2; % -||- - z axis
voxel.FoV_z_1 = round(voxel.fov_cntr_z + c.cut_depth / 2); % the real position in number of pixels in the image
voxel.FoV_z_2 = round(voxel.fov_cntr_z - c.cut_depth / 2);
voxel.csi_slice_1 = round((nfo1.ImagePositionPatient(3,1) - voxel.FoV_z_1) + 1) + 0; %1st - fhead
voxel.csi_slice_last = round((nfo1.ImagePositionPatient(3,1) - voxel.FoV_z_2)) + 0; %2nd - from feet

% check if the values are negative and turn them positive
if voxel.csi_slice_1 < 0
    ttt = voxel.csi_slice_1;
    voxel.csi_slice_1 = abs(voxel.csi_slice_last);
    voxel.csi_slice_last = abs(ttt);
    clear ttt;
end
% cut the voxels-images only in Z direction:
imgs{3,1} = imgs{2,1}(:,:,voxel.csi_slice_1:voxel.csi_slice_last); % remove the additional slices
%%


% look if the csi is rotated about an angle:
voxel.angle_deg = radtodeg(voxel.angle);
% the center point of csi in press box in rotated images has coordinates from the
% left upper corner:
if voxel.angle > 0.0 % clockwise
    voxel.fov_cntr_x = voxel.fov_cntr_x + 4;
    voxel.fov_cntr_y = voxel.fov_cntr_y + 3; % plus posunie ratio hore
elseif voxel.angle < -0.0 % anticlockwise
    voxel.fov_cntr_x = voxel.fov_cntr_x + 3;
    voxel.fov_cntr_y = voxel.fov_cntr_y + 4; % plus posunie ratio hore
else
    voxel.fov_cntr_x = voxel.fov_cntr_x + 3;
    voxel.fov_cntr_y = voxel.fov_cntr_y + 3; % plus posunie ratio hore
end
for ii = 1:numel(imgs{3,1}(1,1,:))
% rotate default pictures:
   imgs{3,1}(:,:,ii) = imrotate(imgs{3,1}(:,:,ii),-voxel.angle_deg,...
       'nearest','crop');
end
imgs_w = imgs{3,1};
voxel.fov_cntr_x_rttd = numel(imgs{3,1}(1,:,1)) / 2 - (-voxel.fov_cntr_x) * cos(-voxel.angle) - ...
     (-voxel.fov_cntr_y) * sin(-voxel.angle) + 0;
voxel.fov_cntr_y_rttd = numel(imgs{3,1}(:,1,1)) / 2 - (-voxel.fov_cntr_y) * cos(-voxel.angle) + ...
     (-voxel.fov_cntr_x) * sin(-voxel.angle) + 0;
clear imgs;
clear slices;
shft_lr = 0;
shft_ud = 0;
voxel.fov_x1 = floor(voxel.fov_cntr_x_rttd - voxel.step_x * voxel.size_x + round(shft_lr * (voxel.FoV_x / voxel.notinterpfov_x))); % voxel values are beeing saved for psf_pics.m!
voxel.fov_y1 = floor(voxel.fov_cntr_y_rttd - voxel.step_y * voxel.size_y + round(shft_ud * (voxel.FoV_y / voxel.notinterpfov_y)));


%% teraz by tu mala ist segmentacia, ale namiesto toho budes len citat
% obrazky na pozadie...

clearvars -except directory cho_ppm bdwtd trnct control field nfo1 directory_spec voxel c SNR;
%%
cd(directory_spec);
% import from simple text file, only amplitude data from jmrui text result
% file named Output_choSNR.txt
%fid = fopen('Output_choSNR.txt','r');
load(strcat(nfo1.PatientName.FamilyName,'_voxel.mat'));
load(strcat(nfo1.PatientName.FamilyName,'_density.mat'));
load(strcat(nfo1.PatientName.FamilyName,'_imgs_w.mat'));
load(strcat(nfo1.PatientName.FamilyName,'_density_wt_cor'));
load(strcat(nfo1.PatientName.FamilyName,'_rtio'));
mtrx1 = SNR.reshaped; %textscan(fid,'%f'); % this is right? not just '%f' ??
%fclose(fid);
aa = 0; % Z direction in siemens MRSI is opposite to image Z direction
for k = voxel.step_z * 2:-1:1
   for j = 1:voxel.step_y * 2
        for i = 1:voxel.step_x * 2
            aa = aa + 1;
            density(j,i,k) = mtrx1(aa);
        end
    end
end

% upsample choline signal information without interpolating
fctr = 2;
f_nas = fctr * voxel.size_x;
B = density;
A = cell(f_nas,1);
for i = 0:f_nas - 1
    A{i + 1} = upsample(B,f_nas,i);
end
B = sum(cat(f_nas,A{:}),f_nas);
B = permute(B,[3 2 1]);
for i = 0:f_nas - 1
    A{i + 1} = upsample(B,f_nas,i);
end
B = sum(cat(f_nas,A{:}),f_nas);
B = permute(B,[2 1 3]);
for i = 0:f_nas - 1
    A{i + 1} = upsample(B,f_nas,i);
end
B = sum(cat(f_nas,A{:}),f_nas);
B = permute(B,[3 1 2]);
signal3D_ = B;
clear B;
%% interpolate images so they have 4 times bigger resolution
c.img_s = size(imgs_w);
[XI,YI,ZI] = meshgrid(1:(c.img_s(2) - 0.9) / (c.img_s(2) * fctr):c.img_s(2),...
    1:(c.img_s(1) - 0.9) / (c.img_s(1) * fctr):c.img_s(1),...
    1:(c.img_s(3) - 0.9) / (c.img_s(3) * fctr):c.img_s(3));
imgs_w = ba_interp3(imgs_w,XI,YI,ZI,'linear');
clear XI; clear YI; clear ZI;
%% define the colormaps for background MRI picture and metabolit map
% than find a maximal value of background and signal images:
h = max(max(max(imgs_w)));
g_test = max(max(max(signal3D_)));
if g_test < 10
    g = max(max(max(signal3D_))) * 100; % amplified signal values
else
    g = g_test;
end
cmp = colormap(gray);
%[X,Y] = meshgrid(1:3,1:64);
[XI,YI] = meshgrid(1:3,1:64 / h:64);

%%
cmp1 = interp2(cmp,XI,YI);
jet1 = colormap(jet);
%%
%[X,Y] = meshgrid(1:3,1:64);
[XI,YI] = meshgrid(1:3,1:64 / g:64);
jet2 = interp2(jet1,XI,YI);

%% draw 8 slices at once:
figure('units','normalized','position',[.05 .1 .9 .9])
for aa = 2:c.pix_depth + 1
    %subplot('Position',[(-0.0 ) 0.05 (2.5 / c.pix_depth) (2.5 / c.pix_depth)]);
    subaxis(round(c.pix_depth / 3 - 1),5 ,aa - 1, 'Spacing', 0, 'Padding', 0, 'Margin', 0);
    subimage(imgs_w(:,:,fix((aa - 1) * (c.cut_depth * fctr) / (c.pix_depth + 1)) + 5),cmp1);
    pic = signal3D_(:,:,fix((aa - 1) * (c.cut_depth * fctr) / (c.pix_depth + 1)) - 0);
    hold on
    if g_test < 10
        hh2 = image((voxel.fov_x1 - 3) * fctr,(voxel.fov_y1 - 3) * fctr,pic .* 100,...
            'AlphaData',0.5);
    else
        hh2 = image((voxel.fov_x1 - 3) * fctr,(voxel.fov_y1 - 3) * fctr,pic,...
            'AlphaData',0.5);
    end

%        'AlphaData',gradient(pic),...
%        'AlphaDataMapping','scaled',...
%        'AlphaData',(abs(pic) * 100),...

    colormap(jet2);
    axis([100 800 100 800]);
    %axis([0 length(imgs_w{1,1}) 0 numel(imgs_w{1,1}(:,1))]);
    caxis([(0 * g) (1 * g)]);
    hold off
    axis off
end






toc