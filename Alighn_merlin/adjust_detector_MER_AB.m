function   adjust_detector_MER_AB(tube_fit,detector)
%alters the theta and phis in the detector.dat file to calibrate for the
%non perfect electronics. Tube_fit is a file containing the fitted results
%of the calibration run on the tubes. This file is produced by running
%tube_calibrate.m
% detector is the non corected detector.dat file as produced by Tobies
% program. THe output file det_corrected.dat is the new corrected detector.dat
% as corrected for the imperfect detector electronics.

tub_fit = load(tube_fit,'r');
fid = fopen(detector,'r');
clobDet = onCleanup(@()fclose(fid));
line1 = fgets(fid)
line2 = fgets(fid)
line3 = fgets(fid)
det = textscan(fid,'%d %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

tube=fix(double(det{1,1}/10000.0));

num_tubes=length(tub_fit(:,1)); %number of tubes
deg_vs_rad=180/pi;% 57.29578;


tube_len_short=1.234; %active length of short packs
tube_len_long =2.900; % active length of a MERLIN detector (m)
low_short_tube_pos   = -1.5; % the position of
high_short_tube_pos  = 0.17;
%z_shift=0.833;   %the vertical shift of the short packs of the centre line (in m)
nDet = 512;
xbin_short=tube_len_short/nDet;% size of each bin for short packs on door 3
xbin=tube_len_long/nDet;% size of each bin in mm on a MERLIN 3m detector

% relative positions of the lines from calibration bin lines on MERLIN detectors
% currently approximated; correct values needed. Boundary lines fitted to
% line-edges positions and internal - to the line centers.
pos_correct = [-0.4800,-0.2100,-0.1400, 0, 0.1500, 0.2800, 0.4800];
% physical positions (in meters) where the ideal lines should be
pos_correct = pos_correct*tube_len_long;
%

for count=1:num_tubes
    tube_fix=tub_fit(count,1); % this specifies the tube to alter
    
    fit_par   = tub_fit(count,3:end);
    tube_type = tub_fit(count,2); % tube type, 0 -- long tube, 1 short upper tube, 2 short lower tube
    
    valid = ~(isnan(fit_par) | fit_par==0);
    if (sum(valid)>=2)  % means fit program could not fit the tube
        %
        index=find(tube==tube_fix);  % indexes all the detectors for that tube
        dist_old=(det{1,3} (index,1));
        theta_old=(det{1,5} (index,1));
        thi_old=(det{1,6} (index,1));
        %
        [x_old,z_old,y_old]=sph2cart((thi_old./deg_vs_rad),((theta_old-90)./deg_vs_rad),dist_old);
        %
        switch(tube_type)
            case 0
                fit_pos = ((fit_par(valid)-256)/512)*tube_len_long ; % position of fitted detectors as function of tube length
                ideal_det_pos = -0.5*tube_len_long + 0.5*xbin+ (0:(nDet-1))*xbin;
            case 1
                fit_pos = high_short_tube_pos+(fit_par(valid)/512)*tube_len_short;
                ideal_det_pos = high_short_tube_pos + 0.5*xbin_short+ (0:(nDet-1))*xbin_short;
            case 2
                fit_pos = low_short_tube_pos+(fit_par(valid)/512)*tube_len_short;
                ideal_det_pos = low_short_tube_pos + 0.5*xbin_short+ (0:(nDet-1))*xbin_short;
            otherwise
                error('ADJUST_DETECTORS:invalid_argument',...
                    ' unknown tube type %d obtained',tube_type);
        end
        z_new = tube_correction(pos_correct(valid),fit_pos,ideal_det_pos);  %new positions along tube
        
        
        [thi_new,theta_new,dist_new] = cart2sph(x_old,z_new',y_old);
        thi_new=thi_new*deg_vs_rad;
        theta_new=(theta_new*deg_vs_rad)+90;
        
        
        
        det{1,3} (index,1)=dist_new;
        det{1,5} (index,1)=theta_new;
        det{1,6} (index,1)=thi_new;
    end
end

fid2 = fopen('det_corrected_2017_3.dat','w');

fprintf(fid2,'%s',line1);
fprintf(fid2,'%s',line2);
fprintf(fid2,'%s',line3);
aa=[];
for i=1:(length(det))
    aa(:,i)=det{1,i};
end

fprintf(fid2,'%7d %8.3f %10.5f %5d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f \n',aa');

fclose(fid2);



function  z_new = tube_correction(pos_correct,bin_measure,ideal_det_pos)
% ideal_det_pos -- the positions where ideal detectors should be
% bin_measure   -- real coordinates of the lines
% pos_correct   -- ideal positions, where lines should be

% figure
%figure(1);errorbar(bin_measure,pos_correct,ones(size(bin_measure)),'o');     % Plot dat
if numel(bin_measure) >3
    p=polyfit(bin_measure,pos_correct,3);
else
    p=polyfit(bin_measure,pos_correct,2);
end
z_new=polyval(p,ideal_det_pos);
%figure(1);hold on;plot(xv,z_new,'r'); % Plot fitted results


