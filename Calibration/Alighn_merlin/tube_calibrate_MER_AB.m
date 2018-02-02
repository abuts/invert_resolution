function fc = tube_calibrate_MER_AB(run,tmin,tmax,varargin)
% takes a calibration run and finds the centre and end points of the
% selected tubes. varargin can take the form of 'bank' 'pack' 'tube'. If
% no arguments then whole detector array is fitted. If bank is used only
% tubes in that bank are fitted and so on. The data are fitted over the
% integrated time range of tmin to tmax.
% An output file is created called results.dat which gives 4 columns. The
% first is the tube index number followed by the fitted start middle and
% end points.


ass run
len=length(varargin);
% target file
fname = 'cal_results';

spec_num=gget('spec');  % spectra number
det_num=gget('udet'); %detecot number
det_tubes=fix(double(det_num)/10000);
packstart=[1 1 1 1 1 1 1 1 1]; % pack to start fitting
packend=[4 4 5 4 4 4 4 3 3]; % defines number of packs per door
tubestart=1;
tubeend=8;
% define the peak structure and functions, which describe peaks, for each
% kind of tube.
peak_pos_map= containers.Map();
peak_pos_map('normal_tube') =  {...
    {@step_fun,@gaus_and_bg,@gaus_and_bg,@gaus_and_bg,@gaus_and_bg,@gaus_and_bg,@step_fun},...
    {@guess_for_step,@guess_for_gauss,@guess_for_gauss,@guess_for_gauss,@guess_for_gauss,@guess_for_gauss,@guess_for_step},...
    [0.03,0.29,0.36,0.5,0.65,0.78,0.98]};
peak_pos_map('upper_short')  = {...
    {@gaus_and_bg,@gaus_and_bg,@step_fun},... % the functions, which describe peaks.
    {@guess_for_gauss,@guess_for_gauss,@guess_for_step},... % the functions, which calculate the guess peak parameters.
    [ 0.2, 0.50, 0.96]};   % approximate relative peak positions along the tube
peak_pos_map('lower_short')  = {...
    {@step_fun,@gaus_and_bg,@gaus_and_bg},...
    {@guess_for_step,@guess_for_gauss,@guess_for_gauss},...
    [ 0.04, 0.63, 0.79]};

if len ==0
    doorstart=1;
    doorend=9;
end
if len > 0
    bank=varargin{1};
    doorstart=bank;
    doorend=bank;
    fname = [fname,'_bank',num2str(doorstart)];
end
if len >1
    pack=varargin{2};
    packstart(bank)=pack;
    packend(bank)=pack;
    fname = [fname,'_pack',num2str(pack)];    
end
if len >2
    tube=varargin{3};
    tubestart=tube;
    tubeend=tube;
    fname = [fname,'_tube',num2str(tube)];        
end
fid = fopen([fname,'.dat'],'w');
clob = onCleanup(@()fclose(fid));

fc = fig_spread('-tight');
n_plot_tubes = 0;
fig_shift = 0;
tube_id_old = -1;
fig_range_c = 0;
for banks=doorstart:doorend
    for packs=packstart(banks):packend(banks)
        for tubes=tubestart:tubeend
            tubeget=str2double([num2str(banks);num2str(packs);num2str(tubes)]');
            fprintf(' tube  num = %d\n',tubeget);
            index=find(det_tubes==tubeget);
            spec_tube=spec_num(index);
            spec_min=double(min(spec_tube));  %min spec for selected tube
            spec_max=double(max(spec_tube));   % max spec for selected tube
            
            spectube=get(sumspec(spec_min,spec_max,tmin,tmax)); %sumspec orders the spectra in accending order
            yval=spectube.y;
            if (max(yval)~= 0)
                % get tube type given tube position in the pack.
                [tube_type,tube_id,fig_range] = get_type(banks,packs,tubes);
                par = peak_pos_map(tube_type);
                %------------------------------------------------------
                % place figures nicely
                if tube_id ~= tube_id_old
                    fig_shift = fig_shift+fig_range_c;
                    tube_id_old = tube_id;
                end
                fig_range_c = fig_range;
                    
                n_plot_tubes  = n_plot_tubes +1;
                if n_plot_tubes > 20 % plot no more than 20 tubes per picture
                    fig_shift  = fig_shift+fig_range;
                    n_plot_tubes  = 0;
                end
                %-------------------------------------------------------
                [peak_pos,fc] =fit_data(par{1},par{2},par{3},yval,20,fig_shift,fc);
                
                if numel(peak_pos) == 3
                    fprintf(' left end= %f; middle = %f; right_end= %f\n',...
                        peak_pos);
                    all_pos = ones(7,1)*NaN;
                    if tube_id == 1
                        all_pos(end-2:end) = peak_pos;
                    elseif tube_id == 2
                        all_pos(1:3) = peak_pos;
                    end
                else
                    fprintf(' peak_pos = %f; %f; %f; %f; %f; %f; %f\n',...
                        peak_pos);
                    all_pos = peak_pos;
                end
                fprintf(fid,'%3.0f %2d %f %f %f %f %f %f %f\n',tubeget,tube_id,all_pos);
            end
            
            
        end
    end
end
end

function [tube_type,tube_id,fig_range] = get_type(banks,pack,tubes)
% given the tube position in the door and the door number get the tube
% type.
%
% to simplify exchange with adjust detectorm, define tube_id as 0 for normal
% tube, 1 for upper and 2 for lower short
%
tube_type = 'normal_tube';
tube_id   = 0;
fig_range = 7;
if banks == 3
    if pack == 1
        tube_type = 'upper_short';
        tube_id  = 1;
        fig_range  = 3;
    elseif pack == 2
        tube_type = 'lower_short';
        fig_range  = 3;
        tube_id  = 2;
    end
end
end

function chisc = err_calc(func,par,x,I,dI)
% calculate chi-squared distance between calculated and measured functions.
yth = func(par,x);
% Calculate chi squared
chisc =sum(((I-yth).^2)./dI.^2);
end
%
function val = gaus_and_bg(p,x)
val =p(4)+p(1)*exp(-(x-p(2)).^2/(2*p(3)*p(3)));
end
%
function par = guess_for_gauss(x,s)
% calculate approximate parameters of a gaussian, describing input data.
par = zeros(4,1);

par(4)=min(s);
s = s - par(4);
av = sum(s);
par(1)=max(s) ;
par(2)=sum(s.*x)/av;   % Used to be peak
%xp = x -par(2);
%par(3)=sqrt(sum(s.*xp.*xp)/av);
par(3) = 3;

%par(5)=(s(1)-s(length(x)))./x(1)-x(length(x));   % Added here!!!

end

%
function val = step_fun(p,x)
if (p(4) <0)
    p(4)=0;  %this variable cannot go below zero
end

if (p(1) > 0)
    val=(p(1)*erfc((p(2)-x)/p(3)))+p(4);
end
if (p(1) < 0)
    val=((-2*p(1))+(p(1)*erfc((p(2)-x)/p(3))))+p(4);
end
end

function par = guess_for_step(x,s)
% calculate approximate parameters of a step function, describing input data.
bg = min(s);
s = s-bg;
par = zeros(4,1);
der = s(2:end)-s(1:end-1);
av  = sum(der);
xp = 0.5*(x(2:end)+x(1:end-1));
par(2) = sum(der.*xp)/av;
ampl = max(s);

par(1) = 0.5*(ampl)*sign(av);
xp = xp - par(2);
sigmaSq = sum(der.*xp.*xp)/av;
if sigmaSq > 0
    par(3) = sqrt(sigmaSq );
else
    par(3) = 6;
end
par(4) = bg;

end




function [peak_pos,fig_pos_controller] =fit_data(func_array,guess_calc_array,pos0,Intensity,search_range,...
    fig_range,fig_pos_controller)

error=sqrt(Intensity);
error(error == 0) = 1;
lx = numel(Intensity);
position=(1:lx)';



peak_pos = zeros(numel(func_array),1);
for i=1:numel(func_array)
    
    range0 = floor(pos0(i)*lx);
    x0 = range0-search_range;
    if x0 <2; x0=2;end;
    x1 = range0+search_range;
    if x1 >lx-1; x1=lx-1;end;
    maxint=max(Intensity(x0:x1));
    if maxint == 0  % if there are no counts in  tube
        continue;
    end
    minim = x0;
    maxim = x1;
    I=Intensity(minim:maxim);% Intensity array
    dI=error(minim:maxim);%  Error in intensity array
    pos = position(minim:maxim);
    
    n_fig = fig_range+i;
    figure(n_fig);errorbar(pos,I,dI,'o');     % Plot data
    guess_calc = guess_calc_array{i};
    % Start Parameters
    p = guess_calc(pos,I);
    func = func_array{i};
    
    options=optimset('MaxIter',1000,'TolFun',...
        1*10^(-6),'MaxFunEvals',5000);           % Set options
    [the_fit,fval,flag]=fminsearch(@(par,x,I,dI)err_calc(func,par,x,I,dI),p,options,pos,I,dI) ;      % Minimize chi_squared
    if (flag ~=1)
        continue
    end
    peak_pos(i) = the_fit(2);
    % Get theory and plot
    x=min(pos):0.5:max(pos);
    Ith=func(the_fit,x);
    h=figure(n_fig);hold on;plot(x,Ith,'r'); % Plot fitted results
    fig_pos_controller=fig_pos_controller.place_fig(h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
