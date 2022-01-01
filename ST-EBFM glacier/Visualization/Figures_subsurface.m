%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURES: SUBSURFACE
%%% Plots the time-evolution of subsurface conditions at a specified
%%% location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER SPECIFICATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdir = '..\Output\';              % output directory
%outdir = 'E:\Output_1999_with_aditions' ; 
coords = [396496 7534700];%[398156 7534170];          % location (UTM)
tsub_start = '15-Sep-1997 0:00';     % start date
tsub_end = '15-Sep-1999 0:00';       % end date
max_depth = 0.2 ;                     % maximum subsurface depth in plot (at 
                                    %   starting time)
vars = {'subT';'subD';'subW'};      % variables to plot (choose from 'subT'
                                    %   , 'subD', 'subW' and/or 'subS')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%folderB = 'D:\Output_1999_with_aditions' ; 
load([outdir '\runinfo.mat']);
%load([folderB '\runinfo.mat']); 


L = length(grid.mask);
time_start = time.ts;
time_end = time.te;
period = IOout.freqout*time.dt;
tvec = datenum(time_start):period:datenum(time_end);
T = (datenum(time_end) - datenum(time_start))/period;

diff_y = grid.y(:)-coords(2);
diff_x = grid.x(:)-coords(1);
diff = sqrt(diff_y.^2 + diff_x.^2);
gp = find(diff == min(diff));

nl = grid.nl;

taxis_model = datetime(datevec(datenum(tsub_start):...
    IOout.freqout*time.dt:datenum(tsub_end)+0.5*IOout.freqout*time.dt));
tind = [];
taxis_plot_datenum = [];
for tt=1:length(taxis_model)
    if min(abs(datenum(IOout.output_times)-datenum(taxis_model(tt))))<...
            IOout.freqout*time.dt
        [~,tind(end+1)] = min(abs(datenum(IOout.output_times)-datenum(...
            taxis_model(tt))));
        taxis_plot_datenum(end+1) = datenum(IOout.output_times(tind(end)));
    end
end
taxis_plot = datetime(datevec(taxis_plot_datenum));

AsubZ = nan(length(tind),nl);
fid = fopen([outdir '/OUT_subZ.bin'],'rb');
for t=1:length(tind)
    fseek(fid,(tind(t)-1)*4*L*nl,'bof');
    Atemp_ = fread(fid,[nl,L],'real*4','l');
    AsubZ(t,:) = Atemp_(:,gp);
end
fclose(fid);
AsurfH = zeros(length(tind),1);
fid = fopen([outdir '/OUT_surfH.bin'],'rb');
for t=1:length(tind)
    fseek(fid,(tind(t)-1)*4*L+(gp-1)*4,'bof');
    AsurfH(t,1) = fread(fid,1,'real*4','l');
end
fclose(fid);

Depth_up = -cumsum(AsubZ,2)+AsubZ+repmat(AsurfH,1,grid.nl);
Depth_down = -cumsum(AsubZ,2)+repmat(AsurfH,1,grid.nl);
Depth_up = Depth_up - min(Depth_down(:));
Depth_down = Depth_down - min(Depth_down(:));
z_int = min(Depth_down(:)):min(AsubZ(:,2:end),[],'all')/5:max(Depth_up(:));
z_int = repmat(z_int,length(tind),1);

figure;
for v=1:length(vars)
    Asub = nan(length(tind),nl);
    fid = fopen([outdir '/OUT_' vars{v} '.bin'],'rb');
    for t=1:length(tind)
        fseek(fid,(tind(t)-1)*4*L*nl,'bof');
        Atemp_ = fread(fid,[nl,L],'real*4','l');
        Asub(t,:) = Atemp_(:,gp);
    end
    fclose(fid);
    
    if any(strcmp(vars{v},{'subW','subS'}))
        Asub = Asub ./ AsubZ;
    end
    
    Asub_regrid = nan(size(z_int));
    for t=1:length(tind)
        for z=1:grid.nl
            ind = find(z_int(t,:)>=Depth_down(t,z) & z_int(t,:)<=Depth_up(t,z));
            Asub_regrid(t,ind) = Asub(t,z);
        end
    end
    zaxis = min(Depth_down(:)):min(AsubZ(:,2:end),[],'all')/5:max(Depth_up(:));
    taxis = datenum(taxis_plot);
    
    s = subplot(length(vars),1,v);
    pcolor(taxis,zaxis-max(zaxis)+max(AsurfH),Asub_regrid'); hold on;
    shading flat;
    datetickzoom('x','dd-mmm-yyyy');
    xlabel('Date');
    ylabel('Height (m)');
    c = colorbar;
    if strcmp(vars{v},'subD')
        colormap(s,cbrewer('seq','YlGnBu',128,'spline'));
        ylabel(c,'Density (kg m^{-3})');
    elseif strcmp(vars{v},'subT')
        colormap(s,flipud(abs(cbrewer('div','RdBu',128,'spline'))));
        ylabel(c,'Temperature (K)');
    elseif strcmp(vars{v},'subW')
       colormap(s,flipud(abs(cbrewer('seq','YlGnBu',128,'spline'))));       % used to be cbrewer('div','RdBu',128)
        ylabel(c,'Irreducible water content (kg m^{-3})');
    elseif strcmp(vars{v},'subS')
        colormap(s,flipud(abs(cbrewer('seq','YlGnBu',128,'spline'))));      % same. call cbrewer() to see available colormaps
        ylabel(c,'Slush water content (kg m^{-3})');
    end
    xlim([datenum(tsub_start) datenum(tsub_end)]);
    ylim([-max_depth max(zaxis-max(zaxis)+max(AsurfH))]);
    plot(taxis,AsurfH,'-k','LineWidth',0.1);
    set(gca,'Layer','top');
end
%%
figure
load([IOout.homedir '\Grid\stor.mat']); 
contour(stor.UTMeast,stor.UTMnorth,stor.elev,[1000:20:2000],'k'),hold on
scatter(coords(1),coords(2),200,'xr','linewidth',1.5)
axis equal, box on 
limeast = stor.UTMeast(~isnan(stor.UTMeast));
limnorth = stor.UTMnorth(~isnan(stor.UTMnorth));
xlim([min(grid.x) max(grid.x)]);
ylim([min(grid.y) max(grid.y)]);
