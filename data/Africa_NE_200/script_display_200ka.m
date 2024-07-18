% Display data
% MLF edits

fig = figure('Position',[0 00 800 1200],...
    'Color',[1 1 1]);


pos = .8;
nr  = 1 ;
d   = data_odp_967_22;
namestr = "ODP 967 (Grant et al., 2022)";
call_plot_200ka

pos = .7;
nr  = 2 ;
d   = data_kl09;
namestr = "KL 09";
call_plot_200ka

pos = .6;
nr  = 3 ;
d   = data_kl11;
namestr = "KL 11";
call_plot_200ka

pos = .5;
nr  = 4 ;
d   = data_kl15;
namestr = "KL 15";
call_plot_200ka

pos = .4;
nr  = 5 ;
d   = data_lake_tana;
namestr = "Lake Tana";
call_plot_200ka

pos = .3;
nr  = 6 ;
d   = data_icdp_chb;
namestr = "Chew Bahir";
call_plot_200ka

pos = .2;
nr  = 7 ;
d   = data_odp721_722_terr;
namestr = "ODP 721/722";
call_plot_200ka

pos = .1;
nr  = 8 ;
d   = data_odp_709;
namestr = "ODP 709";
call_plot_200ka


b(nr+1) = axes('Position',[0.1 pos-.05 0.8 0.05],...
    'XGrid','Off',...
    'YAxisLocation','left', ...
    'YDir', 'normal', ...
    'Color','None', ...
    'Ylim',[0 1], ...
    'YColor','None'); % Set YColor to 'None' to hide the y-axis

xl(10) = xlabel('Time (kyrs BP)');

l(9) = line(0,0);

linkaxes(b,"x")
b(1).XLim = [0 200];

z(1) = zoom(fig);
set(z(1),'Motion','horizontal');

p(1) = pan(fig);
set(p(1),'Motion','horizontal');

hold off
