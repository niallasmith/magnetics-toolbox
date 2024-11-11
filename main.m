close all
clear

plot_flag = 0; % 1 for plot data, 0 for no

% MAGNETICS TOOLBOX
% Niall Smith 
% October 2024 - November 2024

% to do:
% - add ui element allowing user to select core material, view plots etc.
% - interpolate with given granularity etc.
% - validate toroid shape dimensions, core vol etc.
% - calculate inductance with bode 100: just use amplitude permitivity at 0 flux,
%       should be equal
% - calculate magnetising current, amount of slope using curr sense
% - rename variables to make more sense / be clearer
% - temperature/flux/frequency out of bounds error
% - add term of magnetising current to primary current
% - mean length turn calcs: include layer number to increase MLT with more
%       layers
% - title etc for figure 2

% in progress:
% - add new data for different materials, 3C94 etc.
% - generalise code to accept different array sizes, values, etc.
% - make new repo on github and push
% - magnetising current

% done:
% - add flux density calculations based on voltages etc.
% - iterative process for calculating core operating temp
% - fix issue with specific power loss
% - plot power loss vs temperature for core selected
% - figure 2 plot title etc.
% - calculate duty ratio / effective phase shift based on pri, sec turns,
%       RDSon, pri current, 
% - add alternative data input for core volume, area, path length for
%       custom magnetics (not relying on calculated dimensions)
% - think about adding in different core shapes?
% - core size lookup table? 
% - make functions e.g. interpolate mu a vs flux and temp, etc.
% - look into epcos data, is this different? can we import epcos materials?
% - copper losses in power loss / core temp calculations
% - add plot on/off flag to functions

% core shapes
% ETD series:
% ETD34 ETD39 ETD44 ETD49 ETD54 ETD59
%
% EC identical (basically) to ETD, so not included
% ER series very similar to ETD, so not included
%
% PQ series:
% PQ20/16 PQ20/20 
% PQ26/20 PQ26/25 
% PQ32/15 PQ32/20 PQ32/25 PQ32/30 PQ32/35
% PQ35/20 PQ35/30 PQ35/35 
% PQ40/30 PQ40/40 
% PQ50/30 PQ50/35 PQ50/50
%
% path length // minimum area // eff area // volume
% etd 29/16/10  72 mm       71 mm^2     76 mm^2     5470 mm^3
% etd 34/17/11  78.6 mm     91.6mm^2    97.1 mm^2   7640 mm^3
% etd 39/20/13  92.2 mm     123 mm^2    125 mm^2    11500 mm^3
% etd 44/22/15  103 mm      172 mm^2    173 mm^2    17800 mm^3
% etd 49/25/16  114 mm      209 mm^2    211 mm^2    24000 mm^3
% etd 54/28/19  127 mm      270 mm^2    280 mm^2    35500 mm^3
% etd 59/31/22  139 mm      360 mm^2    368 mm^2    51500 mm^3
%
% pq20/16       37.6 mm     59.1 mm^2   61.9 mm^2   2330 mm^3
% pq20/20       45.7 mm     59.1 mm^2   62.6 mm^2   2850 mm^3
%
% pq26/20       45 mm       109 mm^2    121 mm^2    5470 mm^3
% pq26/25       54.3 mm     108 mm^2    120 mm^2    6530 mm^3
%
% pq32/15       44.7 mm     143 mm^2    171 mm^2    7640 mm^3
% pq32/20       55.9 mm     142 mm^2    169 mm^2    9440 mm^3
% pq32/25       64.2 mm     143 mm^2    168 mm^2    10800 mm^3
% pq32/30       74.7 mm     142 mm^2    167 mm^2    12500 mm^3
% pq32/35       83.5 mm     143 mm^2    166 mm^2    13800 mm^3
%
% pq35/20       60.5 mm     163 mm^2    183 mm^2    11100 mm^3
% pq35/30       77 mm       163 mm^2    192 mm^2    14800 mm^3
% pq35/35       86.1 mm     162 mm^2    190 mm^2    16300 mm^3
%
% pq40/30       81 mm       174 mm^2    212 mm^2    17100 mm^3
% pq40/40       102 mm      175 mm^2    201 mm^2    20500 mm^3
%
% pq50/30       74.4 mm     319 mm^2    329 mm^2    24500 mm^3
% pq50/35       84.4 mm     314 mm^2    329 mm^2    27800 mm^3
% pq50/50       113 mm      314 mm^2    328 mm^2    37100 mm^3

% available database cores:
% etd29/16/10
% etd34/17/11
% etd39/20/13
% etd44/22/15
% etd49/25/16
% etd54/28/19
% etd59/31/22 
%
% pq20/16
% pq20/20
% pq26/20
% pq26/25
% pq32/15
% pq32/20
% pq32/25
% pq32/30
% pq32/35
% pq35/20
% pq35/30
% pq35/35
% pq40/30
% pq40/40
% pq50/30
% pq50/35
% pq50/50

% available database materials
% Ferroxcube:
% - 3C90
% - 3C94
%
% EPCOS:
% - N87
% - N97


%
%       USER DATA
%

material = 'N97'; % material
vendor = 'epcos'; % material vendor (epcos, ferroxcube)
core_type = 'database'; % options: toroidal, user, database
core_shape = 'etd29/16/10'; % core shape. core_type must equal "database"


% physical dimensions of the toroidal core
% toroidal core size values core_type must equal "toroidal"
od = 42; % outer diameter in mm
id = 26; % inner diameter in mm
ht = 18; % height in mm


% user defined core size values core_type must equal "user"
user_core_area = 71*1e-6; % area in m^2
user_core_volume = 5470*1e-9; % volume in m^3
user_path_length = 72*1e-3; % length in m


pri_turns = 8; % primary turns of transformer
sec_turns = 7; % secondary turns of transformer
vin = 40; % input voltage
vout = 24; %output voltage
output_current = 14;
operating_frequency = 100000; % fundamental AC frequency of applied transformer voltage
pri_RDS = 10e-3; % Primary MOSFET RDSon
sec_RDS = 15e-3; % Secondary SR MOSFET RDSon
sec_DCR = 5e-3; % Secondary Inductor DCR

rb = 20; % burden resistor
nct = 200; % current transformer turns ratio


wire_diameter = 0.4; % diameter of wire used in mm
pri_num_strands = 5; % number of parallel strands in primary
sec_num_strands = 6; % number of parallel strands in secondary
pri_parallel_layers = 3;
sec_parallel_layers = 2; % secondary has an inherrent 2 parallel windings as each phase winding carries half the load current
% P S1 P S2 P effectively becomes P S P S P when considering total load current over one full cycle
strand_diameter = 1.5; % diameter of bundle of wires in mm

%
%       PRIMARY SCRIPT
%

% retrieve material and core data
[mu_r_matrix, pv_matrix_3d, mu_a_min_temp] = defineMaterialData(material,vendor,plot_flag);
[core_area_m2, core_volume_m3,avg_path_length_m] = defineCoreSize(core_type,od,id,ht,user_core_area,user_core_volume,user_path_length,core_shape);


% calculate skin depth from operating parameters
resistivity = 1.678*1e-8; % resistivity of copper at room temp
rel_perm = 0.999991; % relative permeability of copper at room temp
skin_depth = sqrt(resistivity/(pi()*operating_frequency*4*pi()*1e-7*rel_perm)) * 1e3;
min_wire_dia = 2 * skin_depth;


% calculate pri and sec resistance based on winding 
mean_length_turn = pi()*(11.8+strand_diameter)*1e-3; % mean length per turn in m
pri_wire_length = mean_length_turn * pri_turns * 1.05; % 5% more (ballpark), in m
sec_wire_length = mean_length_turn * sec_turns * 1.05;
pri_resistance_per_wire = resistivity * pri_wire_length / (pi()*(wire_diameter*1e-3)^2/4);
sec_resistance_per_wire = resistivity * sec_wire_length / (pi()*(wire_diameter*1e-3)^2/4);
pri_strand_resistance = 1000 * pri_resistance_per_wire / pri_num_strands; % strand resistance in mOhms
sec_strand_resistance = 1000 * sec_resistance_per_wire / sec_num_strands;
pri_resistance = pri_strand_resistance / pri_parallel_layers; % in mOhms
sec_resistance = sec_strand_resistance / sec_parallel_layers; % in mOhms


% calculate on time, used for calculating peak flux
duty_ratio = (vout/vin)*(pri_turns/sec_turns); % first pass estimate of duty ratio / effective phase shift

for i = 1:5 % use a iterative convergance to converge on true value of duty cycle
    sec_load_current = output_current/duty_ratio; % load current on the secondary when conducting
    pri_load_current = sec_load_current * (sec_turns/pri_turns);
    pri_voltage_drop = pri_load_current*(pri_RDS*2+pri_resistance*1e-3);
    sec_voltage_drop = sec_load_current*(sec_RDS+sec_resistance*1e-3+sec_DCR);
    duty_ratio = ((vout+sec_voltage_drop)/(vin-pri_voltage_drop))*(pri_turns/sec_turns); % duty ratio / effective phase shift - include RDSon voltage drop etc.
end

ton = duty_ratio/(2*operating_frequency); % time of applied voltage to the transformer
% what matters for peak flux is the time that there is an applied voltage.
% (see manitkala) for psfbc or fbc (or any bipolar magnetics topology i
% think), this 'on' time is at twice the switching frequency, e.g. the 
% maximum on time for a 100kHz psfbc is 5us, rather than 10us, as we have 
% 2 phases. This means the peak flux is halved as well. 


% calculate pri and sec copper losses in the transformer
pri_copper_losses = pri_load_current^2 * (pri_resistance*1e-3) * duty_ratio; % magnetising current
sec_copper_losses = sec_load_current^2 * (sec_resistance*1e-3) * duty_ratio;
copper_losses = pri_copper_losses + sec_copper_losses;


% calculate operating flux in the core
operating_flux = round(1000 * (vin*ton) / (pri_turns*core_area_m2)); % calculate peak flux density with applied voltage & frequency, number of primary turns
% rounded to the nearest whole mT


% plot magnetising inductance and power loss over temperature
Lmag_matrix_temps = 25:100; % range of temperatures to plot magnetising inductance over
[Lmag_matrix,Ploss_matrix] = plotLmagPlossOverTemperature(mu_r_matrix,pv_matrix_3d,operating_flux,operating_frequency,core_area_m2,avg_path_length_m,pri_turns,core_volume_m3,Lmag_matrix_temps,mu_a_min_temp,plot_flag);


% call function to calculate core temperature, using iterative convergance
core_temp = calculateCoreTemperature(pv_matrix_3d,operating_frequency,operating_flux,core_volume_m3,copper_losses,plot_flag);


% calculate core loss at operational values
core_loss = calculateCoreLoss(pv_matrix_3d,operating_frequency,operating_flux,core_temp,core_volume_m3);


if core_temp <= Lmag_matrix_temps(size(Lmag_matrix_temps,2)) % if core temperature does not exceed 100 degrees (within range of Lmag matrix)
    operating_lmag = Lmag_matrix((round(core_temp)-25)+1); % look up value of magnetising inductance at the operating temperature

    % calculate magnetising current within the core
    current_slope = 1e-6*vin/(Lmag_matrix((round(core_temp)-25)+1)*1e-6); % A/us
    meas_slope = current_slope * (rb / nct); % V/us
    peak_mag_current = current_slope * ton * 1e6;
    peak_meas = meas_slope * ton * 1e6;
end


% output data to the command window
disp('      MAGNETICS TOOLBOX V0.1      ')
disp(['Data for a ', core_shape, ' core of ', material, ' material'])
disp([num2str(pri_turns),' primary turns, ', num2str(sec_turns),' secondary turns'])
disp([num2str(vin),' V input, ' num2str(vout),' V output, ', num2str(output_current), 'A load',newline])
disp(['Maximum wire diameter: ', num2str(min_wire_dia), ' mm'])
disp(['Used wire diameter: ', num2str(wire_diameter), ' mm'])
disp(['Phase shift ratio (duty ratio) at max load: ', num2str(100 * duty_ratio,3), ' %'])
disp(['Peak flux density: ', num2str(operating_flux),' mT'])
disp(['Magnetising inductance at ',num2str(Lmag_matrix_temps(1)),' degC: ',num2str(Lmag_matrix(1),4),' uH'])
if exist('operating_lmag','var') == 1 % if operating lmag variable exists
    disp(['Magnetising inductance at ',num2str(core_temp,3),' degC: ',num2str(operating_lmag,4),' uH'])
end
if exist('meas_slope','var') == 1 % if operating lmag variable exists
    disp(['Magnetising inductance current sense slope at ',num2str(core_temp,3),' degC: ',num2str(1000*meas_slope,4),' kV/us'])
end
disp(['Core temperature: ', num2str(core_temp,3),' degC']) % display operating core temp
disp(['Core loss at operating temperature: ', num2str(core_loss), ' W'])
disp(['Copper loss at operating temperature: ', num2str(copper_losses), ' W'])
disp(['Total losses: ', num2str(core_loss+copper_losses),' W'])


function [mu_r_matrix,pv_matrix_3d,mu_a_min_temp] = defineMaterialData(material,vendor,plot_flag)

    if vendor == "ferroxcube"

        if material == "3C90" % 3C90 material data
    
            % mu a vs temperature and flux density data input
            mu_a_flux = [0, 50, 100, 150, 200, 250, 300, 350];
            mu_a_temps = [25, 100];
            mu_a_data = [2050, 3025, 3900, 4450, 4850, 4800, 4200, 3250;
                3400, 4200, 4900, 5250, 5300, 4900, 3600, 750];

            % pv vs flux density, frequency, and temperature
            pv_flux = [100, 200]; % known flux densities
            pv_freq = [25000, 50000, 100000, 200000]; % known frequencies to interpolate against
            pv_loss = [10,65;30,200;70,500;250,1700]; % known power losses for each frequency
            pv_temp = [0, 20, 40, 60, 80, 100, 120]; % temperature sample points for normalised pv
            norm_pv_temp = [1.62, 1.43, 1.25, 1.08, 0.99, 1.00, 1.16]; % normalised specific power loss to 100 degC at temperature sample points

        end

        if material == "3C91" % 3C91 material data - to do
        error('not yet implemented')
        end

        if material == "3C92" % 3C92 material data - to do
            error('not yet implemented')
        end

        if material == "3C94" % 3C94 material data

            % mu a vs temperature and flux density data input
            mu_a_flux = [0, 50, 100, 150, 200, 250, 300, 350, 400];
            mu_a_temps = [25, 100];
            mu_a_data = [2400, 2950, 3500, 4000, 4350, 4550, 4600, 4300, 3250;
                4000, 4200, 4350, 4500, 4600, 4500, 4250, 3100, 500];

            % pv vs flux density, frequency, and temperature
            pv_flux = [100, 200]; % known flux densities
            pv_freq = [25000, 100000, 200000]; % known frequencies to interpolate against
            pv_loss = [10,55;60,400;180,1100]; % known power losses for each frequency
            pv_temp = [0, 20, 40, 60, 80, 100, 120]; % temperature sample points for normalised pv
            pv_temp_val = [405,350,295,245,200,185,210]; % specific power loss data for fixed flux and frequency, over temperature
            norm_pv_temp = pv_temp_val/185; % normalised specific power loss to 100 degC at temperature sample points

        end

        if material == "3C95" % 3C95 material data - to do
            error('not yet implemented')
        end

        if material == "3C96" % 3C96 material data - to do
            error('not yet implemented')
        end

        if material == "3C97" % 3C97 material data - to do
            error('not yet implemented')
        end

        if material == "3C98" % 3C98 material data - to do
            error('not yet implemented')
        end

        mu_r_matrix = interpolate_mu_r(mu_a_flux,mu_a_temps,mu_a_data,plot_flag);

        pv_matrix_3d = ferroxcube_interpolate_pv(pv_flux,pv_freq,pv_loss,pv_temp,norm_pv_temp,plot_flag);

        mu_a_min_temp = mu_a_temps(1);
    end

    if vendor == "epcos"

        if material == "N87" % N87 material data
            % mu a vs temperature and flux density data input
            mu_a_flux = [0, 50, 100, 150, 200, 250, 300, 350, 400];
            mu_a_temps = [25,100];
            mu_a_data = [2100,2550,3030,3410,3660,3820,3760,3560,3190;
              4060,4110,4220,4360,4420,4370,4110,3450,2050];
            

            % pv vs flux density, frequency, and temperature
            pv_flux = [100, 200]; % known flux densities
            pv_loss_flux = [53,400]; % known power losses for each flux

            pv_freq = [50000, 100000, 200000]; % known frequencies to interpolate against
            pv_loss_freq = [18,51,171]; % known power losses for each frequency

            pv_temp = [40,60,80,100,120]; % temperature sample points for normalised pv
            pv_temp_val = [104,80,60,52,60]; % specific power loss data for fixed flux and frequency, over temperature
            norm_pv_temp = pv_temp_val/52; % normalised specific power loss to 100 degC at temperature sample points
        end

        if material == "N88" % N88 material data - to do
            error('not yet implemented')
        end

        if material == "N95" % N95 material data - to do
            error('not yet implemented')
        end

        if material == "N96" % N96 material data - to do
            error('not yet implemented')
        end

        if material == "N97" % N97 material data
            % mu a vs temperature and flux density data input
            mu_a_flux = [0, 50, 100, 150, 200, 250, 300, 350, 400];
            mu_a_temps = [25,100];
            mu_a_data = [2300,3880,4730,5160,5500,5720,5760,5580,5160;
              4100,5100,5670,5670,5670,5670,5520,5050,4160];

            % pv vs flux density, frequency, and temperature
            pv_flux = [100, 200]; % known flux densities
            pv_loss_flux = [43,301]; % known power losses for each flux

            pv_freq = [50000, 100000, 200000]; % known frequencies to interpolate against
            pv_loss_freq = [15,43,143]; % known power losses for each frequency

            pv_temp = [40,60,80,100,120,140]; % temperature sample points for normalised pv
            pv_temp_val = [110,80,59,43,48,69]; % specific power loss data for fixed flux and frequency, over temperature
            norm_pv_temp = pv_temp_val/43; % normalised specific power loss to 100 degC at temperature sample points
        end

        mu_r_matrix = interpolate_mu_r(mu_a_flux,mu_a_temps,mu_a_data,plot_flag);

        pv_matrix_3d = epcos_interpolate_pv(pv_flux,pv_loss_flux,pv_freq,pv_loss_freq,pv_temp,norm_pv_temp,plot_flag);

        mu_a_min_temp = mu_a_temps(1);
    end

    if exist('mu_a_flux','var') == 0
        error('error with material data')
    end

    if plot_flag == 1
        figure('Name','Datasheet relative amplitude permitivity','NumberTitle','off');
        plot(mu_a_flux,mu_a_data); % plot this data. should match datasheet.
        title('Relative amplitude permeability vs flux density, for different temperatures')
        xlabel('Flux density [mT]')
        ylabel('Relative Amplitude Permeability')
        legend('25 degC', '100 degC')
    end

end

function mu_r_matrix = interpolate_mu_r(mu_a_flux,mu_a_temps,mu_a_data,plot_flag)

    %
    % relative amplitude permeability vs temperature and flux density
    %

    mu_a_range = mu_a_temps(1):1:mu_a_temps(2); % define range of interpolation temperatures.
    % future, paramterise the interpolation granularity?

    mu_a_flux_range = 1:(max(mu_a_flux)+1);

    mu_a_matrix_temp = interp1(mu_a_temps,mu_a_data,mu_a_range); % linearly interpolate mu r vs temperatures
    mu_r_matrix = interp1(mu_a_flux,mu_a_matrix_temp',mu_a_flux_range); % interpolate mu r vs flux densities across known flux densities. 
    % this could probably be a 2D interpolation?

    if plot_flag == 1
        figure(2)
        surf(25:100,0:max(mu_a_flux),mu_r_matrix, 'EdgeColor', 'interp') % display surface of mu r across flux densities and temperatures
        % needs titles etc.
    end

end

function pv_matrix_3d = ferroxcube_interpolate_pv(pv_flux,pv_freq,pv_loss,pv_temp,norm_pv_temp,plot_flag)

    %
    % specific power loss for flux density, frequency and temperature
    %

    log_pv_flux = log10(pv_flux);
    log_pv_loss = log10(pv_loss);
    pv_range = 40:10:300; % range of flux densities to interpolate / extrapolate
    log_pv_range = log10(pv_range);
 
    log_pv_matrix = interp1(log_pv_flux,log_pv_loss',log_pv_range,'linear','extrap'); % interpolate and extrapolate power loss vs flux densities using given points in loglog space.
    pv_matrix = 10.^log_pv_matrix; % convert back to normal units

    if plot_flag == 1
        figure('Name','Specific Power Loss over Flux Density at Given Frequencies','NumberTitle','off');
        loglog(pv_range,pv_matrix) % plot power loss vs flux densities, for all 4 known freqencies (this should match datasheet)
        xlabel('Flux density [mT]')
        ylabel('Specific power loss [kW/m^3]')
        title('Specific power loss with flux density, for different frequencies')
        %legend('25 kHz','50 kHz','100 kHz', '200 kHz')
    end


    %pv_matrix_freq = zeros(36,17); % new matrix for storing data
    freq_range = 25000:5000:200000; % 25kHz to 200 kHz in 5kHz steps
 
    log_pv_matrix_freq = interp1(log10(pv_freq),log_pv_matrix',log10(freq_range)); % interpolate power loss vs frequency across known frequencies. 
    pv_matrix_freq = 10.^log_pv_matrix_freq; % convert back to normal units


    % figure(4)
    % loglog(freq_range,pv_matrix_freq) % plot power loss vs frequencies across all flux densities
    % xlabel('Frequency [Hz]')
    % ylabel('Specific power loss [kW/m^3]')
    % title('Specific power loss vs frequency, for different flux densities')
    % legend('40 mT','','','','','','','','','','','','','','','','200 mT') % there must be a better way of doing this

    if plot_flag == 1
        figure('Name','Specific Power loss over Flux Density and Frequency','NumberTitle','off');
        surf(pv_range,freq_range,pv_matrix_freq) % plot as surf
        set(gca,'zscale','log') % set all axes to be logarthmic
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        title('Specific power loss vs flux density and frequency')
        zlabel('Specific power loss [kW/m^3]')
        xlabel('Flux density [mT]')
        ylabel('Frequency [Hz]')
        set(gca,'ColorScale','log') % set log scale for color bar
    end


    pv_temp_range = 0:1:120; % range of temperatures to interpolate

    pv_matrix_temp = interp1(pv_temp,norm_pv_temp,pv_temp_range); % interpolate specific power loss at given temperature between given samples

    if plot_flag == 1

        figure('Name','Interpolated Normalised Specific Power Loss Over Temperature','NumberTitle','off');
        plot(pv_temp_range,pv_matrix_temp) % plot interpolated normalised specific power loss vs temperature
        xlabel('Temperature [degC]')
        ylabel('Normalised specific power loss')
        title('Normalised specific power loss vs temperature')
    end

    pv_matrix_3d = zeros(size(freq_range,2),size(pv_range,2),size(pv_temp_range,2)); % make new empty matrix to hold all the data in 4 dimensions (frequency, flux density, temperature, and specific power loss)

    for i = 1:121
        pv_matrix_3d(:,:,i) = pv_matrix_freq * pv_matrix_temp(i); % apply normalised specific power loss vs temperature mapping to 2d matrix of specific power loss vs flux density and frequency, based on temperature (array 3rd dimension index)
    end

    if plot_flag == 1
        figure('Name','Specific Power Loss at Slices of Flux Density, Temperature, and Frequency','NumberTitle','off');
        slice(pv_range,freq_range,pv_temp_range,pv_matrix_3d,[],[25000, 50000, 100000, 200000],25) % plot slices of the 3D matrix (4D data) at specific frequency and temperature slices
        set(gca,'xscale','log') % set log scales for flux density, frequency
        set(gca,'yscale','log')
        hcb=colorbar;
        title(hcb,'Specific Power Loss [kW/m^3]')
        set(hcb,'location','southoutside')
        set(gca,'ColorScale','log') % set log scale for solor bar
        xlabel('Flux density [mT]')
        ylabel('Frequency [Hz]')
        zlabel('Core temperature [degC]')
        title('Specific power loss vs flux density, frequency, and temperature')
    end

end

function pv_matrix_3d = epcos_interpolate_pv(pv_flux,pv_loss_flux,pv_freq,pv_loss_freq,pv_temp,norm_pv_temp,plot_flag)

    %
    % specific power loss for flux density, frequency and temperature
    %

    log_pv_flux = log10(pv_flux);
    log_pv_loss = log10(pv_loss_flux);
    pv_range = 40:10:300; % range of flux densities to interpolate / extrapolate
    log_pv_range = log10(pv_range);
 
    log_pv_matrix = interp1(log_pv_flux,log_pv_loss',log_pv_range,'linear','extrap'); % interpolate and extrapolate power loss vs flux densities using given points in loglog space.
    pv_matrix = 10.^log_pv_matrix; % convert back to normal units

    if plot_flag == 1
        figure('Name','Specific Power Loss over Flux Density at Given Frequencies','NumberTitle','off');
        loglog(pv_range,pv_matrix) % plot power loss vs flux densities (this should match datasheet)
        xlabel('Flux density [mT]')
        ylabel('Specific power loss [kW/m^3]')
        title('Specific power loss with flux density')
    end


    %pv_matrix_freq = zeros(36,17); % new matrix for storing data
    freq_range = 25000:5000:200000; % 25kHz to 200 kHz in 5kHz steps
 
    log_pv_matrix_freq = interp1(log10(pv_freq),log10(pv_loss_freq)',log10(freq_range),'linear','extrap'); % interpolate power loss vs frequency across known frequencies. 
    pv_matrix_freq = 10.^log_pv_matrix_freq; % convert back to normal units
    norm_pv_freq = pv_matrix_freq/pv_matrix_freq(((100000-25000)/5000)+1);

    pv_matrix_fb = zeros(size(pv_matrix,2),size(norm_pv_freq,2));

    for i = 1:size(norm_pv_freq,2)
        pv_matrix_fb(:,i) = pv_matrix * norm_pv_freq(i);
    end


    % figure(4)
    % loglog(freq_range,pv_matrix_freq) % plot power loss vs frequencies across all flux densities
    % xlabel('Frequency [Hz]')
    % ylabel('Specific power loss [kW/m^3]')
    % title('Specific power loss vs frequency, for different flux densities')
    % legend('40 mT','','','','','','','','','','','','','','','','200 mT') % there must be a better way of doing this

    if plot_flag == 1
        figure('Name','Specific Power loss over Flux Density and Frequency','NumberTitle','off');
        surf(freq_range,pv_range,pv_matrix_fb) % plot as surf
        set(gca,'zscale','log') % set all axes to be logarthmic
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        title('Specific power loss vs flux density and frequency')
        zlabel('Specific power loss [kW/m^3]')
        xlabel('Flux density [mT]')
        ylabel('Frequency [Hz]')
        set(gca,'ColorScale','log') % set log scale for color bar
    end


    pv_temp_range = 0:1:120; % range of temperatures to interpolate

    pv_matrix_temp = interp1(pv_temp,norm_pv_temp,pv_temp_range,'linear','extrap'); % interpolate specific power loss at given temperature between given samples

    if plot_flag == 1
        figure('Name','Interpolated Normalised Specific Power Loss Over Temperature','NumberTitle','off');
        plot(pv_temp_range,pv_matrix_temp) % plot interpolated normalised specific power loss vs temperature
        xlabel('Temperature [degC]')
        ylabel('Normalised specific power loss')
        title('Normalised specific power loss vs temperature')
    end

    pv_matrix_3d = zeros(size(pv_range,2),size(freq_range,2),size(pv_temp_range,2)); % make new empty matrix to hold all the data in 4 dimensions (frequency, flux density, temperature, and specific power loss)

    for i = 1:121
        pv_matrix_3d(:,:,i) = pv_matrix_fb * pv_matrix_temp(i); % apply normalised specific power loss vs temperature mapping to 2d matrix of specific power loss vs flux density and frequency, based on temperature (array 3rd dimension index)
    end

    if plot_flag == 1
        figure('Name','Specific Power Loss at Slices of Flux Density, Temperature, and Frequency','NumberTitle','off');
        slice(freq_range,pv_range,pv_temp_range,pv_matrix_3d,[25000, 50000, 100000, 200000],[],25) % plot slices of the 3D matrix (4D data) at specific frequency and temperature slices
        set(gca,'xscale','log') % set log scales for flux density, frequency
        set(gca,'yscale','log')
        hcb=colorbar;
        title(hcb,'Specific Power Loss [kW/m^3]')
        set(hcb,'location','southoutside')
        set(gca,'ColorScale','log') % set log scale for solor bar
        xlabel('Flux density [mT]')
        ylabel('Frequency [Hz]')
        zlabel('Core temperature [degC]')
        title('Specific power loss vs flux density, frequency, and temperature')
    end

end

function [core_area_m2, core_volume_m3, avg_path_length_m] = defineCoreSize(core_type,od,id,ht,core_area,core_volume,path_length,core_shape)

    if core_type == "toroidal"
        avg_dia = (od + id)/2;
    
        avg_path_length_mm = pi * avg_dia;
        avg_path_length_m = avg_path_length_mm * 1e-3;
        core_area_mm2 = (od/2-id/2)*ht;
        core_area_m2 = core_area_mm2 * 1e-6;

        core_volume_mm3 = (pi*od^2/4-pi*id^2/4)*ht;
        core_volume_m3 = core_volume_mm3 * 1e-9;
    end

    if core_type == "user"
        avg_path_length_m = path_length;
        core_area_m2 = core_area;
        core_volume_m3 = core_volume;
    end

    if core_type == "database"
        [core_area_m2, core_volume_m3, avg_path_length_m] = defineDatabaseCore(core_shape);
    end

end

function [core_area_m2, core_volume_m3, avg_path_length_m] = defineDatabaseCore(core_shape)

    %%%% ETD cores

    if core_shape == "etd29/16/10"
        avg_path_length_m = 72*1e-3;
        core_area_m2 = 76*1e-6;
        core_volume_m3 = 5470*1e-9;
    end

    if core_shape == "etd34/17/11"
        avg_path_length_m = 78.6*1e-3;
        core_area_m2 = 97.1*1e-6;
        core_volume_m3 = 7640*1e-9;
    end

    if core_shape == "etd39/20/13"
        avg_path_length_m = 92.2*1e-3;
        core_area_m2 = 125*1e-6;
        core_volume_m3 = 11500*1e-9;
    end

    if core_shape == "etd44/22/15"
        avg_path_length_m = 103*1e-3;
        core_area_m2 = 173*1e-6;
        core_volume_m3 = 17800*1e-9;
    end

    if core_shape == "etd49/25/16"
        avg_path_length_m = 114*1e-3;
        core_area_m2 = 211*1e-6;
        core_volume_m3 = 24000*1e-9;
    end

    if core_shape == "etd54/28/19"
        avg_path_length_m = 127*1e-3;
        core_area_m2 = 280*1e-6;
        core_volume_m3 = 35500*1e-9;
    end

    if core_shape == "etd59/31/22"
        avg_path_length_m = 139*1e-3;
        core_area_m2 = 368*1e-6;
        core_volume_m3 = 51500*1e-9;
    end

    %%%% PQ cores

    if core_shape == "pq20/16"
        avg_path_length_m = 37.6*1e-3;
        core_area_m2 = 61.9*1e-6;
        core_volume_m3 = 2330*1e-9;
    end

    if core_shape == "pq20/20"
        avg_path_length_m = 45.7*1e-3;
        core_area_m2 = 62.6*1e-6;
        core_volume_m3 = 2850*1e-9;
    end

    if core_shape == "pq26/20"
        avg_path_length_m = 45*1e-3;
        core_area_m2 = 121*1e-6;
        core_volume_m3 = 5470*1e-9;
    end

    if core_shape == "pq26/25"
        avg_path_length_m = 54.3*1e-3;
        core_area_m2 = 120*1e-6;
        core_volume_m3 = 6530*1e-9;
    end

    if core_shape == "pq32/15"
        avg_path_length_m = 44.7*1e-3;
        core_area_m2 = 171*1e-6;
        core_volume_m3 = 7640*1e-9;
    end

    if core_shape == "pq32/20"
        avg_path_length_m = 55.9*1e-3;
        core_area_m2 = 169*1e-6;
        core_volume_m3 = 9440*1e-9;
    end

    if core_shape == "pq32/25"
        avg_path_length_m = 64.2*1e-3;
        core_area_m2 = 168*1e-6;
        core_volume_m3 = 10800*1e-9;
    end

    if core_shape == "pq32/30"
        avg_path_length_m = 74.7*1e-3;
        core_area_m2 = 167*1e-6;
        core_volume_m3 = 12500*1e-9;
    end

    if core_shape == "pq32/35"
        avg_path_length_m = 83.5*1e-3;
        core_area_m2 = 166*1e-6;
        core_volume_m3 = 13800*1e-9;
    end

    if core_shape == "pq35/20"
        avg_path_length_m = 60.5*1e-3;
        core_area_m2 = 183*1e-6;
        core_volume_m3 = 11100*1e-9;
    end

    if core_shape == "pq35/30"
        avg_path_length_m = 77*1e-3;
        core_area_m2 = 192*1e-6;
        core_volume_m3 = 14800*1e-9;
    end

    if core_shape == "pq35/35"
        avg_path_length_m = 86.1*1e-3;
        core_area_m2 = 190*1e-6;
        core_volume_m3 = 16300*1e-9;
    end

    if core_shape == "pq40/30"
        avg_path_length_m = 81*1e-3;
        core_area_m2 = 212*1e-6;
        core_volume_m3 = 17100*1e-9;
    end

    if core_shape == "pq40/40"
        avg_path_length_m = 102*1e-3;
        core_area_m2 = 201*1e-6;
        core_volume_m3 = 20500*1e-9;
    end

    if core_shape == "pq50/30"
        avg_path_length_m = 74.4*1e-3;
        core_area_m2 = 329*1e-6;
        core_volume_m3 = 24500*1e-9;
    end

    if core_shape == "pq50/35"
        avg_path_length_m = 84.4*1e-3;
        core_area_m2 = 329*1e-6;
        core_volume_m3 = 27800*1e-9;
    end

    if core_shape == "pq50/50"
        avg_path_length_m = 113*1e-3;
        core_area_m2 = 328*1e-6;
        core_volume_m3 = 37100*1e-9;
    end

    if exist('avg_path_length_m','var') == 0 || exist('core_area_m2','var') == 0 || exist('core_volume_m3','var') == 0
        error('database core not defined')
    end

end

function Rth = calculateThermalResistance(core_volume) 
    % core volume in m^3
    % Rth in K/W

    Rth = 53*(core_volume*1E+6)^-0.54; % calculate thermal resistance of core using manitkala's formula
    % given for EE-EI-ETD-EC cores, so possibly not accurate

end

function [Lmag_matrix,P_loss_matrix] = plotLmagPlossOverTemperature(mu_r_matrix,pv_matrix_3d,operating_flux,operating_frequency,core_area,avg_path_length,pri_turns,core_volume,Lmag_matrix_temps,min_temps,plot_flag)

Lmag_matrix = zeros((max(Lmag_matrix_temps)-min(Lmag_matrix_temps)+1),1);

P_loss_matrix = zeros((max(Lmag_matrix_temps)-min(Lmag_matrix_temps)+1),1);

mu_0 = pi*4e-7; % mu0, permiability of free space

for i = min(Lmag_matrix_temps):max(Lmag_matrix_temps)
    mu_r = mu_r_matrix(operating_flux,i-min_temps+1);
    al = (mu_0*mu_r*core_area)/avg_path_length; % calculate core inductance factor
    
    L_magnetising = al * pri_turns^2; % calculate primary magnetising inductance
    L_magnetising_uh = L_magnetising * 1e6; % convert to uH
    Lmag_matrix(i-min(Lmag_matrix_temps)+1) = L_magnetising_uh;

    P_loss_specific = pv_matrix_3d(round(((operating_frequency-25000)/5000)),round(((operating_flux-40)/10)),round(i+1)); 
    P_loss_matrix(i-min(Lmag_matrix_temps)+1) = 1000 * P_loss_specific * core_volume; % calculate  core power loss at that operating temperature
end

if plot_flag == 1
    figure('Name','Magnetising Inductance Over Temperature','NumberTitle','off');
    plot(Lmag_matrix_temps,Lmag_matrix) % plot magnetising inductance with temperature
    title('Magnetising Inductance vs Temperature')
    xlabel('Temperature [degC]')
    ylabel('Magnetising Inductance [uH]')

    figure('Name','Power loss Over Temperature','NumberTitle','off');
    plot(Lmag_matrix_temps,P_loss_matrix) % plot power loss with temperature
    title('Core loss vs Temperature')
    xlabel('Temperature [degC]')
    ylabel('Core power loss [W]')
end

end

function operating_temp = calculateCoreTemperature(pv_matrix,frequency,flux_density,core_volume,copper_losses,plot_flag)
  %
    % Iterative process for calculating operational core temperature
    %

    new_operating_temp = 25; % start iterative loop with first pass temperature guess (ambient temperature, 25 degC)
    operating_temp_iterations = zeros(10,1);
    operating_temp_iterations(1) = new_operating_temp;

    Rth = calculateThermalResistance(core_volume);

    for i = 2:10 % 10 iterations - enough? probably
        if new_operating_temp >= 120
            error('Core temperature exceeds 120 degrees')
        end
        P_loss_specific = pv_matrix(round(((frequency-25000)/5000)),round(((flux_density-40)/10)),round(new_operating_temp+1)); % retrieve specific power loss based on core temperature
        P_loss = 1000* P_loss_specific * core_volume + copper_losses; % calculate core power loss in W
        core_temp = Rth * P_loss; % calculate core temp. does NOT include wire (copper) losses.
    
        new_operating_temp = (core_temp+new_operating_temp)*0.5; % new starting temperature for next iteration
        % go for the midpoint between the original starting point and the new
        % operating temperature to avoid oscillations/overshoots (does this
        % make sense?)
        operating_temp_iterations(i) = new_operating_temp; % add data to matrix to keep track of how core temmperature chnages with iterations
    end

    operating_temp = new_operating_temp;

    if plot_flag == 1
        figure('Name','Core Operating Temperature Iteration Convergance','NumberTitle','off');
        plot(operating_temp_iterations) % plot how core temperature estimate changes with iterations
    end

end

function core_loss = calculateCoreLoss(pv_matrix,frequency, flux_density, temperature, core_volume)
    P_loss_specific = pv_matrix(round(((frequency-25000)/5000)),round(((flux_density-40)/10)),round(temperature+1)); % retrieve specific power loss based on core temperature
    core_loss = 1000* P_loss_specific * core_volume; % calculate core power loss in W
end