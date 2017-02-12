% thermalFilm.m
% Classical thermal calc for Bismuth Selenide (Bi2Se3) film on sapphire substrate
% Presumes that the temperature rise on the film is uniform and instantaneous
% and that the thermal contact is perfect
% A more sophisticated model with better Physics (e.g. AMM, DMM) will be needed!
% Also, does not include nanoscale phenomena (e.g. phonon mean free path)
% Also, no acoustic propogation is included.  This could be called later.
% Model: see p. 429 (Sec. 10.7) of Hahn, "Heat Conduction" 3rd edition
% Material properties are hard coded
% Based on thermalFilm.m, First written by Eric Landahl, 12.28.2016
% Revised by EL 2.12.2017
%
%% INPUTS:
%   L           Bi2Se3 film thickness in meters (sapphire substrate is infinite)
%   fluence     absorbed laser fluence in mJ/cm^2
%   time        a vector of times to be calculated in seconds
%   max_depth   maximum depth calculated into sapphire, typ. 5 * extinction depth
%% OUTPUTS:
%   T1         temperature of film, size = length(time) x length(z1)
%   T2         temperature of substrate, size = length(time) x length(z2)
%   z1         a vector of film depths in meters
%   z2         a vector of substrate depths in meters
%
%% TYPICAL USAGE
% 
% [T1 T2 z1 z2] = Bi2Se3_thermal (150e-9, 1, (1e-10:2e-10:1e-8), 1e-5);
%
function [T1 T2 z1 z2]=Bi2Se3_thermal (L,fluence,time,max_depth)

% Bi2Se3 film properties.  "1" refers to the film
  C1 = 124.3*1000/654.8; %Specific heat of film in J/(kg K)
  rho1 = 6.82e3; % Film density in kg/m^3
  k1 = 0.75; % Film thermal conductivity in W/(m K)
  %alpha1 = 0.312e-4 % Film thermal diffusivity in m^2/s
  alpha1 = k1/(rho1 * C1); % Film thermal diffusivity in m^2/s
  
% Sahhpire substrate properties.  "2" referes to the substrate
  C2 = 761; % Specific heat of sapphire in J/(kg K)
  rho2 = 3.98e3; % Sapphire density in kg/m^3
  k2 = 23.1; % Sapphire thermal conductivity in W/(m K) parallel to optic axis
  alpha2 = k2/(rho2 * C2); % Sapphire thermal diffusivity in m^2/s
  
  %% Temporary for troubleshooting: make the sampe all semiconductor
%  rho1 = rho2;
%  k1 = k2;
%  C1 = C2;
%  alpha1 = alpha2;
[k1 k2]
[C1 C2]
[rho1 rho2]
 
% Calculate initial temperature rise
  fluence = fluence*10; % Convert from mJ/cm^2 to J/m^2
  T0 = fluence/(L * C1 * rho1); % Initial temperature rise in film
  fprintf('A %d nm thick Bi2Se3 film gives a temperature rise of %.1f K.\n',L*1e9,T0)
  
% Unitless parameters (see Hahn, "Thermal Conductivity", Eqs. 10-135 and 10-138)  
  mu = sqrt(alpha1/alpha2);
  beta = (k1/k2)/mu;
  gamma = (beta - 1)/(beta + 1);
  
% Spatial grid
  num_depths = 1000;  % number of depth points z to be calculated
  dz = max_depth/num_depths;
  z = dz:dz:max_depth;
  
% Meshgrid for calculation speed & ease
  [Time Z] = meshgrid(time,z); % Time and Z are 2D, time and z are 1D

% Calculate temperature profile in bulk
  max_n = 100; % number of terms in series expansion, default 100
  T2a = 0.*Time.*Z; % each term gets added to this, starts at zero
  for n = 0: max_n  % Series expansion solution of heat equation
    T2b = erfc((2*n*L + mu*Z)./(2*(sqrt(alpha1*Time)))); % temporary
    T2c = erfc(((2*n + 2)*L + mu*Z)./(2*sqrt(alpha1*Time))); % temporary
    T2a = T2a + (gamma^n) * (T2b - T2c); % temporary, adding up
  end
  T2 = T0 * (1/2) * (1 + gamma) * T2a; % Temperature at all z and time  
  
% Film calculations 
% Calculate temperature profile in film.  zz and ZZ are the film depths
  dzz = L/100; % Choose 100 depth points in the film by default  
  zz = dzz:dzz:L;
  [Time ZZ] = meshgrid(time,zz); 
  T1a = 0.*Time.*ZZ;
  for n = 0:max_n
    T1b = erfc(((2*n + 1)*L - ZZ)./(2*sqrt(alpha1*Time)));
    T1c = erfc(((2*n + 1)*L + ZZ)./(2*sqrt(alpha1*Time)));
    T1a = T1a + (gamma^n) * (T1b + T1c);
  end
  T1 = T0 - T0 * (1/2) * (1 - gamma) * T1a; % Temperature in film
  T1_end = mean(T1(:,end)); % Average temperature at final timepoint
  fprintf('After %.1f ns, the average film temperature rise is only %.1f deg C.\n', ...
  time(end)*1e9,T1_end);
  fprintf('The maximum bulk temperature rise is %.1f deg C.\n', max(max(T2)));
  
  z1 = zz; % for output
  z2 = z; % for output
  
% Calculate heat in film
% Q = integral dzz of rho*C*T
%  Q1 = trapz(ZZ,T1*rho1*C1);
%  Q2 = trapz(Z,T2*rho2*C2);
%  Q1 = Q1/10; % convert from J/m^2 to mJ/cm^2
%  Q2 = Q2/10; % convert from J/m^2 to mJ/cm^2
  
% Outputs not needed for TRXD

save thermalFilmOut.m; % Save all variables for future use
  
  end
  
  