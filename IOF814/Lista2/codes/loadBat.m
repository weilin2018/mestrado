function [ lon,lat,bat ] = loadBat( batFile )
%loadBat Summary of this function goes here
%   Load bathymetry, longitude and latitude data from .xyz file
%   created by ETOPO.

    % definition of constant values used in the routine
    theta=(latmax+latmin)/2;         % Mean Latitude
    Rt=6371e3;                       % Earth's radius
    Ct=2*pi*Rt;                      % Earth's lenght
    dlat=Ct/360;                     % latitude in meters
    dlon=dlat*cos(deg2rad(theta));   % longitude in meters
    
    % loading and extracting data
    file = load(batFile);
    lon = file(:,1);
    lat = file(:,2);
    bat = file(:,3);
    
    % extract max and min value
    lonmin=min(lon);
    lonmax=max(lon);
    latmin=min(lat);
    latmax=max(lat);
    batmin=min(bat);
    batmax=max(bat);
    
    % generating grid based on data extracted
    nlat = linspace(latmin, latmax, kmax);
    nlon = linspace(lonmin, lonmax, jmax);
    [llat, llon] = meshgrid(nlat, nlon); % grid points of the grid
    
    % interpolating bathymetry data for the new grid
    nbat=griddata(lon,lat,bat,llon,llat);
    
    

end

