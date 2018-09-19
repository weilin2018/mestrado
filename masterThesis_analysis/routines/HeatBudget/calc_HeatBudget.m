%Jessica Benthuysen
%April 27, 2016
%calc_heatbudget calculates the terms in the heat budget from OceanCurrents
%using: dT/dt = - div_H dot (uT) + d/dz(K_V * dT/dz) + R,
%where R is the residual term. The temperature tendency equation is volume
%averaged (see Peter et al 2006, JGR for a mixed layer example).
%
%The code outputs the terms as follows:
%Volume averaged temperature tendency, dT/dt = int_adv_d (horizontal advection) + Qsnet (surface heat flux)

clear all

addpath(genpath('/home/ecoliver/matlab/netcdf_toolbox/'));
%addpath(genpath('/home/ecoliver/matlab/mexcdf/'));
addpath(genpath('/home/danilo/Downloads/delft/delft3d_repository/src/third_party_open/netcdf/matlab/mexnc/'))
addpath(genpath('/home/ecoliver/Desktop/include/'));
addpath(genpath('heat_budget/'));
rehash

% colocar diretorio para os dados do Mercator
header_MERCATOR = '/media/danilo/Danilo/mestrado/ventopcse/data/mercator/';
%header_OMAPS = '/data/MHWs/Tasmania_2015_2016/OMAPS/';
data_shf = 'CFSv2';
