#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

server.retrieve({
    "class": "ea",
    "dataset": "era5",
    "date": "2013-11-01/to/2014-03-31", # Time period
    "expver": "1",
    "levtype": "sfc",
    "param": "2t/sp",           # Parameters. Here we use 2m Temperature (2t) and Surface Pressure (sp). See the ECMWF parameter database, at http://apps.ecmwf.int/codes/grib/param-db
    "stream": "oper",
    "type": "an",
    "time": "00:00:00",
    "step": "0",
    "area": "-20/-50/-30/-40",    # Subset or clip to an area, here to Europe. Specify as North/West/South/East in Geographic lat/long degrees. Southern latitudes and Western longitudes must be given as negative numbers.
    "grid": "0.25/0.25",        # Regrid from the default grid to a regular lat/lon with specified resolution. The first number is east-west resolution (longitude) and the second is north-south (latitude).
    "format": "netcdf",         # Convert the output file from the default GRIB format to NetCDF format. Requires "grid" to be set to a regular lat/lon grid.
    "target": "summer2014.nc",    # The output file name. Set this to whatever you like.
})