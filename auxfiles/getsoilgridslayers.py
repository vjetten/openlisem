# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 10:35:23 2021

@author: vjetten

dowbnloading 6 layers from soilgrids with a mask tif
"""
import os
from osgeo import gdal, osr
from owslib.wcs import WebCoverageService


class GetSoilGridsLayer:
    "downbloading a SOILGRIDS layer from WCS service"
    def __init__(self, mask = "", s="", i=1, j = 1):
        self.mask = mask
        self.varname = s
        self.layer = i
        self.outlayer = j
        self.debug = Debug_

        if self.layer == 1: ID='_0-5cm_mean'
        if self.layer == 2: ID='_5-15cm_mean'
        if self.layer == 3: ID='_15-30cm_mean'
        if self.layer == 4: ID='_30-60cm_mean'
        if self.layer == 5: ID='_60-100cm_mean'
        if self.layer == 6: ID='_100-200cm_mean'

        if self.debug == 1:
            print("Processing layer "+str(self.layer)+": "+self.varname+ID)
        raster=gdal.Open(self.mask)
        proj = osr.SpatialReference(wkt=raster.GetProjection())
        ESPG = 'urn:ogc:def:crs:EPSG::{0}'.format(proj.GetAttrValue('AUTHORITY',1))
        wide = raster.RasterXSize
        high = raster.RasterYSize
        dx = raster.GetGeoTransform()[1]
        dy = raster.GetGeoTransform()[5]
        lx = raster.GetGeoTransform()[0]
        uy = raster.GetGeoTransform()[3]
        rx = lx+wide*dx
        ly = uy+dy*high
        if self.debug == 1: print("Mask ESPG and bounding box:"+ESPG,lx,ly,rx,uy)

        bbox = [lx,ly,rx,uy]

        if self.debug == 1: print("Open SOILGRIDS WCS")
        url = "http://maps.isric.org/mapserv?map=/map/{}.map".format(self.varname)
        wcs = WebCoverageService(url, version='1.0.0')
        # cov_list = list(wcs.contents)
        # mean_covs = [k for k in wcs.contents.keys() if k.find("mean") != -1]
        # print(mean_covs)

        variable = self.varname+ID
        if self.debug == 1: print("Downloading "+variable)
        response = wcs.getCoverage(
            identifier=variable,
            crs=ESPG,
            bbox=bbox,
            resx=dx,
            resy=dx,
            format='GEOTIFF_INT16')

        # Save the data as geotif
        varout = self.varname+str(self.outlayer)
        outputname = "{}.tif".format(varout)
        with open(outputname, 'wb') as file:
            file.write(response.read())
        if self.debug == 1: print("Done.\n")


workingDir = 'C:/data/India/Cauvery/Base/'
os.chdir(workingDir)
maskname_ = 'dem200m.tif'
Debug_ = True

print("downloading SOILGRIDS layers...")
GetSoilGridsLayer(maskname_,'clay',2,1)
GetSoilGridsLayer(maskname_,'sand',2,1)
GetSoilGridsLayer(maskname_,'silt',2,1)
GetSoilGridsLayer(maskname_,'soc',2,1)
GetSoilGridsLayer(maskname_,'bdod',2,1)
GetSoilGridsLayer(maskname_,'cfvo',2,1)

GetSoilGridsLayer(maskname_,'clay',4,2)
GetSoilGridsLayer(maskname_,'sand',4,2)
GetSoilGridsLayer(maskname_,'silt',4,2)
GetSoilGridsLayer(maskname_,'soc',4,2)
GetSoilGridsLayer(maskname_,'bdod',4,2)
GetSoilGridsLayer(maskname_,'cfvo',4,2)
print("Done")