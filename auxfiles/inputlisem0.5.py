# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 09:06:36 2021
attempt lisem input script
@author: vjetten
"""

from pcraster import *
from pcraster.framework import *
import os
from osgeo import gdal, gdalconst, osr
#import gdal
from owslib.wcs import WebCoverageService
import subprocess

workingDir = 'C:/data/India/Cauvery/Base/'
os.chdir(workingDir)

### ---------- class GetSoilGridsLayer ---------- ###

class GetSoilGridsLayer:
    "downbloading a SOILGRIDS layer from WCS service"
    def __init__(self, mask, ESPG="",s="", i=1, j = 1):
        self.mask = mask
        self.varname = s
        self.layer = i
        self.outlayer = j
        self.debug = Debug_
        self.ESPG = ESPG

        if self.layer == 1: ID='_0-5cm_mean'
        if self.layer == 2: ID='_5-15cm_mean'
        if self.layer == 3: ID='_15-30cm_mean'
        if self.layer == 4: ID='_30-60cm_mean'
        if self.layer == 5: ID='_60-100cm_mean'
        if self.layer == 6: ID='_100-200cm_mean'

        #if self.debug == 1:
        print("Processing layer "+str(self.outlayer)+": "+self.varname+ID)

        raster=gdal.Open(self.mask)
        ESPG = 'urn:ogc:def:crs:EPSG::{0}'.format(self.ESPG)
        wide = raster.RasterXSize
        high = raster.RasterYSize
        dx = raster.GetGeoTransform()[1]
        dy = raster.GetGeoTransform()[5]
        llx = raster.GetGeoTransform()[0]
        ury = raster.GetGeoTransform()[3]
        urx = llx+wide*dx
        lly = ury+dy*high
        bbox = [llx,lly,urx,ury]

        if self.debug == 1:
            print("Mask ESPG and bounding box:"+ESPG,llx,lly,urx,ury,dx,dy)

        if self.debug == 1:
            print("Open SOILGRIDS WCS")

        url = "http://maps.isric.org/mapserv?map=/map/{}.map".format(self.varname)
        wcs = WebCoverageService(url, version='1.0.0')
        # show some info:
        # cov_list = list(wcs.contents)
        # mean_covs = [k for k in wcs.contents.keys() if k.find("mean") != -1]
        # print(mean_covs)

        variable = self.varname+ID
        varout = self.varname+str(self.outlayer)
        outputnametif = "{0}.tif".format(varout)
        outputnamemap = "{0}.map".format(varout)
        #outputnametmp = '_temp_.tif'

        if self.debug == 1:
            print("Downloading "+variable)

        # get data as temp geotif and save to disk
        response = wcs.getCoverage(identifier=variable,crs=ESPG,bbox=bbox,
            resx=dx,resy=dx,format='GEOTIFF_INT16')
        with open(outputnametif, 'wb') as file:
             file.write(response.read())

        # warp to some interpolation
        src = gdal.Open(outputnametif, gdalconst.GA_ReadOnly)
        src_proj = src.GetProjection()
        src_geotrans = src.GetGeoTransform()
        dst = gdal.GetDriverByName('PCRaster').Create(outputnamemap, wide, high, 1,
                                   gdalconst.GDT_Float32,["PCRASTER_VALUESCALE=VS_SCALAR"])
        dst.SetGeoTransform( src_geotrans )
        dst.SetProjection( src_proj )
        gdal.ReprojectImage(src, dst, src_proj, src_proj, gdalconst.GRA_Bilinear)
        #gdalconst.GRA_Cubic)

        # brute force convert tif to map by calling pcrcalc !!!!
        # CMD = "pcrcalc.exe"
        # arg = outputnamemap+"="+outputnametif
        # arg = '{0}{1}.map={0}{1}.tif'.format(self.varname,str(self.outlayer))
        # subprocess.run([CMD,arg])

        dst = None
        src = None
        if self.debug == 1:
            print("Done.\n")

### ---------- class GetSoilGridsLayer ---------- ###


### ---------- class PedoTranfer() ---------- ###

class PedoTransfer(StaticModel):
    # creates infiltration input vars: Ksat 1,2; Thetas1,2; Thetai1,2; Psi1,2
    # from SOILGRID.ORG GTiff maps for texture, org matter, gravel and bulkdensity
    # Using Saxton And Rawls 2006 pedotransferfunctions
    # https://hrsl.ba.ars.usda.gov/SPAW/Index.htm
    def __init__(self,mask=0,layer=1,moisture=0.7,sBD=1350.0):
        StaticModel.__init__(self)
    def initial(self):
        standardBD = scalar(standardbulkdensity_)  # standard bulk dens assumed by saxton and rawls. High! 1350 would be better
        fractionmoisture = scalar(initmoisture_)   #inital moisture as fraction between porosity and field capacity
        x = layer_
        mask=mask_

        S1 = readmap("{0}{1}.map".format(SG_names_[0],str(x)))  # sand g/kg
        Si1 = readmap("{0}{1}.map".format(SG_names_[1],str(x))) # silt g/kg
        C1 = readmap("{0}{1}.map".format(SG_names_[2],str(x)))  # clay g/kg
        OC1 = readmap("{0}{1}.map".format(SG_names_[3],str(x)))  # organic carbon in dg/kg
        Gravel1 = readmap("{0}{1}.map".format(SG_names_[4],str(x))) # coarse fragments cm3/dm3,
        bd1 = readmap("{0}{1}.map".format(SG_names_[5],str(x)))   # bulk density in cg/m3

        #output map name strings
        om1 = "om{0}.map".format(str(x))             # organic matter in %
        WP1 = "wilting{0}.map".format(str(x))      	# wilting point moisture content
        FC1 = "fieldcap{0}.map".format(str(x))     	# field capacity moisture content
        PAW1 = "plantAVW{0}.map".format(str(x))    	# plant available water content
        Coh1 = "coh{0}.map".format(str(x))           # soil cohesion (kPa)
        K1 = "k{0}.map".format(str(x))  		        #USLE erodibility
        BD1 = "bulkdens{0}.map".format(str(x))       # bulk density in kg/m3

        Pore1 = "thetas{0}.map".format(str(x))   	#porosity (cm3/cm3)
        Ksat1 = "ksat{0}.map".format(str(x))      	#ksat in mm/h
        initmoist1 = "thetai{0}.map".format(str(x))  # inital moisture (cm3/cm3)
        psi1 = "psi{0}.map".format(str(x))  		    # suction with init moisture in cm, used in LISEM
        Densityfactor1 = "densfact{0}.map".format(str(x))


        print("Creating infil params layer "+str(x))

        S = S1/1000  # from g/kg to fraction
        C = C1/1000
        Si = Si1/1000
        OC = (OC1/10000)*100  # conversion OC from dg/kg to percentage
        OM = OC*1.73  #/2.0   #conversion org carbon to org matter factor 2

        mask = ifthen(S+C+Si > 0.01,mask) # assume areas where sum text is not 1 = MV

        OM = OM*mask
        report(OM, om1)
        Dens = 1.0
        if docover:
           Dens = (1-0.1*cover)  #scalar(1.0)
        # density factor is 1.0, but could be made lower for organic soils and higher for compacted urban areas.
        bdsg = bd1*10            #bulkdensity cg/m3 to kg/m3
        bdsg = ifthenelse(bd1 < 1,standardBD,bdsg) # replace areas with MV bdsg to standard BD
        Gravel = Gravel1/1000  # from cm3/dm3 (1000 cc in a liter)
        Densityfactor = bdsg/standardBD*Dens #(1-0.1*cover)
        report(Densityfactor,Densityfactor1)

        #scalar(1.0)  # range 0.9 to 1.15
        # calculated as the bulk density from soilgrids divided by some standard bd
        # multiple regression from excel

        # wilting point stuff
        M1500 =-0.024*S+0.487*C+0.006*OM+0.005*S*OM-0.013*C*OM+0.068*S*C+0.031  #W18)
        # =-0.024*F18+0.487*G18+0.006*H18+0.005*F18*H18-0.013*G18*H18+0.068*F18*G18+0.031
        M1500adj =M1500+0.14*M1500-0.02  #X18) =W18+0.14*W18-0.02
        # field capacity stuff
        M33  =-0.251*S+0.195*C+0.011*OM+0.006*S*OM-0.027*C*OM+0.452*S*C+0.299  #Y18)
        #=-0.251*F18+0.195*G18+0.011*H18+0.006*F18*H18-0.027*G18*H18+0.452*F18*G18+0.299
        M33adj = M33+(1.283*M33*M33-0.374*M33-0.015)  #Z18) =Y18+(1.283*Y18*Y18-0.374*Y18-0.015)
        # porosity - FC
        PM33    = 0.278*S+0.034*C+0.022*OM-0.018*S*OM-0.027*C*OM-0.584*S*C+0.078  #AA18)
        #=0.278*F18+0.034*G18+0.022*H18-0.018*F18*H18-0.027*G18*H18-0.584*F18*G18+0.078
        PM33adj = PM33+(0.636*PM33-0.107)  #AB18) =AA18+(0.636*AA18-0.107)
        # porosity
        SatPM33 = M33adj + PM33adj  #AC18) =AB18+Z18
        SatSadj = -0.097*S+0.043  #AD18) =-0.097*F18+0.043
        SadjSat = SatPM33  + SatSadj  #AE18) =AC18+AD18
        Dens_om = (1-SadjSat)*2.65  #AF18) =(1-AE18)*2.65
        Dens_comp = Dens_om * Densityfactor  #AG18) =AF18*(I18)
        PORE_comp =(1-Dens_om/2.65)-(1-Dens_comp/2.65)  #AI18) =(1-AG18/2.65)-(1-AF18/2.65)
        M33comp = M33adj + 0.2*PORE_comp  #AJ18)  =Z18+0.2*AI18

        #output maps
        POROSITY = (1-Dens_comp/2.65)*mask  #AH18)
        PoreMcomp = POROSITY-M33comp  #AK18)
        LAMBDA = (ln(M33comp)-ln(M1500adj))/(ln(1500)-ln(33))  #AL18)
        GravelRedKsat =(1-Gravel)/(1-Gravel*(1-1.5*(Dens_comp/2.65)))  #AM18)

        Ksat = mask*max(0.0, 1930*(PoreMcomp)**(3-LAMBDA)*GravelRedKsat)  #AN18)
        BD = Gravel*2.65+(1-Gravel)*Dens_comp* mask     #U18
        WP = M1500adj*mask
        FC = M33adj* mask
        PAW = (M33adj - M1500adj)*(1-Gravel)* mask
        initmoist = fractionmoisture*POROSITY+ (1-fractionmoisture)*FC

        report(POROSITY,Pore1)
        report(Ksat,Ksat1)
        report(BD,BD1)
        report(WP,WP1)
        report(FC,FC1)
        report(PAW,PAW1)
        report(initmoist,initmoist1)

        # A = exp[ln(33) + B ln(T33)]
        # B = [ln(1500) - ln(33)] / [ln(T33) - ln(T1500)]
        bB = (ln(1500)-ln(33))/(ln(FC)-ln(WP))
        aA = exp(ln(33)+bB*ln(FC))
        Psi1= aA * initmoist**-bB *100/9.8
        report(Psi1,psi1)
        report(initmoist/POROSITY,"se1.map")

        ############################################
        ## Generation of soil cohesion map        ##
        ## estimation based on clay content       ##
        ## adapted from Morgan, 2001              ##
        ############################################

        Coh = max(1.0, 4.316*ln(C+1.0) - 6.955)
        # log fit using values below

        #Coh = lookupscalar("claycoh.tbl",C)*mask
        # content of claycoh.tbl
        # [,20>   2
        # [20,35> 3
        # [35,40> 9
        # [40,55> 10
        # [55,60> 11
        # [60,100> 12

        report(Coh, Coh1)

### ---------- class PedoTranfer() ---------- ###

workingDir = 'C:/data/India/Cauvery/Base/'
os.chdir(workingDir)

Debug_ = False #True

# tif file for projection ESPG id
maskname_ = 'dem200m.tif'
mask = gdal.Open(maskname_, gdalconst.GA_ReadOnly)
proj = osr.SpatialReference(wkt=mask.GetProjection())
ESPG = proj.GetAttrValue('AUTHORITY',1)

# pcraster file for exact basin mask
masknamemap_ = 'maskbasin.map'
#maskmap = gdal.Open(masknamemap_, gdalconst.GA_ReadOnly)


print("downloading SOILGRIDS layers...")
SG_names_ = ['sand','silt','clay','soc','cfvo','bdod']
# # texture, doil organic carbon, course fragments, bulk dens
# for x in range(0,6):
#     GetSoilGridsLayer(masknamemap_,ESPG,SG_names_[x],2,1)
# for x in range(0,6):
#     GetSoilGridsLayer(masknamemap_,ESPG,SG_names_[x],4,2)

setclone(masknamemap_)

mask_ = readmap(masknamemap_)
initmoisture_ = 0.7
standardbulkdensity_ = 1450.0
docover = False
cover_ = scalar(0) #readmap("per.map")
# optional, can use for bulkdensity, higher cover is more structure, lower density

PTF = PedoTransfer()
staticModel = StaticFramework(PTF)
layer_ = 1
staticModel.run()
layer_ = 2
staticModel.run()

print("Done")