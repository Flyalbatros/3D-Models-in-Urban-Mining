from __future__ import division
import os
import sys

import gdal, gdalconst
from gdalconst import *

#crs conversion
#from pyproj import Proj, transform


class envi_file(object):

    def __init__(self, file_name):

        '''opens a envi file so that operations can be performed with it'''
        '''code from: https://chris35wills.github.io/python-gdal-raster-io/'''
        #driver.Register()
        self.file_name=file_name
        self.inDs = gdal.Open(file_name, GA_ReadOnly)

        if self.inDs is None:
            print(
            "Couldn't open this file: " + file_name)
            print(
            '\nPerhaps you need an ENVI .hdr file?')
            sys.exit("Try again!")
        else:
            print(
            "%s opened successfully" % file_name)

            print(
            '~~~~~~~~~~~~~~'
            ,'Get image size'
            ,'~~~~~~~~~~~~~~')

            self.cols = self.inDs.RasterXSize
            self.rows = self.inDs.RasterYSize
            self.bands = self.inDs.RasterCount

            print ("columns: %i" % self.cols)
            print ("rows: %i" % self.rows)
            print ("bands: %i" % self.bands)

            print('~~~~~~~~~~~~~~')
            print('Get georeference information')
            print('~~~~~~~~~~~~~~')
            self.geotransform = self.inDs.GetGeoTransform()
            originX = self.geotransform[0]
            originY = self.geotransform[3]
            self.pixelWidth = self.geotransform[1]
            self.pixelHeight = self.geotransform[5]

            print("origin x: %i" % originX)
            print("origin y: %i" % originY)
            print("width: %2.2f" % self.pixelWidth)
            print("height: %2.2f" % self.pixelHeight)

            self.band_x = self.inDs.GetRasterBand(1).ReadAsArray(0, 0, self.cols, self.rows)
            self.band_y = self.inDs.GetRasterBand(2).ReadAsArray(0, 0, self.cols, self.rows)

    def get_array(self, band):
        '''code from: https://chris35wills.github.io/python-gdal-raster-io/'''
        # Set pixel offset.....
        print('~~~~~~~~~~~~~~')
        print('Convert image to 2D array')
        print('~~~~~~~~~~~~~~')
        band = self.inDs.GetRasterBand(band)
        self.image_array = band.ReadAsArray(0, 0, self.cols, self.rows)
        image_array_name = self.file_name
        print(type(self.image_array))
        print(self.image_array.shape)
        return self.image_array
        #self.pixelWidth, (self.geotransform, self.inDs)

    def crop_points(self, wkt_region_file, output_filenm):
        from shapely import wkt as wkt
        from shapely.geometry import Point as Point
        bound_file = open(wkt_region_file, 'r')
        wkt_region_geom = bound_file.readlines()
        bound_file.close()
        print(wkt_region_geom)
        bounding_poly = wkt.loads(wkt_region_geom[0])
        poly_bbox = bounding_poly.bounds
        out_file = open(output_filenm, 'w')
        for row in range(0, self.rows):
            for col in range(0, self.cols):
                coords_center = (self.band_x[row][col], self.band_y[row][col])
                if coords_center[0]>poly_bbox[0] and coords_center[0]<poly_bbox[2] and coords_center[1]>poly_bbox[1] and coords_center[0]<poly_bbox[3]:
                    if bounding_poly.contains(Point(coords_center)):
                        out_file.write(str(row)+','+str(col)+'\n')
                        print("wrote point")
                else:
                    print("rejected point")
        out_file.close()

    def selection_to_pixels(self, input_file, output_file):
        from shapely.geometry import LineString as LineString
        from shapely.geometry import MultiPoint as MultiPoint
        from shapely.geometry import Point as Point
        input = open(input_file, 'r')
        selected_points = input.readlines()
        input.close()
        output = open(output_file, 'w')
        output.write("geometry\n")
        output.close()
        for point_str in selected_points:
            index = point_str.replace('\n','').split(',')
            row = int(index[0])
            col = int(index[1])
            if row==0 or row==self.rows-1 or col==0 or col==self.cols-1:
                #print("end of line")
                continue
            center_coords = (self.band_x[row][col], self.band_y[row][col])
            diagonal_neighbor_coords = [(self.band_x[row-1][col-1], self.band_y[row-1][col-1]), (self.band_x[row+1][col-1], self.band_y[row+1][col-1]), (self.band_x[row+1][col+1], self.band_y[row+1][col+1]), (self.band_x[row-1][col+1], self.band_y[row-1][col+1])]
            pixel_border_pt_coords = []
            for neigbor_coords in diagonal_neighbor_coords:
                mean_coords = ((neigbor_coords[0]+center_coords[0])/2.0, (neigbor_coords[1]+center_coords[1])/2.0)
                pixel_border_pt_coords.append(mean_coords)
            pixel_border_pt_coords.append(pixel_border_pt_coords[0])
            output = open(output_file, 'a')
            output.write(str(LineString(pixel_border_pt_coords))+'\n')
            #print(pixel_border_pt_coords)
            #output.write(str(MultiPoint(diagonal_neighbor_coords))+'\n')
            #output.write(str(Point(center_coords))+'\n')
        output.close()

def to_xy_file(band_x, band_y, output_filenm):
    '''writes the 5000 to 5500th line to a xyz file'''
    file = open(output_filenm, 'w')
    file.write('x y z\n')
    print(len(band_x))
    for line_index in range(5000, 5500):
        for col_index in range(0, len(band_x[0])):
            file.write(str(band_x[line_index][col_index]) + ' ' + str(band_y[line_index][col_index]) + ' ' + '0\n')
    file.close()

def select_to_xy_file(band_x, band_y, selection, output_filenm):
    '''writes the coordinates from the pixels of selection (list of tuples) to a xyz file'''
    file = open(output_filenm, 'w')
    file.write('x y z\n')
    print(len(band_x))
    for coord in selection:
        file.write(str(band_x[coord[1]][coord[0]]) + ' ' + str(band_y[coord[1]][coord[0]]) + ' ' + '0\n')
    file.close()



geocorr_file = envi_file(r"C:\Users\P.A. Ruben\Desktop\Master thesis\12 - coding space\APEX_data\MM097_ROTTE_140917_a031d_calibr_cube000_igm.bsq")
geocorr_x = geocorr_file.get_array(1)
geocorr_y = geocorr_file.get_array(2)
#print(geocorr_x, geocorr_y)
geocorr_file.crop_points('bbox_wkt.txt', 'apex_points_in_bbox_flight_3.txt')
geocorr_file.selection_to_pixels('apex_points_in_bbox_flight_3.txt', 'pixel_geometry_flight_3.txt')
#print(geocorr_x[4144][370],geocorr_y[4144][370])
#to_xy_file(geocorr_x, geocorr_y, 'xy_coord_wgs84.xyz')

#list_pts_luxor=[(228,4144),(369,4145),(376,4140),(379,4143),(371,4147),(288,4045),(317,4032),(298,4058),(271,4072),(463,4125),(475,4076),(475,4086),(449,4104),(479,4113),(552,4279),(657,4251),(327,4057),(331,4053),(332,4061),(359,4145),(388,4143),(381,4137)]
#select_to_xy_file(geocorr_x, geocorr_y, list_pts_luxor, 'xy_coord_wgs84.xyz')

