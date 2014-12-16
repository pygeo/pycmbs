"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
test if polygons or lines interesect
"""


# http://pcjericks.github.io/py-gdalogr-cookbook/geometry.html#calculate-intersection-between-two-geometries

from osgeo import ogr

#~ print 'CASE 1'
#~ wkt1 = "POLYGON ((1208064.271243039 624154.6783778917, 1208064.271243039 601260.9785661874, 1231345.9998651114 601260.9785661874, 1231345.9998651114 624154.6783778917, 1208064.271243039 624154.6783778917))"
#~ wkt2 = "POLYGON ((1199915.6662253144 633079.3410163528, 1199915.6662253144 614453.958118695, 1219317.1067437078 614453.958118695, 1219317.1067437078 633079.3410163528, 1199915.6662253144 633079.3410163528)))"
#~
#~ poly1 = ogr.CreateGeometryFromWkt(wkt1)
#~ poly2 = ogr.CreateGeometryFromWkt(wkt2)
#~
#~ intersection = poly1.Intersection(poly2)
#~
#~ print intersection.ExportToWkt()

# practical example for dateline crossing problem
print ''
print 'CASE 2'

# create a polygon first
# Create ring
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(-10., 20.)
ring.AddPoint(10., 30.)
ring.AddPoint(15., 10.)
ring.AddPoint(-20., 5.)
ring.AddPoint(-10., 20.)

# Create polygon
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(ring)

print poly.ExportToWkt()


# create line
greenwich = ogr.Geometry(ogr.wkbLineString)
greenwich.AddPoint(0., 90.)
greenwich.AddPoint(0., -90.)
print greenwich.ExportToWkt()


#~ # greenwhich line as a polygon
#~ gring = ogr.Geometry(ogr.wkbLinearRing)
#~ gring.AddPoint(-0.00000001,  90.)
#~ gring.AddPoint( 0.00000001,  90.)
#~ gring.AddPoint( 0.00000001, -90.)
#~ gring.AddPoint(-0.00000001, -90.)
#~
#~ # Create polygon
#~ gp = ogr.Geometry(ogr.wkbPolygon)
#~ gp.AddGeometry(gring)


intersection1 = poly.Intersection(greenwich)
if intersection1 is not None:

    print '*** INTERSECTION FOUND!', intersection1.ExportToWkt()
else:
    'Intersection 1 indicates no intersection'

#~ intersection2 = poly.Intersection(gp)
#~ print intersection2.ExportToWkt()



ring1 = ogr.Geometry(ogr.wkbLinearRing)
ring1.AddPoint(10., 20.)
ring1.AddPoint(10., 30.)
ring1.AddPoint(5., 30.)
ring1.AddPoint(5., 20.)
ring1.AddPoint(10., 20.)

# Create polygon
poly1 = ogr.Geometry(ogr.wkbPolygon)
poly1.AddGeometry(ring1)

intersection2 = poly1.Intersection(greenwich)
print intersection2.ExportToWkt()



dateline = ogr.Geometry(ogr.wkbLineString)
dateline.AddPoint(180., 90.)
dateline.AddPoint(180., -90.)

intersection3 = poly1.Intersection(greenwich)
print intersection3.ExportToWkt()



# test for crossing dateline



dateline2 = ogr.Geometry(ogr.wkbLineString)
dateline2.AddPoint(180., 90.)
dateline2.AddPoint(180., -90.)



# ensure coordinates between 0 ... 360 by adding 360 to coordinates which are below zero

ring2 = ogr.Geometry(ogr.wkbLinearRing)
ring2.AddPoint(150., 30.)
ring2.AddPoint(-160+360., 30.)
ring2.AddPoint(-170+360., 20.)
ring2.AddPoint(170., 20.)
ring2.AddPoint(150., 30.)

poly2 = ogr.Geometry(ogr.wkbPolygon)
poly2.AddGeometry(ring2)


print poly2.ExportToWkt()
print dateline2.ExportToWkt()


print 'ref=False', poly2.Intersection(greenwich).ExportToWkt()
print 'ref=True', poly2.Intersection(dateline2).ExportToWkt()






#point in poly

print ''
print ''
print 'POINT IN POLYGON'

point1 = ogr.Geometry(ogr.wkbPoint)
point1.AddPoint(7., 25.)  # should fit in poly1
point2 = ogr.Geometry(ogr.wkbPoint)
point2.AddPoint(70., 25.)  # should NOT fit in poly1

print point1.ExportToWkt()
i1 = poly1.Intersection(point1)
i2 = poly1.Intersection(point2)
print 'ref=True', i1.ExportToWkt()
print 'ref=False', i2.ExportToWkt()






