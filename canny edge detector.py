#code from https://stackoverflow.com/questions/7185655/applying-the-sobel-filter-using-scipy
import numpy
import scipy
import imageio
from scipy import ndimage

im = imageio.imread('edge_detection_test.jpg')
im = im.astype('int32')
dx = ndimage.sobel(im, 0)  # horizontal derivative
dy = ndimage.sobel(im, 1)  # vertical derivative
mag = numpy.hypot(dx, dy)  # magnitude
mag *= 255.0 / numpy.max(mag)  # normalize (Q&D)
imageio.imsave('sobel.jpg', mag)