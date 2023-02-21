#!/usr/bin/env python

### Import required modules
from PIL import Image, ImageFilter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

### return single image
def singleImage(img2d_under):
    
    # Create PIL image for overlay
    pil_under = Image.fromarray(np.uint8(cm.Greys_r(img2d_under, bytes=True)))
    return pil_under

### overlay an img on another
def overlayImage(img2d_under, img2d_over, transparency=0.6, is_edge=0):
    
    # Create PIL images for overlay and underlay images
    pil_under = Image.fromarray(np.uint8(cm.Greys_r(img2d_under, bytes=True)))
    pil_fused = Image.fromarray(np.uint8(cm.Greys_r(img2d_under, bytes=True)))
    
    ## tmpData =  np.uint8(np.tile(np.expand_dims(img2d_over,3),(1,1,4))*255)    ## Old version
    tmpData =  np.uint8(np.tile(np.expand_dims(img2d_over,2),(1,1,4))*255)
    tmpData[:,:,[1,2]] = 0        # Set RED color
    pil_over = Image.fromarray(tmpData)
    
    # Fing edges of overlay img
    if is_edge == 1:
        pil_over = pil_over.filter(ImageFilter.FIND_EDGES)
        pil_over = pil_over.point(lambda p: p > 0 and 255)        # to binarize edge image  

    # Set alpha=0 for the background (intensity=0) of the fg image
    red, green, blue, alpha = pil_under.split()

    if transparency < 1:
        pil_fused = Image.blend(pil_under, pil_over, transparency)
    else:
        pil_fused.paste(pil_over, (0,0), pil_over)
        pil_fused.putalpha(alpha)

    return pil_under, pil_fused
    #return pil_over, pil_fused
    
### overlay two images on another
def overlayImageDouble(img2d_under, img2d_over1, img2d_over2, transparency=0.6, is_edge=0):
    
    # Create PIL images for overlay and underlay images
    pil_under = Image.fromarray(np.uint8(cm.Greys_r(img2d_under, bytes=True)))
    pil_fused = Image.fromarray(np.uint8(cm.Greys_r(img2d_under, bytes=True)))
    
    ## tmpData1 =  np.uint8(np.tile(np.expand_dims(img2d_over1,3),(1,1,4))*255)    ## Old version
    tmpData1 =  np.uint8(np.tile(np.expand_dims(img2d_over1,2),(1,1,4))*255)
    tmpData1[:,:,[1,2]] = 0        # Set RED color

    ## tmpData2 =  np.uint8(np.tile(np.expand_dims(img2d_over2,3),(1,1,4))*255)    ## Old version
    tmpData2 =  np.uint8(np.tile(np.expand_dims(img2d_over2,2),(1,1,4))*255)
    tmpData2[:,:,[0,2]] = 0        # Set GREEN color

    pil_over1 = Image.fromarray(tmpData1)
    pil_over2 = Image.fromarray(tmpData2)
    
    # Fing edges of overlay img
    if is_edge == 1:
        pil_over1 = pil_over1.filter(ImageFilter.FIND_EDGES)
        pil_over1 = pil_over1.point(lambda p: p > 0 and 255)        # to binarize edge image  

        pil_over2 = pil_over2.filter(ImageFilter.FIND_EDGES)
        pil_over2 = pil_over2.point(lambda p: p > 0 and 255)        # to binarize edge image  

    # Set alpha=0 for the background (intensity=0) of the fg image
    red, green, blue, alpha = pil_under.split()

    transparency = float(transparency)
    if transparency < 1:
        pil_fused = Image.blend(pil_under, pil_over1, transparency)
        pil_fused = Image.blend(pil_fused, pil_over2, transparency)
    else:
        pil_fused.paste(pil_over1, (0,0), pil_over1)
        pil_fused.paste(pil_over2, (0,0), pil_over2)
        pil_fused.putalpha(alpha)

    return pil_under, pil_fused

