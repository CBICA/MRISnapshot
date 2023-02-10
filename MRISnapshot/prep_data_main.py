#!/usr/bin/env python
import numpy as np
import pandas as pd
import pathlib

## Set logger  ## FIXME to be updated
import logging
format='%(levelname)-8s [%(filename)s : %(lineno)d - %(funcName)20s()] %(message)s'
format='%(levelname)-8s %(message)s'
logging.basicConfig(level=logging.DEBUG, format = '\n' + format, datefmt='%Y-%m-%d:%H:%M:%S')
logger = logging.getLogger(__name__)

##logger.setLevel(logging.DEBUG)      ## While debugging
logger.setLevel(logging.INFO)    ## FIXME Debug comments will be removed in release version
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"


### Import modules
def prep_data(params):

    ## Find underlay images
    ll = list(pathlib.Path(params.in_dir).glob('*' + params.suff_ulay))

    logger.info(ll)
    
    
    
    

    
    
    
