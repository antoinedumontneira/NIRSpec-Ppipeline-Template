### LOAD NECESARY PACKAGES #####
import numpy as np
from glob import glob
import shutil
from pathlib import Path
import json
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import vstack
from show_image import show_image
import jwst
import traceback
import configparser



import traceback
import configparser
def mk_stpipe_log_cfg(out_dir,log_name):
    """
    Create a log file with the name log_name, where
    the pipeline will write all output.
    Args:
        log_name: str, name of the log to record screen output
    Returns:
        nothing
    """
    config = configparser.ConfigParser()
    config.add_section("*")
    config.set("*", "handler", "file:" + out_dir +log_name)
    config.set("*", "level", "INFO")
    pipe_log_config = os.path.join(out_dir,"pipeline-log.cfg")
    config.write(open(pipe_log_config, "w"))

#### STEP 1 PIPEKLINE 
from jwst.pipeline import calwebb_detector1

def mk_config_det1_before_snowblind(out_dir,uncal_file):
    ### CREATE AN INSTANCE OF THE PIPELINE TO MAKE A SPECIFIC CONFIG FILE 
    print("######################### CREATING CONFIGURATION FILE ###############################")
    config = calwebb_detector1.Detector1Pipeline.get_config_from_reference(uncal_file)
    detector1 = calwebb_detector1.Detector1Pipeline.from_config_section(config) 
    # Set some parameters that pertain to the entire pipeline
    detector1.output_dir = out_dir  ### We are still saving this temporary files in the working output directory
    detector1.save_results = False  ### NOT SAVING THE RESULTS AS WE ARE ONLY RUNNNG UNTIL JUMP STEP
    # Set some parameters that pertain to some of the individual steps
    detector1.jump.save_results = True
    # ALL THESE STEPS ARE SKIP, AS THEY COME AFTER JUMP DETECTION
    detector1.ramp_fit.skip = True 
    detector1.gain_scale.skip= True
    ##### TURNS OFF SNOWBALL FLAGGING IN THE PIPELINE
    detector1.jump.expand_large_events = False ##### TURNS OFF SNOWBALL FLAGGING IN THE PIPELINE
    detector1.export_config(out_dir+"config_detector1_before_snowblind.asdf")
    print("######################### CONFIGURATION FILE SAVED ###############################")

def PIPELINE_DETECTOR1_BEFORE_SNOWBLIND(uncal_file,out_dir):
    """
    Run the Detector1 pipeline on the given file.
    Args:
        uncal_file: str, name of uncalibrated file to run
        out_dir: str, path of the output directory
    Returns:
        nothing
    """
    log_name = os.path.basename(uncal_file).replace('.fits', '')
    mk_stpipe_log_cfg(out_dir,log_name+'.log')
    detector1 = calwebb_detector1.Detector1Pipeline()
    config_file = out_dir + "config_detector1_before_snowblind.asdf"
    logcfg = out_dir + "pipeline-log.cfg"
    detector1.call(uncal_file,config_file=config_file,logcfg=logcfg)

### RUN SNOWBLIND    
from snowblind import SnowblindStep
def run_snowblind(out_dir,jump_file):
    """ Remove large comic hits called snoballs of the jump files. 

    Args:
        out_dir (str): output directory. Same as where the jump files are located, 
        the rest of the pipeline is saving the outputs
        jump_file (fits):  jump file
    """
    snowblind_instance = SnowblindStep()
    snowblind_instance.output_dir = out_dir
    snowblind_instance.save_results = True
    snowblind_instance.suffix =  "snowblind"
    snowblind_instance.min_radius = 3 ### DEFAULT VALUE IS 4
    snowblind_instance.growth_factor = 1.5 ## DEFAULT VALUE IS 2
    run_snowblind = snowblind_instance.run(jump_file)
    ### REMOVE TEMPORARY JUMPs FILES OUTPUT DIRECTORY ################
    os.remove(jump_file)

def mk_config_det1_after_snowblind(out_dir,uncal_file):
    ### CREATE AN INSTANCE OF THE PIPELINE TO MAKE A SPECIFIC CONFIG FILE 
    print("######################### CREATING CONFIGURATION FILE ###############################")
    config = calwebb_detector1.Detector1Pipeline.get_config_from_reference(uncal_file)
    detector1 = calwebb_detector1.Detector1Pipeline.from_config_section(config) 
    detector1.output_dir = out_dir  ### We are still saving this temporary files in the working output directory
    detector1.save_results = True  ### 
    # Set some parameters that pertain to some of the individual steps
    # ALL THESE STEPS ARE SKIP AS THEY COME BEFORE RAMP FITTING
    detector1.group_scale.skip = True
    detector1.dq_init.skip = True
    detector1.saturation.skip = True
    detector1.ipc.skip = True
    detector1.superbias.skip = True
    detector1.refpix.skip = True
    detector1.rscd.skip = True
    detector1.firstframe.skip = True
    detector1.lastframe.skip = True
    detector1.linearity.skip = True
    detector1.dark_current.skip = True
    detector1.reset.skip = True
    detector1.persistence.skip = True
    detector1.jump.skip = True
    detector1.export_config(out_dir+"config_detector1_after_snowblind.asdf")
    print("######################### CONFIGURATION FILE SAVED ###############################")

def PIPELINE_DETECTOR1_AFTER_SNOWBLIND(jump_file,out_dir):
    """
    Run the Detector1 pipeline on the given file.
    Args:
        uncal_file: str, name of uncalibrated file to run
        out_dir: str, path of the output directory
    Returns:
        nothing
    """
    log_name = os.path.basename(jump_file).replace('.fits', '')
    mk_stpipe_log_cfg(out_dir,log_name+'.log')
    detector1 = calwebb_detector1.Detector1Pipeline()
    config_file = out_dir + "config_detector1_after_snowblind.asdf"
    logcfg = out_dir + "pipeline-log.cfg"
    detector1.call(jump_file,config_file=config_file,logcfg=logcfg)

####### STEP 2 PIPELINE
from jwst.pipeline import Spec2Pipeline
def mk_config_det2(out_dir,asn_file):
    ### CREATE AN INSTANCE OF THE PIPELINE TO MAKE A SPECIFIC CONFIG FILE 
    print("######################### CREATING CONFIGURATION FILE ###############################")
    config2 = Spec2Pipeline.get_config_from_reference(asn_file)
    ### CREATE AN INSTANCE OF THE PIPELINE WITH THE SPECIFIC CONFIG FILE 
    spec2 = Spec2Pipeline.from_config_section(config2) 
    spec2.save_results = True
    spec2.bkg_subtract.skip = False
    spec2.nsclean.skip = False
    #output_folder = os.mkdir(out_dir+'/'+asn_data["products"][0]['name'])
    spec2.output_dir = out_dir 
    spec2.input_dir = out_dir
    spec2.export_config(out_dir+"config_detector2.asdf")
    print("######################### CONFIGURATION FILE SAVED ###############################")

def PIPELINE_DETECTOR2(asn_file,out_dir):
    """
    Run the Detector2 pipeline on the given configuration file.
    Args:
        asn_file: str, name of association file to run
        out_dir: str, path of the output directory
    Returns:
        nothing
    """
    log_name = os.path.basename(asn_file).replace('.fits', '')
    mk_stpipe_log_cfg(out_dir,log_name+'.log')
    spec2 = Spec2Pipeline()
    config_file = out_dir + "config_detector2.asdf"
    logcfg = out_dir + "pipeline-log.cfg"
    spec2.call(asn_file,config_file=config_file,logcfg=logcfg)

####### STEP 3 PIPELINE
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline
def mk_config_det3(out_dir,asn_file):
    ### CREATE AN INSTANCE OF THE PIPELINE TO MAKE A SPECIFIC CONFIG FILE 
    print("######################### CREATING CONFIGURATION FILE ###############################")
    # Create an instance of the pipeline class
    config3 = Spec3Pipeline.get_config_from_reference(asn_file)
    spec3 = Spec3Pipeline.from_config_section(config3) 
    # Set some parameters that pertain to the entire pipeline
    spec3.output_dir = out_dir ############################## SAVING PATH ###################
    spec3.save_results = True
    # Set some parameters that pertain to some of the individual steps
    #spec3.extract_1d.bkg_fit = 'poly' # Fit a polynomial to the background values for each column or row
    spec3.cube_build.weighting='drizzle' # 'emsm'#'drizzle' 
    spec3.cube_build.single = False
    spec3.cube_build.coord_system= "ifualign"  ##### INSTRUMENT ALIGN INSTEAD OF SKY ALIGN
    spec3.export_config(out_dir+"config_detector3.asdf")
    print("######################### CONFIGURATION FILE SAVED ###############################")

def PIPELINE_DETECTOR3(asn_file,out_dir):
    """
    Run the Detector3 pipeline on the given configuration file.
    Args:
        asn_file: str, name of association file to run
        out_dir: str, path of the output directory
    Returns:
        nothing
    """
    log_name = os.path.basename(asn_file).replace('.fits', '')
    mk_stpipe_log_cfg(out_dir,log_name+'.log')
    spec3 = Spec3Pipeline()
    config_file = out_dir + "config_detector3.asdf"
    logcfg = out_dir + "pipeline-log.cfg"
    spec3.call(asn_file,config_file=config_file,logcfg=logcfg)