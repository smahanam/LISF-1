#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.8
#
# Copyright (c) 2026 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
#------------------------------------------------------------------------------
#
# SCRIPT: generate_ldtconfig_files_ensrst_nrt.py
#
# PURPOSE: Generates ldt.config files for adjusting NoahMP401 ensemble restart
# files for the six different NMME models.  Also copies/renames restart files
# for the 1-member HyMAP2 restart for use as ICs for the forecast runs.  Based
# on generate_ldtconfig_files_ensrst_nrt.sh by Kristi Arsenault/SAIC.
#
# REQUIREMENTS as of 14 Oct 2021:
# * Python 3.8 or higher
#
# REVISION HISTORY:
# * 14 Oct 2021: Eric Kemp/SSAI, first version.
# * 01 Nov 2021: Eric Kemp/SSAI, addressed pylint format complaints.
# * 15 Nov 2021: K. Arsenault/SAIC, made minor fixes and added config entries.
# * 18 Jan 2022: K. Arsenault/SAIC, fixed prev to current year in naming convention
# * 27 Oct 2026: Sarith Mahanama/SSAI, replaced LDT with Python in order to process 
#                HYMAP RST files.
#------------------------------------------------------------------------------
"""


# Standard modules
import os
import datetime
import argparse
import shutil
import yaml
import numpy as np
import netCDF4 as nc
from ghis2s.bcsd.bcsd_library.nmme_module import NMMEParams
from ghis2s.shared.logging_utils import TaskLogger

# Global Constants
_LSM_NAME = ""
_ROUTING_NAME = ""
_INPUT_NUMFCSTMONS = 9
_ROUTING_SUBDOMAINS = False
_INPUT_ENS_MEMS = 12
_NMME_MODELS = []

INPUTDIR = ''
WORKDIR = ''
FORECAST_YEAR = 0
FORECAST_MONTH = 0
_LOGGER = None
_SUBTASK = ''

def _recursive_chmod(path, mode):
    """Recursively runs chmod"""
    os.chmod(path, mode)
    for dirpath, dirnames, filenames in os.walk(path):
        for dirname in dirnames:
            os.chmod(os.path.join(dirpath, dirname), mode)
        for filename in filenames:
            os.chmod(os.path.join(dirpath, filename), mode)

def _handle_dates():
    """Collect and return date information."""
    currentdate = datetime.date(FORECAST_YEAR, FORECAST_MONTH, 1)
    _LOGGER.info(f"Current year / month: {currentdate.year:04d} / {currentdate.month:02d}", subtask=_SUBTASK)
    rstdate = datetime.date(currentdate.year, currentdate.month, 1)
    _LOGGER.info(f" Rst Start Month Name: {rstdate.strftime('%b')}", subtask=_SUBTASK)
    prevdate = rstdate - datetime.timedelta(days=1)
    _LOGGER.info(f" Rst year / month / day: {prevdate.year:04d} / {prevdate.month:02d} / {prevdate.day:02d}"
                 , subtask=_SUBTASK)
    return currentdate, prevdate

def _scale_lis_restart(in_file, out_file, input_ens, target_ens):
    """scales a LIS restart file's ntiles dimension by reshaping the grouped ensembles."""

    # pylint: disable=no-member
    with nc.Dataset(in_file, 'r') as src, nc.Dataset(out_file, 'w') as dst:
        total_tiles = src.dimensions['ntiles'].size
        actual_tiles = total_tiles // input_ens
        target_tiles = actual_tiles * target_ens

        # Copy dimensions safely
        for dimname, dim in src.dimensions.items():
            if dimname == 'ntiles':
                dst.createDimension(dimname, target_tiles)
            elif dim.isunlimited():
                dst.createDimension(dimname, None)
            else:
                dst.createDimension(dimname, len(dim))

        # Copy global attributes
        dst.setncatts(src.__dict__)

        # 3. Copy, reshape, and scale variables dynamically
        for name, variable in src.variables.items():
            # extract _FillValue if it exists
            fill_val = variable._FillValue if hasattr(variable, '_FillValue') else None
            dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_val)

            # Copy all other attributes
            attrs = {k: variable.getncattr(k) for k in variable.ncattrs() if k != '_FillValue'}
            if attrs:
                dst[name].setncatts(attrs)

            if 'ntiles' in variable.dimensions:
                axis = variable.dimensions.index('ntiles')
                shp = list(variable.shape)
                data = src[name][:]

                # (actual_tiles, input_ens) preserves geography!
                new_shape = list(shp)
                new_shape[axis:axis+1] = [actual_tiles, input_ens]
                data = data.reshape(new_shape)

                if target_ens < input_ens:
                    # Drop extra ensembles from EVERY tile
                    slices = [slice(None)] * len(new_shape)
                    slices[axis+1] = slice(0, target_ens)
                    data = data[tuple(slices)]
                else:
                    # UPSCALE: Replicate member 1
                    slices = [slice(None)] * len(new_shape)
                    slices[axis+1] = slice(0, 1)
                    base_data = data[tuple(slices)]
                    data = np.repeat(base_data, target_ens, axis=axis+1)

                # Compress back into 1D ntiles
                final_shape = list(shp)
                final_shape[axis] = target_tiles
                dst[name][:] = data.reshape(final_shape)

            else:
                dst[name][:] = src[name][:]

def _driver():
    """Main driver"""
    currentdate, prevdate = _handle_dates()

    # Form input restart directory and filenames together
    lsm_rstdir = f"{INPUTDIR}/SURFACEMODEL/{prevdate.year:04d}{prevdate.month:02d}"
    hymap_rstdir = f"{INPUTDIR}/ROUTING/{prevdate.year:04d}{prevdate.month:02d}"

    for nmme_model in _NMME_MODELS:
        num_ensmems = NMMEParams(nmme_model).ens_num
        _LOGGER.info(f" Number of members: {num_ensmems:d} for model {nmme_model}",
                     subtask=_SUBTASK)

        rst_date = f"{prevdate.year:04d}{prevdate.month:02d}{prevdate.day:02d}"
        rst_monname = currentdate.strftime("%b")

        # ---------------------------------------------------------
        # 1. PROCESS LSM RESTART
        # ---------------------------------------------------------
        lsm_inrstfile = f"{lsm_rstdir}/LIS_RST_NOAHMP401_{rst_date}2345.d01.nc"
        lsm_outrstfile = (
            f"./{nmme_model}/LIS_RST_NOAHMP401_{rst_date}2345.ICS_{rst_monname}"
            f"{currentdate.year:04d}.ens{num_ensmems:d}.nc")

        input_ensmems = _INPUT_ENS_MEMS
        target_ensmems = num_ensmems

        _LOGGER.info(f" Processing NOAHMP401 restart: {input_ensmems} -> {target_ensmems} members"
                     , subtask=_SUBTASK)
        if input_ensmems == target_ensmems:
            shutil.copy(lsm_inrstfile, lsm_outrstfile)
            _LOGGER.info(f" Copied NOAHMP401 restart file to {lsm_outrstfile}", subtask=_SUBTASK)
        else:
            _scale_lis_restart(lsm_inrstfile, lsm_outrstfile, input_ensmems, target_ensmems)
            _LOGGER.info(f" Scaled NOAHMP401 restart to {target_ensmems} members.",
                         subtask=_SUBTASK)

        # ---------------------------------------------------------
        # 2. PROCESS ROUTER RESTART
        # ---------------------------------------------------------
        router_upper = _ROUTING_NAME.upper() if _ROUTING_NAME else "HYMAP2"
        hymap_inrstfile = f"{hymap_rstdir}/LIS_RST_{router_upper}_router_{rst_date}2345.d01.nc"
        hymap_outrstfile = (
            f"./{nmme_model}/LIS_RST_{router_upper}_router_{rst_date}2345.ICS_{rst_monname}"
            f"{currentdate.year:04d}.ens{num_ensmems:d}.nc")

        if _ROUTING_SUBDOMAINS:
            _LOGGER.info((f" Processing {router_upper} restart: {input_ensmems} ->"
                          f"{target_ensmems} members")
                         , subtask=_SUBTASK)
            if input_ensmems == target_ensmems:
                shutil.copy(hymap_inrstfile, hymap_outrstfile)
                _LOGGER.info(f" Copied {router_upper} restart file to {hymap_outrstfile}",
                             subtask=_SUBTASK)
            else:
                _scale_lis_restart(hymap_inrstfile, hymap_outrstfile, input_ensmems, target_ensmems)
                _LOGGER.info(f" Scaled {router_upper} restart to {target_ensmems} members.",
                             subtask=_SUBTASK)
        else:
            # Output Global router as exactly 1 member
            hymap_outrstfile_glob = (
                f"./{nmme_model}/LIS_RST_{router_upper}_router_{rst_date}2345.ICS_{rst_monname}"
                f"{currentdate.year:04d}.ens1.nc")
            shutil.copy(hymap_inrstfile, hymap_outrstfile_glob)
            _LOGGER.info(f" Copied global {router_upper} restart file to {hymap_outrstfile_glob}",
                         subtask=_SUBTASK)

    print("[INFO] Done generating restart files.")


if __name__ == "__main__":

    # input arguments
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-y', '--year', required=True, help='forecast start year')
    PARSER.add_argument('-m', '--month', required=True, help='forecast start month')
    PARSER.add_argument('-i', '--INPUTDIR', required=True, help='lisda_run/output directory')
    PARSER.add_argument('-w', '--WORKDIR', required=True, help='working directory')
    PARSER.add_argument('-s', '--CONFIGFILE', required=True, help='s2s.config file')
    ARGS = PARSER.parse_args()

    FORECAST_YEAR = int(ARGS.year)
    FORECAST_MONTH = int(ARGS.month)
    INPUTDIR = ARGS.INPUTDIR
    WORKDIR = ARGS.WORKDIR
    CONFIGFILE = ARGS.CONFIGFILE

    # Read s2s.config
    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        CONFIG = yaml.safe_load(file)

    # update global params
    _INPUT_NUMFCSTMONS = CONFIG['EXP']['lead_months']
    _LSM_NAME = CONFIG['EXP']['lsm']
    _ROUTING_NAME = CONFIG['EXP']['routing_model'].get('model')
    _ROUTING_SUBDOMAINS = CONFIG['EXP']['routing_model'].get('subdomains', False)
    if _ROUTING_SUBDOMAINS:
        _SUBTASK = os.path.basename(os.getcwd())
    _NMME_MODELS = CONFIG['EXP']['NMME_models']
    _INPUT_ENS_MEMS = CONFIG['BCSD'].get('nof_raw_ens', 12)
    _TASK_NAME = os.environ.get('SCRIPT_NAME')
    _LOGGER = TaskLogger(_TASK_NAME,
                         os.getcwd(),
                         f'ldt_ics/generate_ldtconfig_files_ensrst_nrt.py creating ICs files')

    _driver()
