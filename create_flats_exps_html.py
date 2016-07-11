#!/usr/bin/env python

"""
Creates an html page with all the png images in current directory,
ordered by night.
"""

import os
import sys
import subprocess
from coreutils.desdbi import *
import argparse

parser = argparse.ArgumentParser(description="Check if flats is appropiate for using it")
parser.add_argument("-d", "--dir", required=False, dest="directory", help="Directory where the png files are")

args = parser.parse_args()

path = args.directory()


def main():
    html_path = '/home/ricardoc/public_html/Flats_pngs'
    
    png_path = '/home/ricardoc/DES-SNe/Manifest_database/DomeFlatsExposures'
    
    filenames = os.listdir(png_path)
    
    for png in filenames:
        print png
    


    return 0

    
if  __name__ == '__main__':
    status = main()
    sys.exit(status)

    