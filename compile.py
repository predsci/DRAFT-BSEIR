#!/usr/bin/python
import os
import string
import sys

os.system("/bin/rm -rf DRAFT_*tar.gz draft/src/*.o draft/src/*.so")
os.system("R CMD build draft")
os.system("R CMD INSTALL draft")

## For a local installation - comment the above line and uncomment the lines starting with 'mkdir', 'export' and 'os.system'
## Create a directory in the home directory, called for example R_libs:
## mkdir /home/your_username/R_libs
## set a variable to point R at that directory:
## export R_LIBS="/home/your_username/R_libs"
## and then install the package ih this directory:
#os.system("R CMD INSTALL $HOME/R_LIBS bSEIR_1.0.tar.gz")

