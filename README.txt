
*******************************************************************************
Copyright and License Information
*******************************************************************************
Copyright (C) 2014 Benjamin E. Decato
  
Author: Benjamin E. Decato, University of Southern California
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


*******************************************************************************
Building and Installing 
*******************************************************************************
This software package has been designed to operate in a UNIX-like environment.
It has been tested on Ubuntu 12.04 LTS, and will soon be tested on Mac OS X
Lion. g++ version 4.8.1 or higher is required: if you do not have g++ installed,
it can be installed by getting xcode for Mac OS X, or by googling win64/win32
(depending on your architecture) g++ and looking for a download.

Step 1
------
  To build the binaries, type the following, where '>' is your prompt and the
  CWD is the root of the distribution  
  
  > make all 
  
Step 2
------
  To wipe the directory clean for a fresh build, type the following, where '>'
  is your prompt and the CWD is the root of the distribution
  
  > make clean
 
  This will remove all binaries and intermediate object files. 
  
*******************************************************************************
Usage
*******************************************************************************
Read alignmentbd.pdf in the docs directory.

************************
Contacts and bug reports
*******************************************************************************
Benjamin E. Decato
decato@usc.edu
