#    Copyright (C) 2014 Benjamin E. Decato
#
#    Authors: Benjamin E. Decato, University of Southern California
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifeq "$(shell uname)" "Darwin"
	CXXFLAGS += -arch x86_64
endif

ifdef DEBUG
	CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT 
	CXXFLAGS += $(OPTFLAGS)
endif

all: galignbd lalignbd banded_galignbd malignbd bslalignbd

galignbd: galignbd.cpp
	$(CXX) $(CXXFLAGS) -o $@ galignbd.cpp

lalignbd: lalignbd.cpp 
	$(CXX) $(CXXFLAGS) -o $@ lalignbd.cpp

banded_galignbd: banded_galignbd.cpp
	$(CXX) $(CXXFLAGS) -o $@ banded_galignbd.cpp

malignbd: malignbd.cpp
	$(CXX) $(CXXFLAGS) -o $@ malignbd.cpp

bslalignbd:  bslalignbd.cpp
	$(CXX) $(CXXFLAGS) -o $@ bslalignbd.cpp

clean: 
	@-rm -f *.o *~
.PHONY: clean

