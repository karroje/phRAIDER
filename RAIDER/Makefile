# Makefile is part of RAIDER.
#
# phRAIDER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# phRAIDER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with phRAIDER.  If not, see <http://www.gnu.org/licenses/>.

# Created by Nathan Figueroa, Carley Schaeffer, and John Karro
CC = g++ -std=c++11

CXXFLAGS+=-Wall -W -O3 -DNDEBUG

CPPINCL+=-I../seqan/include
CPPINCL+=-I../seqan/extras/include

###########################################################
# RAIDER
###########################################################
all: main

main: Family.h LmerVector.h SeedChain.h main.cpp
	@echo compiling RAIDER; 
	$(CC) $(CPPINCL) $(CXXFLAGS) $^ -o raider

