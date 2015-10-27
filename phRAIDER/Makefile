# Makefile is part of phRAIDER.
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

# Created by Carly Schaeffer, Nathan Figueroa, and John Karro

CC = g++ -std=c++11 -Wno-deprecated-declarations -Wno-header-guard -Wno-deprecated-register

OPTS = -O3
CXXFLAGS+=-Wall -W -DNDEBUG $(OPTS)

CPPINCL+=-I../seqan/include
CPPINCL+=-I../seqan/extras/include

###########################################################
# RAIDER
###########################################################
phRAIDER: Family.h LmerVector.h SeedChain.h AppOptions.h main.cpp
	@echo compiling phRAIDER
	$(CC) $(CPPINCL) $(CXXFLAGS) main.cpp -o phRAIDER


clean:
	rm *.o ../phRAIDER tests

