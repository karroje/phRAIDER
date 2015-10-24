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
all: ../../phRAIDER ../../pre-phRAIDER

../../phRAIDER: phRAIDER
	cp phRAIDER ../../

phRAIDER: Family.h LmerVector.h SeedChain.h AppOptions.h main.cpp
	@echo compiling phRAIDER
	$(CC) $(CPPINCL) $(CXXFLAGS) main.cpp -o phRAIDER


../../pre-phRAIDER: pre-phRAIDER
	cp pre-phRAIDER ../../

pre-phRAIDER: Family.h LmerVector.h SeedChain.h AppOptions.h main.cpp
	@echo compiling pre-phRAIDER
	$(CC) $(CPPINCL) $(CXXFLAGS) main.cpp -o pre-phRAIDER -DPRE



############################################################
# GOOGLE TEST SETUP
############################################################

GTEST_DIR = ../../../gtest-1.6.0
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

TESTFLAGS+=-Wall -Wextra
TESTFLAGS+=-I$(GTEST_DIR)/include

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(TESTFLAGS) -I$(GTEST_DIR) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(TESTFLAGS) -I$(GTEST_DIR) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

#############################################################
# UNIT TESTS
#############################################################

tests.o: Family.h LmerVector.h SeedChain.h tests.cpp $(GTEST_HEADERS)
	@echo compiling tests.o; $(CC) $(TESTFLAGS) $(CPPINCL) -c tests.cpp

tests: Family.h LmerVector.h SeedChain.h tests.o gtest_main.a
	@echo compiling tests; $(CC) $(TESTFLAGS) $(CPPINCL) $^ -lpthread -o tests

clean:
	rm *.o ../raider tests


check-syntax:
	gcc -o nul -S ${CHK_SOURCES}
