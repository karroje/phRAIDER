# Nathan Figueroa

CC = g++ -std=c++11

CXXFLAGS+=-Wall -W -O3 -DNDEBUG

CPPINCL+=-I../../seqan/include
CPPINCL+=-I../../seqan/extras/include

###########################################################
# RAIDER
###########################################################
all: main

main: Family.h LmerVector.h SeedChain.h main.cpp
	@echo compiling RAIDER; 
	$(CC) $(CPPINCL) $(CXXFLAGS) $^ -o raider
#	$(CC) $(CPPINCL) $(CXXFLAGS) $^ -o ../../raider


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
