.KEEP_STATE:

# source code
# EXES: individual binary executables
# OBJS: non-executable object files
EXE = MetaPhase
OBJS = TextFileParsers.o LinkMatrix.o HierarchicalClustering.o TrueMapping.o ClusteringResult.o MetaAssembly.o MetaPhase.o

# compiler commands

CC      = g++
RM      = /bin/rm -rf
BACKUPS = *~ \\\#*\\\#


# compiler flags
CFLAGS = -Wall -g -O3
CFLAGS += -ansi
CFLAGS += -pedantic
#CFLAGS += -pg

# Directories for includes & linking flags
INCLUDE_DIR=$(HOME)/include
EXTERN_DIR=/net/shendure/vol10/jnburton/extern
BOOST_DIR=$(EXTERN_DIR)/boost_1_47_0
SAMTOOLS_DIR=$(EXTERN_DIR)/samtools-0.1.16

# Includes
INCLUDES=-I$(INCLUDE_DIR) -I$(SAMTOOLS_DIR)

# Linking flags
# If linking works but you get the runtime error "cannot open shared object file", make sure the boost lib directory is in your $LD_LIBRARY_PATH.
INC_LIBS=-L$(INCLUDE_DIR) -lJtime -lJgtools -lJmarkov
BOOST_LIBS=-L$(BOOST_DIR)/stage/lib -lboost_filesystem #-lboost_system #-lboost_regex
SAMTOOLS_LIBS=-L$(SAMTOOLS_DIR) -lbam
LFLAGS = -lz $(INC_LIBS) $(BOOST_LIBS) $(SAMTOOLS_LIBS)


# dependencies

.cc.o:  .cc
	$(CC) -c $< $(CFLAGS) $(INCLUDES)

all:    $(EXE)

libs:
	cd $(INCLUDE_DIR); make all; cd -

MetaPhase:  $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o MetaPhase $(LFLAGS)

clean:
	$(RM) $(OBJS) core .make.state gmon.out
	#cd $(INCLUDE_DIR); make clean; cd -

clobber: clean
	$(RM) $(BACKUPS) $(EXE)

                                                                               
