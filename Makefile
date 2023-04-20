MCC=mpicc  
CXX=mpic++
#ARCH=ARCH_MPE
ARCH=ARCH_SW
MCFLAGS=-O3  -w -D$(ARCH) -DNO_BLAS -mieee
CPPFLAGS=-O3 -w -D$(ARCH) -DNO_BLAS -mieee
CC=swgcc 
CCFLAGS=-mslave -O3 -D$(ARCH) -msimd -w -mieee -faddress_align=32
FC=mpif90 -mieee 
FFLAGS=-O3 -mhost -w
LDFLAGS=-Wl,--wrap=GPTLstart
DEP=cpp
DEPFLAGS=-MM

SQLITEHOME := $INSTALL_PATH_SQL
PROJHOME := $INSTALL_PATH_PROJ
NETCDFHOME := $INSTALL_PATH_NETCDF

GPTL := $INSTALL_PATH_GPTL


LIBS := -Wl,-Bstatic -L$(PROJHOME)/lib -lproj
INCS := -I$(PROJHOME)/include  

LIBS += -Wl,-Bstatic -L$(SQLITEHOME)/lib -lsqlite3
INCS += -I$(SQLITEHOME)/include  

LIBS += -Wl,-Bstatic -L$(NETCDFHOME)/lib -lnetcdf
INCS += -I$(NETCDFHOME)/include

LIBS += -L$(GPTL)/lib -lgptl  
INCS += -I$(GPTL)/include

LIBS += -Wl,-Bdynamic -L/usr/lib64 #-lpthread

BIN_PATH=./bin
BIN_NAME=cgfd3d
BUILD_PATH=./build
SRC_PATH=.
CPE_PATH=$(SRC_PATH)/src.cpe
SRC_C =$(wildcard $(SRC_PATH)/*.c)
OBJS = $(SRC_C:$(SRC_PATH)/%.c=$(BUILD_PATH)/%.o)
SRC_CPP = $(wildcard $(SRC_PATH)/*.cpp)
OBJS += $(SRC_CPP:$(SRC_PATH)/%.cpp=$(BUILD_PATH)/%.o)
SRC_CPE = $(wildcard $(CPE_PATH)/*.c)
OBJS += $(SRC_CPE:$(CPE_PATH)/%.c=$(BUILD_PATH)/%.o)

all: $(OBJS) 
	$(CXX)  $(OBJS) $(LDFLAGS) $(LIBS)  -o $(BIN_PATH)/$(BIN_NAME)
$(BUILD_PATH)/%.o : $(SRC_PATH)/%.cpp
	$(CXX) $(CPPFLAGS) $(INCS) -c $< -o $@
$(BUILD_PATH)/%.o : $(SRC_PATH)/%.c
	$(MCC) $(MCFLAGS) $(INCS) -c $< -o $@
$(BUILD_PATH)/%.o : $(CPE_PATH)/%.c
	$(CC) $(CCFLAGS) $(INCS) -c $< -o $@

.PHONY: print
print:
	@echo $(SRC_C)
	@echo $(SRC_CPP)
	@echo $(OBJS)

.PHONY: clean
clean:
	-rm -rf $(BIN_PATH)/$(BIN_NAME)
	-rm -rf $(BUILD_PATH)/*.o

.PHONY: install
install:
	cp ../bin/$(BIN_NAME) $(BIN_PATH)/$(BIN_NAME)
