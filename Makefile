ifdef OP2_INSTALL_PATH
  OP2_INC = -I$(OP2_INSTALL_PATH)/include
  OP2_LIB = -L$(OP2_INSTALL_PATH)/lib
endif

default: build
	echo "Start Build"

CXX = g++
EXE = miniAero.host
KOKKOS_ARCH = "SNB"
CXXFLAGS = -O3 -g -DATOMICS_FLUX
LINK = ${CXX}
LINKFLAGS = -lop2_seq

DEPFLAGS = -M

SRC = $(wildcard *.C)
OBJ = $(SRC:.C=.o)
LIB =

CXXFLAGS += -Isrc

build: ${EXE}

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) $(OP2_LIB) $(OP2_INC) -lop2_seq -o $(EXE)

clean: 
	rm -f *.o *.cuda *.host

test: $(EXE)
	cd tests; \
	./run_tests.sh ../$(EXE)

# Compilation rules

%.o:%.C $(KOKKOS_CPP_DEPENDS)
	$(CXX) -I./ -std=c++14 -mavx $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) $(OP2_LIB) $(OP2_INC) -c $<