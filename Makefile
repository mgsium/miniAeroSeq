default: build
	echo "Start Build"

CXX = g++
EXE = miniAero.host
KOKKOS_ARCH = "SNB"
CXXFLAGS = -O3 -g -DATOMICS_FLUX
LINK = ${CXX}
LINKFLAGS =  

DEPFLAGS = -M

SRC = $(wildcard *.C)
OBJ = $(SRC:.C=.o)
LIB =

CXXFLAGS += -Isrc

build: ${EXE}

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) $(OP2_LIB) $(OP2_INC) -I./ ./*.c -o $(EXE)

clean: rm -f *.o *.cuda *.host

test: $(EXE)
	cd tests; \
	./run_tests.sh ../$(EXE)

# Compilation rules

%.o:%.C $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) $(OP2_LIB) $(OP2_INC) -c $<