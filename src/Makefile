# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS= -lm
#c flags
ZFLAGS= -lz -DHAVE_ZLIB
#c flags
CFLAGS= -DNII2MESH
NOTES=
ifeq ($(CXX),g++)
	#not valid for c++
	# CFLAGS += -std=gnu99
	NOTES=gcc slow marching cubes (to compile with clang/llvm 'make CXX=clang')
endif
ifeq ($(CXX),clang)
	CNAME=clang
endif
#run "OMP=1 make -j" for both OpenMP
ifeq "$(OMP)" "1"
	CFLAGS+= -fopenmp -flto
	LFLAGS+= -static-libgcc -static-libstdc++ -L.
	NOTES=OpenMP will demand all CPUs!
endif
#run "JSON=1 make -j" 
ifeq "$(JSON)" "1"
	CFLAGS+=  -DHAVE_JSON cJSON.c
endif

# C code
CFILES= isolevel.c meshify.c quadric.c bwlabel.c radixsort.c nii2mesh.c 

# output name
OFILES= -o nii2mesh

# Web Assembly only
WASM := -s DEMANGLE_SUPPORT=1 -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap", "FS_createDataFile", "FS_readFile", "FS_unlink"]' -s ALLOW_MEMORY_GROWTH=1 -s WASM=1 -s EXPORTED_FUNCTIONS='["_simplify"]'

# Modern marching cubes (vs classic)
# MFILES= MarchingCubes.c
MFILES= -DUSE_CLASSIC_CUBES oldcubes.c

# GIfTI, obj, vtk, stl, etc support
GFILES= -DHAVE_FORMATS base64.c

all: release

release:
	$(info $(NOTES))
	$(CNAME) -O3 $(CFLAGS) $(CFILES) $(GFILES) $(MFILES) $(LFLAGS) $(ZFLAGS) $(OFILES)

obj2:
	$(CNAME) -O3 $(CFLAGS) $(CFILES) $(GFILES) $(MFILES) $(LFLAGS) $(ZFLAGS) -o obj2mesh

classic:
	$(CNAME) -O3 $(CFLAGS) $(CFILES) $(GFILES) -DUSE_CLASSIC_CUBES oldcubes.c $(LFLAGS) $(ZFLAGS) $(OFILES)

# tiny: old cubes without zlib/GIFTI: create a command line executable that emulates WASM
tiny:
	$(CNAME) -O3 $(CFLAGS) $(CFILES) -DUSE_CLASSIC_CUBES oldcubes.c $(LFLAGS) $(OFILES)

wasm:
	em++ -O3 $(CFLAGS) $(CFILES) -DUSE_CLASSIC_CUBES oldcubes.c $(LFLAGS) ${WASM} -o nii2mesh.js

sanitize:
	g++ -O1 -g -fsanitize=address -fno-omit-frame-pointer $(CFLAGS) $(CFILES) $(GFILES) $(MFILES) $(LFLAGS) $(ZFLAGS) $(OFILES)


	
