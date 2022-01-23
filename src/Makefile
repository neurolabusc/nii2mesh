# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS=-lz -lm
#c flags
CFLAGS=-O3 -DHAVE_ZLIB -DHAVE_JSON
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
	CFLAGS+=-fopenmp -flto
	LFLAGS+= -static-libgcc -static-libstdc++ -L.
	NOTES=OpenMP will demand all CPUs!
endif


all:
	$(info $(NOTES))
	$(CNAME) $(CFLAGS) nii2mesh.c cJSON.c isolevel.c meshify.c quadric.c base64.c bwlabel.c radixsort.c -o nii2mesh $(LFLAGS)
	
