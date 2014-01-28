CC	:=	gcc
CFLAGS	:=	-O2
LDFLAGS	:=	-lm
SOURCEDIR	:= src
HEADERS	:=	$(SOURCEDIR)/def_param.h
SOURCES	:=	$(SOURCEDIR)/main_reg.c $(SOURCEDIR)/read_protein.c $(SOURCEDIR)/initialize1.c $(SOURCEDIR)/rforce.c $(SOURCEDIR)/force_2poly.c $(SOURCEDIR)/iteration.c $(SOURCEDIR)/update.c $(SOURCEDIR)/write_protein.c $(SOURCEDIR)/identify.c $(SOURCEDIR)/tqli_dub.c $(SOURCEDIR)/tred2_dub.c
BIN	:=	run_myo

$(BIN): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
	
clean:
	rm -r Struct_data Data $(BIN) output
