
SRC = $(wildcard $(FEWARPDIR)FEWarp/*.cpp)
OBJ = $(patsubst $(FEWARPDIR)FEWarp/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(FEWARPDIR)FEWarp/%.cpp, %.d, $(SRC))


SO = libfewarp_$(PLAT).$(SFX)
LIB = $(FEWARPDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore_$(PLAT).a

FEBIOMECH = $(FEBLIB)/libfebiomech_$(PLAT).a

NUMCORE = $(FEBLIB)/libnumcore_$(PLAT).a

FEBIOLIBS = $(NUMCORE) $(FECORE) $(FEBIOMECH)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS)
endif

%.o: $(FEWARPDIR)FEWarp/%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
