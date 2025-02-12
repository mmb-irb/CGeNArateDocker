define obj_files
$(addprefix $1/, $(addsuffix .o, $(notdir $(basename $(SRC)))))
endef

define target_rule
$1: $2
	$(CC) -o $$@ $2 $(CFLAGS) $3 $(LDFLAGS)
endef

define obj_rule
$1/%.o: $(SRC_PATH)/%.c
	$(CC) -c $(CFLAGS) $2 -o $$@ $$<
endef

# Common flags:
CFLAGS += -Wall -Wno-extra-tokens -I$(SRC_PATH) -D_XOPEN_SOURCE=500
LDFLAGS = -lm

ifneq ($(filter glogin%,$(HOSTNAME)),) # marenostrum login node
	# For OpenMP-V:
	OMP_FLAG ?= -fopenmp=libompv
else
	OMP_FLAG ?= -fopenmp
endif

# Compiler
STD = c11
CC = clang # CC ?= clang
CC += -std=$(STD)

SRC_PATH := src
SRC := $(foreach x, $(SRC_PATH), $(wildcard $(addprefix $(x)/*,.c)))

OPT = -O3
DEBUG_FLAGS = -g -DDEBUG

REL_FLAGS = $(OMP_FLAG) $(OPT)
DEB_FLAGS = $(OMP_FLAG) $(DEBUG_FLAGS)
SEQ_FLAGS = $(OPT)
DBS_FLAGS = $(DEBUG_FLAGS)

REL_TARGET = CGeNArate.exe
DEB_TARGET = CGeNArate.debug.exe
SEQ_TARGET = CGeNArate.seq.exe
DBS_TARGET = CGeNArate.debug.seq.exe

REL_OBJ_PATH := build_rel
DEB_OBJ_PATH := build_deb
SEQ_OBJ_PATH := build_seq
DBS_OBJ_PATH := build_dbs

REL_OBJ := $(call obj_files, $(REL_OBJ_PATH))
DEB_OBJ := $(call obj_files, $(DEB_OBJ_PATH))
SEQ_OBJ := $(call obj_files, $(SEQ_OBJ_PATH))
DBS_OBJ := $(call obj_files, $(DBS_OBJ_PATH))

all: prep release debug seq seq-debug

$(eval $(call target_rule, $(REL_TARGET), $(REL_OBJ), $(REL_FLAGS)))
$(eval $(call target_rule, $(DEB_TARGET), $(DEB_OBJ), $(DEB_FLAGS)))
$(eval $(call target_rule, $(SEQ_TARGET), $(SEQ_OBJ), $(SEQ_FLAGS)))
$(eval $(call target_rule, $(DBS_TARGET), $(DBS_OBJ), $(DBS_FLAGS)))

$(eval $(call obj_rule, $(REL_OBJ_PATH), $(REL_FLAGS)))
$(eval $(call obj_rule, $(DEB_OBJ_PATH), $(DEB_FLAGS)))
$(eval $(call obj_rule, $(SEQ_OBJ_PATH), $(SEQ_FLAGS)))
$(eval $(call obj_rule, $(DBS_OBJ_PATH), $(DBS_FLAGS)))

release: $(REL_TARGET)
debug: $(DEB_TARGET)
seq: $(SEQ_TARGET)
seq-debug: $(DBS_TARGET)

check: $(REL_TARGET)
	OMP_NUM_THREADS=8 ./$(REL_TARGET) simulations/1j5n.toml
	./check.py

prep:
	@mkdir -p $(REL_OBJ_PATH) $(DEB_OBJ_PATH) $(SEQ_OBJ_PATH) $(DBS_OBJ_PATH)

.PHONY : clean release debug seq seq-debug prep check
clean:
	-rm -f $(REL_TARGET) $(REL_OBJ) $(DEB_TARGET) $(DEB_OBJ) $(SEQ_TARGET) $(SEQ_OBJ) $(DBS_OBJ)
