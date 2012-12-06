ALG ?= naive

SOMMELIER := sommelier
SOMMELIER-SESC := sommelier.sesc

SOMMELIER-SRCS := sommelier.c smat.c $(ALG).c 
SOMMELIER-DEPS := sommelier.h smat.h job_api.h strassen.c strassenParallel.c stassenParallel.h
SOMMELIER-OBJS := $(SOMMELIER-SRCS:%.c=%.o)
SOMMELIER-SESC-OBJS := $(SOMMELIER-SRCS:%.c=%.sesc.o)

MIPS-CC ?= ../utils/bin/mipseb-linux-gcc
MIPS-CFLAGS := -g -mips2 -mabi=32 -Wa,-non_shared -mno-abicalls
MIPS-CFLAGS += -DSESC -I../utils/lib/gcc/mipseb-linux/3.4.4/include
MIPS-CFLAGS += -I../utils/mipseb-linux/sys-include
MIPS-CFLAGS += -I../utils/mipseb-linux/include
MIPS-CFLAGS += -L../utils/lib/gcc/mipseb-linux/3.4.4 -nostdinc
MIPS-CLINK := -static -Wl,--script=../utils/mipseb-linux/lib/ldscripts/mint.x
SOMMELIER-SESC-OBJS += sesc_thread.sesc.o sesc_events.sesc.o

CC ?= gcc
CFLAGS := -Wall -g
CLINK := -lpthread

SCRIPTS := ../scripts

EXAMPLE := job_api_example
EXAMPLE-SESC := job_api_example.sesc

EXAMPLE-SRCS := job_api_example.c
EXAMPLE-OBJS := $(EXAMPLE-SRCS:%.c=%.o)
EXAMPLE-SESC-OBJS := $(EXAMPLE-SRCS:%.c=%.sesc.o)
EXAMPLE-SESC-OBJS += sesc_thread.sesc.o sesc_events.sesc.o

all: $(SOMMELIER) $(SOMMELIER-SESC) $(CLINK)

$(SOMMELIER): $(SOMMELIER-OBJS)
	$(CC) -o $@ $^ $(CLINK)

%.o: %.c $(SOMMELIER_DEPS)
	$(CC) $(CFLAGS) -o $@ -c $<

$(SOMMELIER-SESC): $(SOMMELIER-SESC-OBJS)
	$(MIPS-CC) $(MIPS-CFLAGS) -o $@ $^ $(MIPS-CLINK)

%.sesc.o: %.c $(SOMMELIER_DEPS)
	$(MIPS-CC) $(MIPS-CFLAGS) -o $@ -c $<

checkpatch:
	$(SCRIPTS)/checkpatch.pl --no-tree --no-signoff --quiet -f $(SOMMELIER-SRCS) $(SOMMELIER-DEPS)
	@echo "Style OK for ALG=$(ALG)"

sparse:
	sparse $(SOMMELIER-SRCS) $(SOMMELIER-DEPS)
	@echo "Sparse-clean pass for ALG=$(ALG)"

test: $(SOMMELIER)
	$(MAKE) -C $@

example: $(EXAMPLE) $(EXAMPLE-SESC)

$(EXAMPLE): $(EXAMPLE-OBJS)
	$(CC) -o $@ $^ $(CLINK)

$(EXAMPLE-SESC): $(EXAMPLE-SESC-OBJS)
	$(MIPS-CC) $(MIPS-CFLAGS) -o $@ $^ $(MIPS-CLINK)

clean:
	$(MAKE) -C test clean
	rm -rf *.o $(SOMMELIER) $(SOMMELIER-SESC)
	rm -rf $(EXAMPLE) $(EXAMPLE-SESC)

.PHONY: all checkpatch test clean example
