CC     ?= cc
# Optimisation level. For a debug build: make clean && make OPT="-O0 -g"
OPT    ?= -O2
CFLAGS += -std=c99 -Wall $(OPT)
LDLIBS += -lz -lm

BUILD = build
SRC   = src/dude.c src/iondude.c src/seq.c
HDR   = src/dude.h src/seq.h lib/klib/kseq.h lib/klib/khash.h

GIT_BRANCH      := $(shell git symbolic-ref --short HEAD 2>/dev/null || echo unknown)
GIT_COMMIT_HASH := $(shell git describe --always --dirty 2>/dev/null || echo "none yet")

.PHONY: all test clean FORCE

all: $(BUILD)/iondude

$(BUILD)/iondude: $(SRC) $(HDR) $(BUILD)/version.h
	$(CC) $(CFLAGS) -I$(BUILD) $(SRC) -o $@ $(LDFLAGS) $(LDLIBS)

# rewritten only when the branch or the commit changes, so make rebuilds after a new commit
$(BUILD)/version.h: src/version.h.in FORCE
	@mkdir -p $(BUILD)
	@sed -e 's|@GIT_BRANCH@|$(GIT_BRANCH)|' -e 's|@GIT_COMMIT_HASH@|$(GIT_COMMIT_HASH)|' $< > $@.tmp
	@if cmp -s $@.tmp $@; then rm $@.tmp; else mv $@.tmp $@; fi

test: $(BUILD)/iondude
	$(BUILD)/iondude -t

clean:
	rm -rf $(BUILD)
