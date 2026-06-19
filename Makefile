MINCO_BIN := bin/minco
STAGE0_BIN := bin/minco_stage0
STAGE1_BIN := bin/minco_stage1
STAGE2_BIN := bin/minco_stage2
STAGE3_BIN := bin/minco_stage3
STAGE3_NATIVE_BIN := bin/minco_stage3_native
STAGE3_SPARSE_NATIVE_BIN := bin/minco_stage3_sparse_native
MINCO_CORE_DIR := minco_core
MINCO_CORE_CFLAGS := -std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp

.PHONY: all clean test stage0 stage1 stage2 stage3 stage3_native stage3_sparse_native minco stages install_env

all: minco

stage0: | bin
	$(MAKE) -C $(MINCO_CORE_DIR) PRONAME=minco_stage0 BINDIR=$(abspath bin)

stage1: | bin
	$(MAKE) -C $(MINCO_CORE_DIR) PRONAME=minco_stage1 BINDIR=$(abspath bin) OBJDIR=$(abspath $(MINCO_CORE_DIR)/obj_stage1) CFLAGS="$(MINCO_CORE_CFLAGS) -DMINCO_HASH_BOTTOMK=1 -DMINCO_KEEP_SOURCE_FILTER=1"

stage2: | bin
	$(MAKE) -C $(MINCO_CORE_DIR) PRONAME=minco_stage2 BINDIR=$(abspath bin) OBJDIR=$(abspath $(MINCO_CORE_DIR)/obj_stage2) CFLAGS="$(MINCO_CORE_CFLAGS) -DMINCO_HASH_BOTTOMK=1 -DMINCO_KEEP_SOURCE_FILTER=0"

stage3: | bin
	$(MAKE) -C $(MINCO_CORE_DIR) PRONAME=minco_stage3 BINDIR=$(abspath bin) OBJDIR=$(abspath $(MINCO_CORE_DIR)/obj_stage3) CFLAGS="$(MINCO_CORE_CFLAGS) -DMINCO_HASH_BOTTOMK=1 -DMINCO_KEEP_SOURCE_FILTER=0 -DMINCO_STREAM_BOTTOMK=1"

stage3_native: | bin
	$(MAKE) -C $(MINCO_CORE_DIR) PRONAME=minco_stage3_native BINDIR=$(abspath bin) OBJDIR=$(abspath $(MINCO_CORE_DIR)/obj_stage3_native) CFLAGS="$(MINCO_CORE_CFLAGS) -march=native -mavx2 -mbmi2 -DMINCO_HASH_BOTTOMK=1 -DMINCO_KEEP_SOURCE_FILTER=0 -DMINCO_STREAM_BOTTOMK=1"

stage3_sparse_native: | bin
	$(MAKE) -C $(MINCO_CORE_DIR) PRONAME=minco_stage3_sparse_native BINDIR=$(abspath bin) OBJDIR=$(abspath $(MINCO_CORE_DIR)/obj_stage3_sparse_native) CFLAGS="$(MINCO_CORE_CFLAGS) -march=native -mavx2 -mbmi2 -DMINCO_HASH_BOTTOMK=1 -DMINCO_KEEP_SOURCE_FILTER=0 -DMINCO_STREAM_BOTTOMK=1 -DMINCO_HASH_SPARSE_CTX=1"

minco: stage3_native
	cp $(STAGE3_NATIVE_BIN) $(MINCO_BIN)

stages: stage0 stage1 stage2 stage3 stage3_native stage3_sparse_native minco

install_env:
	@printf '%s\n' 'export PATH="$(abspath bin):$$PATH"'
	@printf '%s\n' 'Add the line above to your shell profile if you want to run minco directly.'

bin:
	mkdir -p bin

test: minco
	bash tests/smoke.sh
	bash tests/full_cli.sh

clean:
	rm -rf bin $(MINCO_CORE_DIR)/obj*
