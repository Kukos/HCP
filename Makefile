export

CC := gcc
RM := rm -rf

CCWARNINGS := 	-Wall -Wextra -pedantic -Wcast-align \
		-Winit-self -Wlogical-op -Wmissing-include-dirs \
		-Wredundant-decls -Wshadow -Wstrict-overflow=5 -Wundef  \
		-Wwrite-strings -Wpointer-arith -Wmissing-declarations \
		-Wuninitialized -Wold-style-definition -Wstrict-prototypes \
		-Wmissing-prototypes -Wswitch-default -Wbad-function-cast \
		-Wnested-externs -Wconversion -Wunreachable-code

CFLAGS := -std=gnu99 $(CCWARNINGS) -O3

PROJECT_DIR := $(shell pwd)
EDIR := $(PROJECT_DIR)/external
SUBDIR := $(PROJECT_DIR)/submodules
LDIR := $(EDIR)/libs
EIDIR := $(EDIR)/include
ESDIR := $(EDIR)/libs

ifeq ("$(origin V)", "command line")
  VERBOSE = $(V)
endif

ifndef VERBOSE
  VERBOSE = 0
endif

ifeq ($(VERBOSE),1)
  Q =
else
  Q = @
endif

define print_info
	$(if $(Q), @echo "$(1)")
endef

define print_make
	$(if $(Q), @echo "[MAKE]     $(1)")
endef

define print_cc
	$(if $(Q), @echo "[CC]      $$(1)")
endef

define print_bin
	$(if $(Q), @echo "[BIN]     $$(1)")
endef

all: l1

libs:
	$(Q)if [ ! -d $(LDIR) ]; then \
	cd $(SUBDIR)/MyLibs/scripts && \
	./install_libs.sh $(EDIR) && \
	cd $(PROJECT_DIR) ;fi

l1: libs
	$(call print_make,list1)
	$(Q)$(MAKE) -f $(PROJECT_DIR)/list1/Makefile --no-print-directory

clean:
	$(call print_info,Cleaning)
	$(Q)$(MAKE) -f $(PROJECT_DIR)/list1/Makefile clean --no-print-directory

clean_libs:
	$(Q)$(RM) $(EDIR)/*
	$(Q)cd $(SUBDIR)/MyLibs && $(MAKE) clean --no-print-directory
