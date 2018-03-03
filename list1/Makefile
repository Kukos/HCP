THIS_DIR := $(PROJECT_DIR)/list1

IDIR := $(THIS_DIR)/include
SDIR := $(THIS_DIR)/src
LDIR :=

SRCS := $(wildcard $(SDIR)/*.c)
SRCS += $(ESDIR)/log.c
OBJS := $(SRCS:%.c=%.o)
DEPS := $(wildcard $(IDIR)/*.h)
DEPS += $(wildcard $(EIDIR)/*.h)

LIBS := -lgmp

EXEC := $(THIS_DIR)/pollard.out

%.o: %.c
	$(call print_cc, $<)
	$(Q)$(CC) $(CFLAGS) -I$(IDIR) -I$(EIDIR) -c $< -o $@

$(EXEC): $(OBJS)
	$(call print_bin, $@)
	$(Q)$(CC) $(CFLAGS) -L$(LDIR) -I$(IDIR) -I$(EIDIR) $(OBJS) $(LIBS) -o $@

clean:
	$(call print_info,Cleaning)
	$(Q)$(RM) $(OBJS)
	$(Q)$(RM) $(T_OBJS)
	$(Q)$(RM) $(EXEC)
	$(Q)$(RM) $(T_EXEC)
