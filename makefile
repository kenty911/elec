# Define compiler and flags
CC = gcc
CFLAGS = -cpp -fPIC -shared -lm -O3

# Targets and their corresponding object files
TARGETS = calc_adv.so calc_charge.so calc_gauss.so flash.so lookup.so

# Phony target to avoid conflicts with files named clean or log
.PHONY: all clean log

# Default target
all: $(TARGETS)

# Pattern rule to compile .c files to shared objects
%.so: %.c
	$(CC) $(CFLAGS) $< -o $@

# Clean up generated files
clean:
	rm -f $(TARGETS)

# Log the time of compilation
log:
	@mv nohup.out $$(date +%s).log
