CC=mpicc

%.o:%.c
	$(CC) $(CFLAGS) $< -o $@
hwk5make: main.c
	#module load gnu openmpi # module is a function type module
	# be sure to call "make debug" in order to see homework 4 print debug messages
	$(CC) -o hwk5 -lm main.c 
debug:
	#module load gnu openmpi # module is a function type module
	#$(call module,load gnu openmpi)
	$(CC) -o hwk5 -lm main.c -DDEBUG
clean:
	rm hwk5

#$(call module,load gnu openmpi)
define module
    $(shell eval $(${LMOD_CMD} bash load gnu openmpi); [ $? = 0 ] && eval $(${LMOD_SETTARG_CMD:-:} -s sh))
endef 
