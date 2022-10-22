CC = mpic++
EXEC = mpiexec
BIN_DIR = ./bin
HOSTS = ../hosts
OUTS := ${BIN_DIR}/Smooth

${BIN_DIR}/%: %.cpp
	mkdir -p $(BIN_DIR)
	$(CC) $< -o $@

build: $(OUTS)

test: $(OUTS)
	$(EXEC) -f $(HOSTS) -n 2 $<