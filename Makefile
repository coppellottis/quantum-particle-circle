# Compilatore

CC = gcc

# Opzioni compilazione

CFLAGS = -Wall -Wextra -O2 -Iinclude -fopenmp

# -Wall : attiva tutti i warning
# -Wextra : aggiunge altri controlli
# -O2 : ottimizzazione del codice
# -Iinclude : cerca i file .h nella cartella /include/

# Nome eseguibile

TARGET = quantum

# File sorgenti

SRC = src/main.c \
      src/io.c \
      src/polymer.c \
      src/pcg_basic.c \
      src/metropolis.c \
      src/observables.c \
      src/simulation.c \
      src/progress.c

# Cartella degli oggetti

BUILDDIR = build

# Conversione src/*.c -> build/*.o

OBJ = $(patsubst src/%.c,$(BUILDDIR)/%.o,$(SRC))

# Regola principale: link degli oggetti nell'eseguibile

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) -lm -fopenmp
# -lm linka la libreria matematica <math.h>

# Compilazione dei singoli file .c -> .o

$(BUILDDIR)/%.o: src/%.c
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Target opzionale con ottimizzazioni piu' aggressive, specifiche per la
# macchina su cui si compila (-march=native non e' portabile: il binario
# risultante va eseguito sulla stessa macchina/microarchitettura). Usare
# "make fast" invece di "make" quando si vuole spremere le simulazioni piu'
# lunghe e non serve distribuire l'eseguibile altrove.
fast: CFLAGS = -Wall -Wextra -O3 -march=native -flto -Iinclude -fopenmp
fast: clean $(TARGET)

# Pulizia

clean:
	rm -rf $(BUILDDIR) $(TARGET)

.PHONY: clean fast