SET  = set4
OPEN = xdg-open

CC = g++
FLAGS   = -g -ggdb -fPIC
# uncomment if linker flags are needed
#LFLAGS

$(SET).pdf: $(SET).tex Makefile molec.hpp Molec.py visualize.py
	pdflatex -shell-escape $(SET).tex
	pdflatex -shell-escape $(SET).tex

molec: molec.o
	$(CC) -o molec $(LFLAGS) molec.o

molec.so: molec.o
	$(CC) -shared -o $@ $^

molec.o: molec.cpp molec.hpp Makefile
	$(CC) -c $(FLAGS) molec.cpp

animated_gif/%.gif:
	sh mkgif.sh $(subst animated_gif/,,$(subst .gif,,$@))

.PHONY: clean view gifs
gifs: animated_gif/temp0.gif animated_gif/temp1.gif animated_gif/temp10.gif animated_gif/temp100.gif

view: $(SET).pdf
	$(OPEN) $(SET).pdf

clean:
	rm -f $(SET).pdf $(SET).log $(SET).aux
	rm -f molec molec.o molec.so
