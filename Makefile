SET  = set4
OPEN = xdg-open

FLAGS = -g -ggdb

$(SET).pdf: $(SET).tex Makefile
	pdflatex -shell-escape $(SET).tex
	pdflatex -shell-escape $(SET).tex

molec: molec.o
	g++ -o molec molec.o

molec.o: molec.cpp molec.hpp
	g++ -c $(FLAGS) molec.cpp

.PHONY: clean view
view: $(SET).pdf
	$(OPEN) $(SET).pdf

clean:
	rm -f $(SET).pdf $(SET).log $(SET).aux
	rm -f molec molec.o
