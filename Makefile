
all: doc classico

classico: tmp/classico.o
	@ mkdir -p bin
	$(LD) -l c -o bin/classico tmp/classico.o

doc: main.pdf

tmp/%.o: src/%.c
	@ mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $<

%.pdf: doc/%.tex
	@ mkdir -p tmp
	cd tmp && \
	TEXINPUTS=../src: BIBINPUTS=../src: latexmk \
	   -pdflatex \
	   -silent -quiet \
	   -f -g \
	   -logfilewarninglist \
	   ../$<

clean:
	rm -Rf tmp
