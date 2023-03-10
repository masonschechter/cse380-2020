PDF := hw04.pdf
PDFLATEX := pdflatex -interaction nonstopmode
LOG := $(wildcard *.log)
AUX := $(wildcard *.aux)

$(PDF):
	$(PDFLATEX) main.tex

.PHONY: clean

clean:
	$(RM) $(AUX) $(LOG)
