PROJECT=signed_tensor_miccai
RST2LATEX=rst2latex
RUBBER=rubber

all: $(PROJECT).pdf

$(PROJECT).tex: $(PROJECT).rst
	$(RST2LATEX) $(PROJECT).rst $(PROJECT).tex

$(PROJECT).pdf: $(PROJECT).tex
	$(RUBBER) --pdf $(PROJECT).tex

