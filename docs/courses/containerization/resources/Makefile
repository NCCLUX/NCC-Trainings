export SHELL=/usr/bin/bash

TARGETS=$(patsubst sources/%.tex, %, $(wildcard sources/*.tex))
OUTPUT_RESOLUTION?=300

.PHONY: all
all: $(TARGETS)

.PHONY: $(TARGETS)
$(TARGETS): %: %.png

$(patsubst %, %.png, $(TARGETS)): %.png: sources/build/output/%.pdf
	convert -units PixelsPerInch -density $(OUTPUT_RESOLUTION) $< $@

# Add explicit dependency and ensure that PDF timestamp is updated when there is no change in the TEX file
sources/build/output/%.pdf: sources/%.tex
	$(MAKE) --directory=sources $(patsubst sources/build/output/%.pdf, %, $@)
	touch $@
sources/build/output/%.tex:
	;

.PHONY: clean deapclean
deapclean: clean
	$(MAKE) --directory=sources $<
clean:
	@echo "Removing files:"
	@find . -type f -print0 | \
		grep -zE '\.png$$' | \
		xargs -0I % bash -c '{ echo "%"; rm "%"; }'
