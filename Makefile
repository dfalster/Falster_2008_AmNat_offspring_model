all: output/fig2.pdf

output/fig2.pdf: analysis.R
	Rscript -e "source('$<')"

clean:
	rm -r output

.PHONY: all clean
