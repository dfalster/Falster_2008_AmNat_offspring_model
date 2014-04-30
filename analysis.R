# Daniel Falster - offspring size model

library(smatr)

source("R/data-fun.R")
source("R/model-fun.R")
source("R/utils.R")
source("R/Figs.R")

build.dataset.mammals()

if(!file.exists("output"))
	dir.create("output")

to.pdf(Fig2(), "output/fig2.pdf", height = 6, width = 10, pointsize = 10)
to.pdf(Fig3(), "output/fig3.pdf", height = 6, width = 6, pointsize = 10)
to.pdf(Fig4(), "output/fig4.pdf", height = 6, width = 6, pointsize = 10)
to.pdf(Fig5(), "output/fig5.pdf", height = 6, width = 6, pointsize = 10)
to.pdf(Fig6(), "output/fig6.pdf", height = 6, width = 6, pointsize = 10)
