build.dataset.mammals <- function(){

    file <- "data/Ernest2003_Mammal_lifehistories_v2.txt"
    if(!file.exists(file))
        download.file("http://www.esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt", file)

    data <- read.table(file, header = TRUE, sep = "\t", strip.white = TRUE, blank.lines.skip = TRUE, na.strings = c("NA", "-999.00", "-999"))

    dry.mass.frac <- 1-0.605

    ernest <- data.frame(data[,c("order", "family", "Genus", "species")],
        W0 = data$wean.mass.g./1000 * dry.mass.frac,
        Wa =  data$mass.g./1000 * dry.mass.frac,
        A = 0.5*data$litter.size*data$litters.year*data$wean.mass.g./1000 * dry.mass.frac,
        stringsAsFactors=FALSE)
    ernest <- ernest[paste(ernest$Genus,ernest$species) != "Balaenoptera musculus", ]

    purvis  <- read.table("data/Purvis1995.txt", header = TRUE, sep = "\t")
    purvis$E <- 1/(1-exp(-purvis$M))/12
    mammals <- merge(ernest, purvis, c("Genus", "species"), all=TRUE)
    mammals <- mammals[!is.na(mammals$Wa*mammals$W0), ]

    write.table(mammals, "data/mammals.txt", quote=FALSE, row.names=FALSE, sep="\t")
}
