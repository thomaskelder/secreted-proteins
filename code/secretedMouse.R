#################################################
## Generate list of secreted proteins in mouse ##
#################################################
require(org.Mm.eg.db)

## SPD (http://spd.cbi.pku.edu.cn/)
spdData = read.delim('data/spd.name.nr90.species-20121016.txt', as.is = T, header = F, col.names = c("id", "name", "rank", "species"))
spdNames = spdData[spdData$species == 'mouse', 'name']
spdNames = gsub("(^.{1})(.+)", "\\U\\1\\L\\2", spdNames, perl=TRUE)
spdEntrez = mget(spdNames, revmap(org.Mm.egSYMBOL), ifnotfound=NA)
spdEntrez = unique(as.character(spdEntrez[!is.na(spdEntrez)]))

## Sys-BodyFluid (http://www.biosino.org/bodyfluid/)
sbfData = read.delim('data/sysBodyFluid-20121016.txt', as.is=T, header = F, sep=' ')
sbfData = sbfData[, 7:8]
colnames(sbfData) = c("IPI", "PubMed")
require(org.Hs.ipi.db)
ipiToEntrez = as.list(org.Hs.ipiGENEID)
names(ipiToEntrez) = gsub("\\.[0-9]+$", "", names(ipiToEntrez))
sbfEntrezHs = ipiToEntrez[sbfData$IPI]
sbfEntrezHs = unique(as.character(unlist(sbfEntrezHs[!is.null(sbfEntrezHs) & !is.na(sbfEntrezHs)])))
library("hom.Hs.inp.db")
library("org.Mm.eg.db")
library("org.Hs.eg.db")
sbfEntrez = mapEntrezToSpecies(sbfEntrezHs, inProt = org.Hs.egENSEMBLPROT, outProt = org.Mm.egENSEMBLPROT2EG, homology = hom.Hs.inpMUSMU)
#inProt = org.Hs.egENSEMBLPROT, outProt = org.Mm.egENSEMBLPROT2EG, hom = hom.Hs.inpMUSMU


## Uniprot keyword secreted:
## organism:Mus musculus keyword:KW-0964 NOT keyword:KW-0134 NOT keyword:KW-0272 NOT keyword:KW-0767
upSecr = read.delim('data/uniprot-secreted-mouse.tab', as.is = T)
upEntrez = mget(upSecr$Entry, revmap(org.Mm.egUNIPROT), ifnotfound=NA)
upEntrez = unique(unlist((upEntrez[!is.na(upEntrez)])))

library(Vennerable)
v = Venn(list(Uniprot = upEntrez, SPD = spdEntrez, SBF = sbfEntrez))
plot(v)

secretedMouse = list(all = unique(c(upEntrez, spdEntrez, sbfEntrez)), Uniprot = upEntrez, SPD = spdEntrez, SBF = sbfEntrez)

save(secretedMouse, file = 'data/secretedMouse.RData')