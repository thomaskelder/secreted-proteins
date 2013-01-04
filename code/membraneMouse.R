#################################################
## Generate list of membrane proteins in mouse ##
#################################################
require(org.Mm.eg.db)

#### Get from GO
plasmaMembraneTerm = "GO:0005886"
goEntrez = org.Mm.egGO2ALLEGS[["GO:0006629"]]

#### Get from Uniprot
## search terms: annotation:(type:location "plasma membrane") AND organism:"Mus musculus [10090]"
## http://www.uniprot.org/uniprot/?query=annotation%3a(type%3alocation+%22plasma+membrane%22)+AND+organism%3a%22Mus+musculus+%5b10090%5d%22&force=yes&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length
upMembrane = read.delim('data/uniprot-membrane-mouse.tab', as.is = T)
upEntrez = mget(upMembrane$Entry, revmap(org.Mm.egUNIPROT), ifnotfound=NA)
upEntrez = unique(unlist((upEntrez[!is.na(upEntrez)])))


#### Overlap
library(Vennerable)
v = Venn(list(Uniprot = upEntrez, GO = goEntrez))
plot(v)

#### Combine and save
membraneMouse = list(all = unique(c(goEntrez, upEntrez)), Uniprot = upEntrez, GO = goEntrez)
save(membraneMouse, file = 'data/membraneMouse.RData')