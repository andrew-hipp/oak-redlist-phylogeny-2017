#redlist.tree.script.2017-05-15.R

library(ggtree)
library(caper)
library(phytools)
library(dismo)


do.pdf = TRUE
paint.branches = T
add.area = TRUE
do.phylo.signal = TRUE

## read data
tr <- read.tree('../data/tr.spp.4c.discreteClock.noOG.tre')
dat <- read.csv('../data/redlist.data.2017-04-07.csv', as.is = T, row.names = 1)
dat.trans <- read.csv('../data/redlist.codes.2017-04-07.csv', as.is = T, row.names = 1)
dat.coord <- read.csv('../data/specimenData_ahSubsetted_edited-2017-02-21b.csv', as.is = T)
dat.coord <- dat.coord[dat.coord$use, ]
dat.coord$Species <- paste('Quercus', dat.coord$Species)

## variables -- these can be edited to change what clades are colored and labeled
sections = list("White oaks\nsect. Quercus" = c('Quercus sadleriana', 'Quercus arizonica'),
                "Southern live oaks\nthe Virentes group" = c('Quercus fusiformis', 'Quercus minima'),
                "Intermediate oaks\nsect. Protobalanus" = c('Quercus palmeri', 'Quercus tomentella'),
                "Red oaks\nsect. Lobatae" = c('Quercus kelloggii', 'Quercus myrtifolia')
                )

section.colors = c('seagreen', 'seagreen3', 'sienna2', 'blue'); names(section.colors) <- names(sections)

## weld on taxon names -- actually, I'm not doing that now, so just tells you who's missing
whoIsMissing = setdiff(row.names(dat), tr$tip.label)
message('Missing these taxa:')
message(paste(sort(whoIsMissing), collapse = ', '))
message(paste("Total of", length(whoIsMissing), "taxa missing."))

## format data
tips <- intersect(tr$tip.label, row.names(dat))
tr.cleaned <- drop.tip(tr, which(!tr$tip.label %in% tips))
dat.cleaned <- data.frame("IUCN Red List Category" = dat[tr.cleaned$tip.label, ], row.names = tr.cleaned$tip.label)
row.names(dat.cleaned) <- tr.cleaned$tip.label <- gsub("_", " ", tr.cleaned$tip.label, fixed = T)
tr.cleaned$node.label <- rep(NA, tr.cleaned$Nnode + length(tr.cleaned$tip.label))
tr.mrca <- mrca(tr.cleaned)
tr.colors <- rep('black', tr.cleaned$Nnode + length(tr.cleaned$tip.label))
for(i in names(sections)) {
  tr.cleaned$node.label[tr.mrca[sections[[i]][1], sections[[i]][2]]] <- i
  if(paint.branches) tr.colors[getDescendants(tr.cleaned, tr.mrca[sections[[i]][1], sections[[i]][2]])] <- section.colors[i]
  }

## area
tips.range <- row.names(dat.cleaned)[which(row.names(dat.cleaned) %in% dat.coord$Species)]
ranges <- sapply(tips.range, function(x) raster:::area(convHull(dat.coord[dat.coord$Species == x, c('longitude', 'latitude')])@polygons))
latitudes <- sapply(tips.range, function(x) mean(dat.coord[dat.coord$Species == x, 'latitude']))
## plot data
if(do.pdf) pdf(format(Sys.time(), "redlist.phylogeny.%Y-%m-%d.16.pdf"), 4.5, 5.5) #4.5" for two columns, 2.125" for one column
par(mar = c(0,0,0,0))
p <- ggtree(tr.cleaned, layout = 'rectangular', size = 0.25, color = tr.colors)
p <- p + geom_tiplab(fontface='italic', size = 1.7)
p <- gheatmap(p, dat.cleaned, width = 0.025, font.size = 0, offset = ifelse(add.area, 15, 12))
p <- p + scale_fill_manual("IUCN Red List Category",
                  labels = dat.trans$category,
                   breaks=row.names(dat.trans),
                   values=structure(dat.trans$colors2, names = row.names(dat.trans)))
if(add.area) p <- p + geom_tippoint(size = (log(ranges[tr.cleaned$tip.label])-20)/5.5, position = position_nudge(x = 13.5))
p <- p + theme(legend.position = c(0.15,0.90),
       legend.title = element_text(size = 7),
       legend.text = element_text(size = 5),
       legend.key.size = unit(0.20, 'cm'),
       legend.box.background = element_rect(color = NA)
       )
      #p <- ggtree(tr.cleaned, layout = 'rectangular', size = 0.25, color = 'black')
p = p + geom_label(aes(x=branch),
                    label = tr.cleaned$node.label,
                    size = 1.7,
                    fill= section.colors[tr.cleaned$node.label],
                    label.padding = unit(0.18, "lines"),
                    label.r = unit(0.1, "lines"),
                    color = 'white')
print(p)
if(do.pdf) dev.off()

# phylogenetic signal
if(do.phylo.signal){
  area.cont <- phylosig(tr.cleaned, ranges, test = TRUE, nsim = 5000)
  area.cont.lambda <- phylosig(tr.cleaned, ranges, test = TRUE, method = 'lambda')
  dd = row.names(dat.cleaned)[dat.cleaned$IUCN.Red.List.Category == 'DD']
  if(length(dd) > 0){
    dat.cleaned <- dat.cleaned[which(!row.names(dat.cleaned) %in% dd), , drop = FALSE]
    tr.cleaned <- drop.tip(tr.cleaned, dd)
  }
  dat.cleaned$binary1 <- ifelse(dat.cleaned$IUCN.Red.List.Category %in% c('CR'), 1, 0)
  dat.cleaned$binary2 <- ifelse(dat.cleaned$IUCN.Red.List.Category %in% c('EN', 'CR'), 1, 0)
  dat.cleaned$binary3 <- ifelse(dat.cleaned$IUCN.Red.List.Category %in% c('EN', 'CR', 'VU'), 1, 0)
  dat.cleaned$binary4 <- ifelse(dat.cleaned$IUCN.Red.List.Category %in% c('EN', 'CR', 'VU', 'NT'), 1, 0)
  dat.cleaned$continuous <- 0
  dat.cleaned$continuous[which(dat.cleaned$IUCN.Red.List.Category == 'CR')] <- 5
  dat.cleaned$continuous[which(dat.cleaned$IUCN.Red.List.Category == 'EN')] <- 4
  dat.cleaned$continuous[which(dat.cleaned$IUCN.Red.List.Category == 'VU')] <- 3
  dat.cleaned$continuous[which(dat.cleaned$IUCN.Red.List.Category == 'NT')] <- 2
  dat.cleaned$continuous[which(dat.cleaned$IUCN.Red.List.Category == 'LC')] <- 1
  dat.cleaned$tip <- row.names(dat.cleaned)
  tr.cleaned$node.label <- NULL
  dat.bin <- comparative.data(tr.cleaned, dat.cleaned, tip)
  analysis.bin <- list(
    bin1 = phylo.d(dat = dat.cleaned, phy = tr.cleaned, names.col = tip, binvar = binary1, permut = 5000),
    bin2 = phylo.d(dat = dat.cleaned, phy = tr.cleaned, names.col = tip, binvar = binary2, permut = 5000),
    bin3 = phylo.d(dat = dat.cleaned, phy = tr.cleaned, names.col = tip, binvar = binary3, permut = 5000),
    bin4 = phylo.d(dat = dat.cleaned, phy = tr.cleaned, names.col = tip, binvar = binary4, permut = 5000),
    continuous = phylosig(tr.cleaned, dat.cleaned[tr.cleaned$tip.label, 'continuous'], test = TRUE, nsim = 5000),
    continuous.lambda = phylosig(tr.cleaned, dat.cleaned[tr.cleaned$tip.label, 'continuous'], test = TRUE, method = 'lambda'),
    area.continuous = area.cont,
    area.continuous.lambda = area.cont.lambda
  )
  sapply(analysis.bin, print)
}
