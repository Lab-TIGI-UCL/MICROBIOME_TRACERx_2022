library(dplyr)
library(reshape2)
library(cowplot)
library(gridExtra)
library(deconstructSigs)
library(parallel)
library(Biostrings)
library(RColorBrewer)
library(ggsci)

source('~/work/ucl/scripts/misc/functions.R')

make.key = function(df){
    paste(df$CHROM, df$POS, df$REF, df$ALT, sep = '_')
}

contam = read.csv('~/work/ucl/data/ucl.projects/microbiome/krupa.contamination.csv')

contam = contam %>%
    filter(Keep == 'NO')

contam

############################
#INITIAL DATA WRANGLING 
############################

taxa = read.tsv('~/work/ucl/data/ucl.projects/microbiome/rep82.tax', header  = F)
taxa$domain = multi.str.split(taxa$V2, ';', 1)
taxa$genus = multi.str.split(taxa$V2, ';', 6)

unique(taxa$genus)
head(taxa)
colnames(taxa)[1] = 'CHROM'
unique(taxa$domain)
taxa = taxa %>%
    filter(
        domain %in% c('k__Bacteria', 'k__Archaea',
            'k__BacteriaPlasmid', 'k__ArchaeaPlasmid')
    )

snps = read.tsv('~/work/ucl/data/ucl.projects/microbiome/new.callset.wgs.snps.w.tri.tsv')

colnames(snps) = c(
    'CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
    'FILTER',
    'INFO',
    'FORMAT',
    'CALLS',
    'sample',
    'tri'
)

############################
#FILTER SNPS IN VIRUSES & CONTAMINANTS
############################

contaminants = apply(contam, 1, function(x){
    unique(taxa$genus[grepl(x['Contaminant'], taxa$genus)])
}) %>% unlist

contaminants = taxa %>%
    filter(genus %in% contaminants) %>%
    pull(CHROM)

taxa = taxa %>%
    filter(!CHROM %in% contaminants) #contaminant filtering

snps = snps[snps$CHROM %in% taxa$CHROM, ]

############################
#MORE WRANGLING 
############################

snps$sample = gsub('.calls.txt$', '', snps$sample)
snps$key = make.key(snps)
snps$ad = gsub('.*?:.*?:.*?:.*?:.*?:(.*)', '\\1', snps$CALLS)

snps$dp = gsub('(DP=[0-9]*).*', '\\1', snps$INFO)
snps$dp = gsub('DP=', '', snps$dp)
snps$dp = as.numeric(snps$dp)

snps$mq = gsub('.*(MQ=[0-9]*).*', '\\1', snps$INFO)
snps$mq = gsub('MQ=', '', snps$mq)
snps$mq = as.numeric(snps$mq)

gel.inventory = read.tsv('~/work/ucl/data/ucl.projects/microbiome/ucl.gel100.inventory.tsv', header = F)
gel.inventory = gel.inventory[, 1:5]
colnames(gel.inventory)[1] = 'patient'
gel.inventory$sample = gsub('_[0-9]$', '', unique(gel.inventory$V5))

sampledb = read.csv('~/work/ucl/data/ucl.projects/microbiome/clinical_meta_data_tracerx.csv')

sampledb = sampledb[c('Hospital_ID', 'SmokedPastStop', 'SmokedPastStrt', 'DEMSmkYrs', 'DEMSmkPerDay',
    'SmokeAgeStrt', 'CurrentlySmoke', 'DEMSmkStatus.R_TRACERx_BaselineForm', 'PathHist.R_TRACERx_Lesion1Form')]
#sampledb = sampledb[, c('Hospital_ID', colnames(sampledb)[grepl('Smoke|Smk', colnames(sampledb))])]

sampledb = dplyr::rename(sampledb, hist = PathHist.R_TRACERx_Lesion1Form)


sampledb$patient = gsub('.*_', '', sampledb$Hospital_ID)

gel.inventory = left_join(gel.inventory, sampledb)
gel.inventory = distinct(gel.inventory, sample, .keep_all = T)

gel.inventory$DEMSmkYrs[is.na(gel.inventory$DEMSmkYrs)] = 0
gel.inventory$DEMSmkPerDay[is.na(gel.inventory$DEMSmkPerDay)] = 0

gel.inventory$DEMSmkYrs[gel.inventory$DEMSmkYrs == -9] = 1
gel.inventory$DEMSmkPerDay[gel.inventory$DEMSmkPerDay == -9] = 1

gel.inventory$packyear = gel.inventory$DEMSmkYrs * gel.inventory$DEMSmkPerDay

snps = left_join(snps, distinct(gel.inventory[c('V4', 'sample')]))
snps$patient = gsub('_.*$', '', snps$V4)
snps = dplyr::rename(snps, tumor.key = V4)

snps = snps[!is.na(snps$tumor.key), ] #remove samples where we don't have a mapping yet

############################
#GERMLINE FILTERING 
############################
# 1. filter SNPs found in germline samples
# 2. filter germline samples

snps$germline = F
snps$germline[grepl('_GL$', snps$tumor.key)] = T

germline.summary = snps %>%
    group_by(key) %>%
    summarise(
        n.gl = sum(germline == T),
        n.tumour = sum(germline == F),
        n = n()
    )

germline.summary = arrange(germline.summary, n)

germline.snps = germline.summary %>%
    filter(n.gl >= 1) %>%
    pull(key) %>%
    unique

germline.summary
germline.snps

#snps = snps[!snps$key %in% germline.snps, ] #filter GL SNPs

snps = snps %>% #filter GL samples
    filter(germline == F)

############################
#COLLAPSE REGIONS OF TUMOURS 
############################

gel.inventory$tumor.single = gsub('(_T1)-R[0-9]', '\\1', gel.inventory$V4)
snps$tumor.single = gsub('(_T1)-R[0-9]', '\\1', snps$tumor.key)

snps = snps %>%
    distinct(tumor.single, key, .keep_all = T)

sort(table(snps$tumor.single))

#snps = snps %>%
    #filter(tumor.single != 'LTX0033_T1')

############################
#RECURRANCE ANALYSIS
############################

g = snps %>%
    group_by(key) %>%
    summarise(
        n = n(),
        npatient = nunique(patient)
    )

snps.excl = snps %>%
    filter(tumor.key != 'LTX0033_T1-R4')

snps.freq = table(snps.excl$key)

df1 = data.frame(
    cumsum = cumsum(table(snps.freq)) / max(cumsum(table(snps.freq))),
    patient = 1:nunique(snps.excl$tumor.single)
)

plot1 = ggplot(df1, aes(x = patient, y = cumsum)) + 
    geom_point() +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black', size = 0.2)
    ) +
    ylab('Cumulative proportion of SNPs') +
    xlab('# of patients with SNP') +
    coord_cartesian(ylim = c(0, 1)) 
    #ggtitle('LTX0033 Excluded')

#plot with hypermutator
snps.freq = table(snps$key)
df2 = data.frame(
    cumsum = cumsum(table(snps.freq)) / max(cumsum(table(snps.freq))),
    patient = 1:nunique(snps$tumor.single)
)

plot2 = ggplot(df2, aes(x = patient, y = cumsum)) + 
    geom_point() +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black', size = 0.2)
    ) +
    ylab('Cumulative proportion of SNPs') +
    xlab('# of patients with SNP') +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle('LTX0033 Included')

#pl(plot1, 3, 3)

num.muts.multi = snps %>%
    group_by(tumor.key) %>%
    summarise(
        n = n()
    ) %>% arrange(n)

tail(num.muts.multi)

#the LTX0033 high mutations only occur in one region, not the entire tumour.
#it must be erroneous.
num.muts.multi %>%
    filter(grepl('LTX0033', tumor.key))

snps = snps %>%
    filter(tumor.key != 'LTX0033_T1-R4')



head(snps)

recurrent.ids = names(snps.freq)[snps.freq > 1]
snps$recurrent = F
snps$recurrent[snps$key %in% recurrent.ids] = T

head(gel.inventory)
head(snps)

snps = left_join(
    snps,
    distinct(gel.inventory[c('tumor.single', 'packyear')])
)

head(snps)
snps$smoker = 'ever'
snps$smoker[snps$packyear == 0] = 'never'

recurrent.by.smoker = snps %>%
    group_by(smoker) %>%
    summarise(
        num.recurrent = sum(recurrent),
        num.singleton = sum(!recurrent)
    ) 

recurrent.by.smoker

recurrent.by.smoker %>%
    summarise(
        smoker = smoker,
        prop = num.singleton / num.recurrent
    )

recurrent.by.smoker %>%
    select(num.recurrent, num.singleton) %>%
    chisq.test







############################
#HISTOLOGY BOXPLOT 
############################

num.muts = snps %>%
    group_by(tumor.single) %>%
    summarise(
        n = n()
    )

#head(as.data.frame(num.muts))
#head(gel.inventory)

num.muts = left_join(
    num.muts,
    distinct(gel.inventory[c('tumor.single', 'hist')])
) %>%
filter(hist %in% c('Invasive adenocarcinoma', 'Squamous cell carcinoma'))

hist.mut.plot = ggplot(num.muts, aes(x = hist, y = n)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black', size = 0.2)
    ) +
    rot.lab() +
    ylab('# Mutations') +
    xlab('Histology')




t.test(
    num.muts$n[num.muts$hist == 'Invasive adenocarcinoma'],
    num.muts$n[num.muts$hist == 'Squamous cell carcinoma']
)

############################
#GENERA MUTATION ANALYSIS 
############################

genera.muts = snps %>%
    group_by(CHROM) %>%
    summarise(
        n = n()
    ) %>%
    arrange(n)

as.data.frame(genera.muts)


genera.muts = left_join(
    genera.muts,
    distinct(taxa[c('CHROM', 'genus')])
)

genera.muts
tail(genera.muts, n = 10)

cutoff1 = tail(genera.muts, n = 10) %>% pull(n) %>% min

genera.muts$genus[genera.muts$n < cutoff1] = 'Other'
genera.muts = dplyr::rename(genera.muts, Genus = genus)
genera.muts$Genus = gsub('^g__', '', genera.muts$Genus)

other.muts = genera.muts %>%
    filter(Genus == 'Other') %>%
    pull(n) %>%
    sum()

genera.muts = genera.muts %>%
    filter(Genus != 'Other') %>%
    bind_rows(data.frame(CHROM = '', n = other.muts, Genus = 'Other'))

levels1 = sort(genera.muts$Genus)
levels1 = unique(c(levels1[!levels1 == 'Other'], 'Other'))

genera.muts$Genus = factor(genera.muts$Genus, levels = levels1) 

pieplot = ggplot(genera.muts, aes(x = "", y = n, fill = Genus)) +
    geom_bar(stat = 'identity', width = 1) +
    coord_polar('y', start = 0) +
    theme_void() +
    scale_fill_ucscgb()

############################
#FILTER RECURRENT SNPS 
############################

#currently filtering where SNP is present in more than 1 sample
recurrent.ids = names(snps.freq)[snps.freq > 2]
snps$recurrent = F
snps$recurrent[snps$key %in% recurrent.ids] = T
snps = snps[snps$recurrent == F, ] # filter recurrent SNPs

############################
#PARSING 
############################

num.genera = with(snps, table(sample, CHROM))
num.genera
apply(num.genera, 1, sum)

mq.thres = 20
dp.thres = 10

############################
#FURTHER PROCESSING 
############################

#snps = snps[snps$dp > dp.thres & snps$mq > mq.thres & snps$recurrent == F, ]
#snps = snps[!is.na(snps$dp), ]


############################
#PROCESS TRINUC 
############################

snps = snps[!nchar(snps$REF) > 1, ]
snps = snps[!nchar(snps$ALT) > 1, ]
snps = snps[!grepl('N', snps$tri), ]
snps = snps[!grepl('M|Y|W|R|K|S', snps$tri), ]

snps.g.a = snps[grepl('G|A', snps$REF), ]
snps.t.c = snps[grepl('T|C', snps$REF), ]

snps.g.a$tri = as.character(reverseComplement(DNAStringSet(snps.g.a$tri)))
snps.g.a$REF = as.character(reverseComplement(DNAStringSet(snps.g.a$REF)))
snps.g.a$ALT = as.character(reverseComplement(DNAStringSet(snps.g.a$ALT)))

snps = bind_rows(snps.t.c, snps.g.a)

data(signatures.cosmic)
data(sample.mut.ref)

dcs.example = mut.to.sigs.input(sample.mut.ref)

snps$tri2 = paste0(
    substr(snps$tri, 1, 1),
    '[',
    snps$REF,
    ">",
    snps$ALT,
    ']',
    substr(snps$tri, 3, 3)
)

snps$tri2 = factor(snps$tri2, levels = colnames(dcs.example))
mb.sigs = as.data.frame(rbind(with(snps, table(tumor.single, tri2))))
mb.sigs.df = as.data.frame(mb.sigs)
mb.sigs.df$sample = rownames(mb.sigs.df)

mb.sigs.df.melt = melt(mb.sigs.df)

mb.sigs.plot = ggplot(mb.sigs.df.melt, aes(x = variable, y = value)) + geom_col() +
    facet_grid(rows = vars(sample), scales = 'free') +
    rot.lab(size = 4)

#pl(mb.sigs.plot, 10, 40)

############################
#SPECTRA PLOTTING 
############################

head(mb.sigs.df.melt)
mb.sigs.df.melt$patient = gsub('_.*', '', mb.sigs.df.melt$sample)

mb.sigs.df.melt = left_join(mb.sigs.df.melt, distinct(gel.inventory[c('patient', 'packyear')]))
mb.sigs.df.melt$group = 'Never Smokers'
mb.sigs.df.melt$group[mb.sigs.df.melt$packyear > 0] = 'Ever Smokers'
mb.sigs.df.melt$sub = gsub('[A-Z]\\[(.*)\\][A-Z]', '\\1', mb.sigs.df.melt$variable)
mb.sigs.df.melt$variable2 = 
    gsub('([A-Z])\\[([A-Z]).*\\]([A-Z])', '\\1\\2\\3', mb.sigs.df.melt$variable)

mb.sigs.df.melt = mb.sigs.df.melt %>% dplyr::rename(Substitution = sub)

head(mb.sigs.df.melt)

mb.prop.all = mb.sigs.df.melt %>%
    group_by(variable) %>%
    summarise(
        n.mut = sum(value)
    ) %>%
    mutate(
        prop = n.mut / sum(n.mut)
    )

mb.prop = mb.sigs.df.melt %>%
    group_by(variable, group) %>%
    summarise(
        n.mut = sum(value)
    ) %>%
    arrange(group, variable)


mb.split = split(mb.prop, mb.prop$group)

mb.prop2 = lapply(mb.split, function(x){
    x$prop = x$n.mut / sum(x$n.mut)
    x
}) %>% bind_cols %>% reset.colnames

all(mb.prop2$V1 == mb.prop2$V5)

mb.prop2$diff = mb.prop2$V4 - mb.prop2$V8
head(mb.prop2)
mb.diff = mb.prop2[c('V1', 'V2', 'diff')]

mb.prop.all$cat = 'Combined Spectrum'
mb.diff$cat = 'Differential Spectrum'

mb.prop.all = mb.prop.all[c('variable', 'prop', 'cat')]
mb.diff = mb.diff[c('V1', 'diff', 'cat')]
colnames(mb.diff) = c('variable', 'prop', 'cat')

mb.revised = bind_rows(mb.prop.all, mb.diff)
mb.revised$sub = gsub('[A-Z]\\[(.*)\\][A-Z]', '\\1', mb.revised$variable)
mb.revised$variable2 = 
    gsub('([A-Z])\\[([A-Z]).*\\]([A-Z])', '\\1\\2\\3', mb.revised$variable)
mb.revised = mb.revised %>% dplyr::rename(Substitution = sub)


never.spec.plot = ggplot(mb.revised, aes(x = variable2, y = prop, fill = Substitution)) +
    geom_col() +
    facet_grid(rows = vars(cat), cols = vars(Substitution), scales = 'free') +
    theme_bw() +
    rot.lab(size = 8) +
    ylab('# Mutations') +
    xlab('Trinucleotide Context') +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black', size = 0.2),
        panel.margin.x = unit(0, 'lines'),
        strip.background = element_rect(fill = 'grey90', color = 'grey70')
    ) 
    #scale_fill_brewer(palette = 'Set1')

#pl(never.spec.plot, 17, 6)

############################
#ALL PLOTS 
############################

g = plot_grid(plot1, hist.mut.plot, ncol = 2)
g2 = plot_grid(never.spec.plot, g, ncol = 1, axis = 'lr', align = 'hv')

pl(g2, 16, 9)

head(snps)

snps %>%
    group_by(smoker) %>%
    summarise(
        n = n()
    )



#snps
#contaminants


#con.chrom = taxa %>%
    #filter(genus %in% contaminants) %>%
    #pull(CHROM)

#snps[which(snps$CHROM %in% con.chrom), ]

browser()


#never.spec.plot = ggplot(mb.sigs.df.melt, aes(x = variable2, y = value, fill = Substitution)) +
    #geom_col() +
    #facet_grid(rows = vars(group), cols = vars(Substitution), scales = 'free') +
    #theme_bw() +
    #rot.lab(size = 8) +
    #ylab('# Mutations') +
    #xlab('Trinucleotide Context') +
    #theme(
        #panel.grid = element_blank(),
        #panel.border = element_blank(),
        #axis.line = element_line(color = 'black', size = 0.2),
        #panel.margin.x = unit(0, 'lines'),
        #strip.background = element_rect(fill = 'grey90', color = 'grey70')
    #) 
    ##scale_fill_brewer(palette = 'Set1')

pl(never.spec.plot, 17, 6)



left_join(mb.sigs.df.melt, gel.inventory[c('patient', '')])
head(gel.inventory)


num.muts = mb.sigs.df.melt %>%
    group_by(sample) %>%
    summarise(
        num.mut = sum(value)
    )  %>%
    arrange(num.mut)

de.sigs.all = mclapply(unique(snps$tumor.single), function(x){
    print(paste0('RUNNING ', x))

    test.sigs = whichSignatures(
        tumor.ref = as.data.frame(mb.sigs),
        signatures.ref = signatures.nature2013,
        sample.id = x,
        contexts.needed = T,
        tri.counts.method = 'default'
    )
}, mc.cores = 4)

de.sigs.all[[1]]
sigs.weights = lapply(de.sigs.all, function(x){
    x$weights
})

sigs.weights2 = bind_rows(sigs.weights)

sigs.weights2


sigs.weights2 = sigs.weights2[order(sigs.weights2$Signature.4, decreasing = T), ]
sigs.weights2$sample = rownames(sigs.weights2)

sig4.samps = sigs.weights2 %>%
    select(sample, Signature.4) %>%
    reset.rownames %>%
    tibble

sig4.samps2 = left_join(
    sig4.samps,
    distinct(
        gel.inventory[c('sample', 'V4')]
    )
)

as.data.frame(sig4.samps2)
head(sig4.samps2)






head(gel.inventory)

dim(gel.inventory)
head(sig4.samps)

sig4.samps3 = left_join(sig4.samps, gel.inventory) %>% left_join(num.muts)

sig4.samps3
num.muts

wt(sig4.samps3)
sort(unique(gel.inventory$sample))






sig.sum = t(as.data.frame(lapply(sigs.weights2, sum)))
sig.sum = as.data.frame(sig.sum)
sig.sum$sig = rownames(sig.sum)
sig.sum = reset.rownames(sig.sum)

sig.sum$sig = factor(sig.sum$sig, levels = unique(sig.sum$sig))

sig.plot = ggplot(sig.sum, aes(x = sig, y = V1)) + geom_col() + rot.lab()
pl(sig.plot)


############################
#COMPARE TRANSCRIPTOME TO WGS 
############################

#bwa.stats = slp(snps.bwa, 'sample', function(x){
    #data.frame(
        #num.snps = nrow(x),
        #sample = unique(x$sample)
    #)
#})


#wgs.stats = slp(snps, 'sample', function(x){
    #data.frame(
        #num.snps = nrow(x),
        #sample = unique(x$sample)
    #)
#})

#mean(bwa.stats$num.snps)
#sd(bwa.stats$num.snps)

#mean(wgs.stats$num.snps)
#sd(wgs.stats$num.snps)


############################
#MORE PLOTTING 
############################

sig4df = as.data.frame((t(as.data.frame((signatures.nature2013['Signature.4', ])))))
rownames(sig4df)
sig4df$tri = rownames(sig4df)
head(sig4df)
#head(wgs.df)
#sort(as.character(wgs.df$Var1))
colnames(signatures.nature2013)
rownames(sig4df)
sig4df$tri
sig4df$tri = factor(sig4df$tri, levels = unique(sig4df$tri))

plot1 = ggplot(sig4df, aes(x = tri, y = Signature.4)) + geom_col() +
    rot.lab(size = 4) +
    ggtitle('Signature 4 trinuc')
#pl(plot1)

wgs.snps.all = table(snps$tri2)
wgs.df = as.data.frame(wgs.snps.all)
head(wgs.df)

wgs.df$Var1 = factor(wgs.df$Var1, levels = colnames(signatures.nature2013))

plot2 = ggplot(wgs.df, aes(x = Var1, y = Freq)) + geom_col() +
    rot.lab(size = 4) +
    ggtitle('WGS trinuc')
#pl(plot2)

plot3 = plot_grid(plot1, plot2, align = 'v', axis = 'lr', ncol = 1)
#pl(plot3)
#pl(grid.arrange(plot1, plot2, ncol = 1))


head(mb.sigs)
mb.sigs.df = as.data.frame(mb.sigs)
mb.sigs.df$sample = rownames(mb.sigs.df)

head(mb.sigs.df)
mb.sigs.df.melt = melt(mb.sigs.df)
head(mb.sigs.df.melt)

mb.sigs.plot = ggplot(mb.sigs.df.melt, aes(x = variable, y = value)) +
    geom_col() +
    facet_grid(rows = vars(sample), space = 'free', scales = 'free') +
    rot.lab(size = 4)

#pl(mb.sigs.plot, 10, 200)


############################
#GENERATE SAMTOOLS TVIEW COMMANDS 
############################

snps$tview = paste0(
    '$samtools tview -p ',
    snps$CHROM,
    ':',
    snps$POS,
    ' ',
    paste0(snps$sample, '.bam.unmapped.txt.fastq.sam.sorted.bam'),
    ' /cluster/scratchpo1/acoulton/shogun.ref/rep82.fna'
)


snps$tview = paste0(
    '$samtools tview -p ',
    snps$CHROM,
    ':',
    snps$POS,
    ' ',
    paste0(snps$sample, '.bam.unmapped.txt.fastq.sam.sorted.bam'),
    ' /cluster/scratchpo1/acoulton/shogun.ref/rep82.fna'
)

wt(snps)


############################
#FILTERING BY 16S 
############################

taxa.map = read.tsv('~/work/ucl/data/ucl.projects/microbiome/rep82.tax', header = F)
taxa.map = taxa.map[!grepl('Viruses', taxa.map$V2), ]

data16s = read.csv('~/work/ucl/data/ucl.projects/microbiome/new_samples_data_16s.csv')
data16s = data16s[!apply(data16s, 1, function(x){
    all(as.numeric(x[2:length(x)]) == 0)
}), ]

data16s = data16s[order(apply(data16s, 1, function(x){
    sum(as.numeric(x[2:length(x)]))
}), decreasing = T), ]

#get the reference genomes (chromosome ids) for each genus
top20genera = lapply(data16s[1:20, ]$X, function(x){
    g = taxa.map[grepl(x, taxa.map$V2), ]
    if(nrow(g) > 0) g$genus = x
    g
}) 

summary.snps = lapply(top20genera, function(x){
    genus = unique(x$genus)
    genus.snps = snps[snps$CHROM %in% x$V1, ]
    genus.snps = genus.snps[genus.snps$dp > 20, ]
    genus.snps = genus.snps[genus.snps$mq > 20, ]
    if(nrow(genus.snps) > 0){
        data.frame(
            genus = genus,
            num.snps = nrow(genus.snps),
            mean.dp = mean(genus.snps$dp, na.rm = T),
            sd.dp = sd(genus.snps$dp, na.rm = T),
            mean.mq = mean(genus.snps$mq, na.rm = T)
        )
    } else {
        data.frame(
            genus = NA,
            num.snps = NA,
            mean.dp = NA,
            sd.dp = NA,
            mean.mq = NA
        )
    }
}) %>% bind_rows
summary.snps

#genera.snps = lapply(top20genera, function(x){
    #genus = unique(x$genus)
    #genus.snps = snps[snps$CHROM %in% x$V1, ]
    #genus.snps
#}) %>% bind_rows

#mb.sigs = matrix(
    #nrow = nunique(genera.snps$sample),
    #ncol = ncol(dcs.example)
#)

#mb.sigs[] = 0
#colnames(mb.sigs) = colnames(dcs.example)
#rownames(mb.sigs) = unique(genera.snps$sample)


#count = 1
#apply(genera.snps, 1, function(x){
    ##if(count == 232) browser()
    #mb.sigs[x['sample'], ]
    #mb.sigs[x['sample'], x['tri2']] <<- 
        #mb.sigs[x['sample'], x['tri2']] + 1
    #count <<- count + 1
#})

#head(genera.snps)

#all.sigs = lapply(unique(genera.snps$sample), function(x){
    #genera.sigs = whichSignatures(
        #tumor.ref = as.data.frame(mb.sigs),
        #signatures.ref = signatures.nature2013,
        #sample.id = x,
        #contexts.needed = T,
        #tri.counts.method = 'default'
    #)
#})

#lapply(all.sigs, function(x) x$weights)

#mb.sigs.df = as.data.frame(mb.sigs)
#mb.sigs.df$sample = rownames(mb.sigs.df)

#head(mb.sigs.df)
#mb.sigs.df.melt = melt(mb.sigs.df)
#head(mb.sigs.df.melt)

#mb.sigs.plot = ggplot(mb.sigs.df.melt, aes(x = variable, y = value)) +
    #geom_col() +
    #facet_grid(rows = vars(sample), space = 'free', scales = 'free') +
    #rot.lab(size = 4)

#pl(mb.sigs.plot)

############################
#POSITIVE CONTROL 
############################

#human.sigs = read.delim('~/work/ucl/data/ucl.projects/microbiome/michelled/SBS_includingSBS92/signature_weights_perTumour.txt', sep = ' ', header = T)

#head(human.sigs)

#head(human.sigs[order(human.sigs$SBS4, decreasing = T), ], n = 20)


############################
#ASSOCIATE TOP SIG 4 MUTS WITH PACK YEAR
############################


head(snps)
snps4 = snps[snps$tri2 %in% c(
    'C[C>A]A',
    'C[C>A]C',
    'C[C>A]T',
    'T[C>A]C'
), ]

snps4.split = split(snps4, snps4$sample)







mutation.distribution = lapply(snps4.split, function(x){
    data.frame(
        num.mut = nrow(x),
        sample = unique(x$sample)
    )
}) %>% bind_rows

mutation.distribution = arrange(mutation.distribution, num.mut)

#wt(sampledb)

samples = read.tsv('~/work/ucl/data/ucl.projects/microbiome/ucl.gel100.inventory.tsv', header = F)
samples = samples[c('V5', 'V4')]
samples$V5 = gsub('_[0-9]$', '', samples$V5)
colnames(samples) = c('sample', 'patient.id')

head(samples)
samples$patient = multi.str.split(samples$patient.id, '_', 1)
samples
head(sampledb$Hospital_ID)
sampledb$patient = multi.str.split(sampledb$Hospital_ID, '_', 2)
sampledb = sampledb[sampledb$patient %in% samples$patient, ]

sampledb

mutation.distribution$patient = samples$patient.id[match(mutation.distribution$sample, samples$sample)]
mutation.distribution

mutation.distribution2 = mutation.distribution[!is.na(mutation.distribution$patient), ]
mutation.distribution2
mutation.distribution2$patient.short = multi.str.split(mutation.distribution2$patient, '_', 1)

head(sampledb)
sampledb
mut2 = mutation.distribution2
mut2$smk = sampledb$DEMSmkYrs[match(mut2$patient.short, sampledb$patient)]
mut2$smk.per.day = sampledb$DEMSmkPerDay[match(mut2$patient.short, sampledb$patient)]
mut2 = mut2[!grepl('GL', mut2$patient), ]

mut2
mut2$smk[mut2$smk < 0] = 0
mut2$smk[is.na(mut2$smk)] = 0
mut2$smk.per.day[mut2$smk.per.day < 0] = 0
mut2$smk.per.day[is.na(mut2$smk.per.day)] = 0
head(mut2)
mut2$weighted.smk = mut2$smk * mut2$smk.per.day

dim(mut2)
nunique(mut2$patient.short)
#pl(plot(mut2$num.mut, mut2$smk))
#pl(plot(mut2$num.mut, mut2$weighted.smk))

mut2.split = split(mut2, mut2$patient)

# Are the same mutations shared across samples?
snps$key = paste0(snps$CHROM, '_', snps$POS, '_', snps$REF, '_', snps$ALT)
key.split = split(snps, snps$key)
# No, not really, each mutation seems to be unique to an individual sample, except for a 
# small minority of mutations that are present in many samples
#pl(hist(sort(uulapply(key.split, nrow))))



#pl(hist(mut2$num.mut, breaks = 50))

###########################
#INVESTIGATION INTO ~200 SNP SAMPLE
###########################

s5 = snps4.split[order(uulapply(snps4.split, nrow), decreasing = T)]

snps4.summary = lapply(snps4.split, function(x){
    data.frame(
        num.snps = nrow(x),
        sample = unique(x$sample)
    )
}) %>% bind_rows

s5.tab = lapply(s5, function(x) table(x$CHROM))

table(s5[[1]]$CHROM)
head(s5[[1]])

#wt(s5[[1]])

#write.tsv(s5[[1]], '~/work/ucl/data/ucl.projects/microbiome/LP3000893-DNA_C01.snps')

uulapply(s5, function(x) length(which(x$CHROM == 'NZ_CP008918.1')))

taxa.map[taxa.map$V1 == 'NZ_CP008918.1', ]

length(grep('Pasteurella', taxa.map$V2))

unique(s5[[1]]$tri2)

arrange(snps4.summary, num.snps)

snps.review = slp(snps, "sample", function(x){
    data.frame(
        num.snps = nrow(x),
        sample = unique(x$sample)
    )
})

snps.review = arrange(snps.review, num.snps)
snps.review
arrange(snps4.summary, num.snps)

snps.comb = left_join(snps.review, snps4.summary, by = 'sample')
head(snps.comb)

plot.review = ggplot(snps.comb, aes(x = num.snps.x, y = num.snps.y)) +
    geom_point() +
    xlab('total snps') +
    ylab('top C->A snps')

#pl(plot.review, 5, 5)
