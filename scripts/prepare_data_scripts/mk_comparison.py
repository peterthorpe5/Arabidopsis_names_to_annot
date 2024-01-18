groups = set("""alkbh10c_high
alkbh10c_high
alkbh10c_high
alkbh10c_high
alkbh10c_low
alkbh10c_low
alkbh10c_low
alkbh10c_low
col0_high
col0_high
col0_high
col0_high
col0_low
col0_low
col0_low
col0_low
vir1_high
vir1_high
vir1_high
vir1_high
vir1_low
vir1_low
vir1_low
vir1_low""".split())

#################
#################

for i in groups:
    for j in groups:
        if i == j: continue
        #print("%s_vs_%s = condition%s - condition%s, " % (i, j, i, j))
        w = """# %s_vs_%s\n```{r}

# %s_vs_%s
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"%s_vs_%s"])

tTags = topTags(qlf,n=NULL)

result_table = tTags$table

result_table = data.frame(sampleA="%s", sampleB="%s", result_table)
result_table = merge(result_table, cts,by="row.names",all.x=TRUE)
result_table <- result_table[order(result_table$logFC),]




write.table(result_table, file='%s_vs_%s.GLM.edgeR.DE_results',
            sep='\t', quote=F, row.names=T, col.names = NA)

```\n""" % (i, j, i, j, i, j, i, j,  i, j)
        #print(w)

#my.contrasts <- makeContrasts(cortex_vs_striatum = conditioncortex - conditionstriatum, levels=design)



for i in groups:
    for j in groups:
        if i == j: continue
        #print("%s_vs_%s = condition%s - condition%s, " % (i, j, i, j))
        r = """\n```{r}
Pthreshold <- 0.05
logFC_threshold <- 1.6

# make a new column if the gene is sig yes or no
temp_data <- ''

temp_data <- read.table('results/%s_vs_%s.GLM.edgeR.DE_results', header=T, sep='\t')

temp_data$Significant<-ifelse(temp_data$PValue <=Pthreshold, 'Yes', 'No')

temp_data <- as_tibble(temp_data)

temp_data HHH>HHH 
  rename(gene_id = Row.names)
  
names(temp_data)[names(temp_data) == 'Row.names'] <- 'gene_id'

# assign ggplot to var called volc
volc <- NA

volc <- ggplot(data=temp_data, aes(logFC, -log10(PValue), label=gene_id)) +
   geom_point(aes(col=Significant))  +
  geom_text_repel(aes(label=ifelse(abs(logFC >= logFC_threshold) 
                       & PValue < Pthreshold, as.character(gene_id), ''))) +
  theme_grey() +
  geom_vline(xintercept=c(logFC_threshold,-logFC_threshold),
             linetype = 'dashed') +
  geom_hline(yintercept = -log10(Pthreshold), linetype = 'dashed') +
  scale_color_manual(values = c('black','red'))

# can save voc as a standard image if required 
# this is the interactive package
ggplotly(volc)

# save it to html
library(htmlwidgets)

g_plotly<-ggplotly(volc)
saveWidget(g_plotly, 'plots/%s_%s.GLM.edgeR.DE_results.html')



df <- temp_data HHH>HHH
  filter(Significant == 'Yes')

print('number DE genes:')
nrow(df)

write.table(df, file='results/%s_vs_%s.GLM.edgeR.LOGFC1.6_FDR0.05.txt',
            sep='\t', quote=F, row.names=T,col.names = NA)


```\n""" % (i, j, i, j, i, j)
        print(r)