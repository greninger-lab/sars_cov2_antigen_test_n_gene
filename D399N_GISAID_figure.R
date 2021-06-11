library(ggplot2)
library(tidyverse)
library(dplyr)
library(maps)
library(viridis)
library(stringr)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(plyr)
library(UpSetR)

setwd("/Users/gerbix/Documents/Michelle/swift_testing/quidel_discordant_fastas/")

# Generated from GISAID master list of sequences + variants, pulled on 4/18/21
D399N_gisaid <- read_tsv("variant_surveillance.tsv") %>% filter(grepl("D399N",`AA Substitutions`))

# Format GISAID-specified locations into usable format
D399N_gisaid <- separate(D399N_gisaid, Location, into=c("continent","region","subregion"), sep=" / ") %>%
                mutate(region = recode(region,"United Kingdom" = "UK"), subregion = recode(subregion,"England" = "Great Britain"))

# Get world map
world_map <- map_data("world") %>% filter(region!="Antarctica")

# Find mean latitude and longitudes for each subregion
region_coords <- aggregate(cbind(long, lat) ~ region, data=world_map, 
                     FUN=function(x)mean(range(x)))
subregion_coords <- aggregate(cbind(long, lat) ~ subregion, data=world_map, 
                        FUN=function(x)mean(range(x)))
colnames(region_coords) <- c("region","mean_long","mean_lat")
colnames(subregion_coords) <- c("subregion","mean_long","mean_lat")

## Do US and UK separately for more subregion data
cleaned_uk_map <- world_map %>% filter(region=="UK")
cleaned_uk_map <- merge(cleaned_uk_map, subregion_coords, by=c("subregion"))

cleaned_us_map <- map_data("state")
cleaned_us_map$subregion <- str_to_title(cleaned_us_map$region)
cleaned_us_map <- merge(cleaned_us_map, subregion_coords, by="subregion")
cleaned_us_map$region <- "USA"

# Do the rest of the world by region
cleaned_world_map <- world_map %>% filter(region!="USA" & region!="UK")
cleaned_world_map <- merge(cleaned_world_map, region_coords, by="region")


# Now merge with our variants data
D399N_gisaid_uk <- D399N_gisaid %>% filter(region=="UK")
D399N_gisaid_us <- D399N_gisaid %>% filter(region=="USA")

#D399N_gisaid_noukus <- D399N_gisaid %>% filter(region!="USA" & region!="UK")
D399N_geo_uk <- inner_join(cleaned_uk_map, D399N_gisaid_uk, by=c("region","subregion"))
D399N_geo_us <- inner_join(cleaned_us_map, D399N_gisaid_us, by=c("region","subregion"))
D399N_geo <- subset(inner_join(cleaned_world_map, D399N_gisaid, by=c("region")), select=-c(subregion.x))
names(D399N_geo)[names(D399N_geo)=="subregion.y"] <- "subregion"
D399N_geo_final <- rbind(D399N_geo, D399N_geo_uk,D399N_geo_us)
D399N_geo_final <- D399N_geo_final %>% group_by(`Accession ID`) %>% top_n(lat, n=1) %>% top_n(long,n=1)

# Make a color scale for clades
myColors <- c(brewer.pal(7, "Set2"),brewer.pal(8,"Accent"),brewer.pal(7, "Dark2"),brewer.pal(8,"Set3"))
names(myColors) <- levels(D399N_geo_final$`Pango lineage`)
colScale <- scale_color_manual(name = , values=myColors)

## A: Plot geographic locations of GISAID consensuses with D399N
D399N_geo_final_loc <- D399N_geo_final %>% group_by(region,subregion,mean_lat,mean_long,`Pango lineage`) %>% 
                            dplyr::summarise(num_samples = n_distinct(`Accession ID`))
worldplot <- ggplot(data = D399N_geo_final_loc, 
                    aes(x=mean_long, y=mean_lat, color=`Pango lineage`, size=num_samples)) + 
              #scale_size_continuous(limits=c(1,13),breaks=c(1,5,9,13)) + 
              scale_size_continuous(limits=c(1,30),breaks=c(1,10,20,30)) + 
              #scale_color_brewer(palette = "Dark2", direction=-1) + 
              #scale_color_brewer(palette = "Set2",direction=-1) + 
              colScale + 
            geom_polygon(data = world_map, aes(x=long, y= lat, group= group),fill="lightgray",colour="white", size=0.1) + 
            coord_fixed(1.3)  + 
            theme_void() + scale_x_continuous(limits=c(-170,175)) + geom_jitter(alpha=0.6, width=1, height=1) + 
            theme(legend.box="vertical",legend.position="bottom", legend.key.size = unit(0.7, "lines")) + labs(color="Lineage",size="# Samples")
            


# Format collection date
D399N_geo_final$`Collection date` <- as.POSIXct.Date(as.Date(D399N_geo_final$`Collection date`))
# Some regions don't have subregions (e.g. Japan) - just replace these with region names
D399N_geo_final$subregion <- ifelse(is.na(D399N_geo_final$subregion), D399N_geo_final$region, D399N_geo_final$subregion)

## B: Plot GISAID consensuses with D399N over time
timemuts <- ggplot(data=D399N_geo_final, aes(x=`Collection date`,y=subregion, color=`Pango lineage`)) + #color=`Pango lineage` 
        #scale_color_brewer(palette="Dark2",direction=-1) +
        colScale + 
        geom_point() + facet_grid(rows = vars(region), scales="free",switch="y",space="free_y") +
          scale_x_datetime(date_breaks = "1 month") + theme_clean() + 
      theme(strip.placement="outside",strip.text.y.left = element_text(angle=0), strip.background.y = element_rect(fill="gray90"),
            axis.title.y = element_blank(), legend.position = "None", plot.background = element_blank(), 
      axis.text.x=element_text(angle=45,vjust=1,hjust=1), axis.title.x = element_text(vjust=-0.2)) + labs(x="Collection Date") 

## C: Find co-occurence of D399N with other N mutations
full_n_muts_list <- list()
D399N_upset <- data.frame()

# Get rid of parentheses
D399N_geo_final[] <- lapply(D399N_geo_final, gsub, pattern=')',replacement='')
D399N_geo_final[] <- lapply(D399N_geo_final, gsub, pattern='\\(',replacement='')

# Go through and format GISAID N mutations into their own column.
for(i in 1:nrow(D399N_geo_final)) {
  row <- D399N_geo_final[i,]
  muts_list <- unlist(strsplit(row$`AA Substitutions`,","))
  n_muts_list <- (sort(muts_list[grepl("^N_.*.",muts_list)]))
  for(n_mut in n_muts_list) {
    row[, n_mut] <- 1
  }
  D399N_upset <- rbind.fill(D399N_upset, row)
}
# Fill NAs with 0
D399N_upset[is.na(D399N_upset)] <- 0

# These mutations have <10 sequences, get thrown into other N mutations
mutations_to_combine <- D399N_upset %>% select_if(negate(function(col) is.numeric(col) && sum(col) >10))
mutations_to_combine <- mutations_to_combine[c(27:ncol(mutations_to_combine))]
D399N_upset <- D399N_upset %>% select_if(negate(function(col) is.numeric(col) && sum(col) <=10))
D399N_upset$Other_N_mutations <- rowSums(mutations_to_combine)
D399N_upset$Other_N_mutations[D399N_upset$Other_N_mutations>=1] <- 1

# D399N_upset <- D399N_upset %>% unite("Other_N_mutations", c("N_I292T","N_P13T","N_S193N","N_I130V","N_A220V","N_P365S","N_D3Y","N_S197L",
#                                                                       "N_A376T","N_M234I","N_P80L","N_R385K","N_N345S","N_K388R","N_Q389K","N_Q389L",
#                                                                       "N_Q406stop","N_H59R", "N_P364L","N_A152V","N_I94V","N_T391I","N_S183P"))


# Formatting for upset plot
unique_n_muts <- unlist(colnames(D399N_upset))
unique_n_muts <- c(unique_n_muts[grepl("^N_.*.",colnames(D399N_upset))],"Other_N_mutations")

# saved separately and combined later
upsetplot <- upset(D399N_upset, sets = unique_n_muts, order.by = "freq", set_size.show = TRUE, set_size.scale_max = 250,
      mainbar.y.label = "Frequency of N Mutation Intersections", sets.x.label = "# Sequences with Mutation",
      text.scale = c(1.3,1,1.3,1,1.3,1), mb.ratio = c(0.6,0.4),
      queries = list(list(query = intersects, params = list("N_D399N"), color = "darkorange", active = T)))

#combined <- plot_grid(worldplot, timemuts, labels=c("A","B"), ncol=1, rel_widths = c(1,0.5),rel_heights=c(0.6,0.7))

ggsave("doublemut_gisaid_locations_onlyD399N_0418.pdf",combined,units="in",width=8.5,height=16)
