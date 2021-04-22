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

setwd("/Users/gerbix/Documents/Michelle/swift_testing/quidel_discordant_fastas/")

## Generated from GISAID master list of sequences + variants, pulled on 4/18/21
# doublemuts_gisaid <- read_tsv("bothmuts.tsv")
# doublemuts_gisaid <- read_tsv("only_T205I.tsv")
doublemuts_gisaid <- read_tsv("onlyD399n_0418.tsv")

# Format GISAID-specified locations into usable format
doublemuts_gisaid <- separate(doublemuts_gisaid, Location, into=c("continent","region","subregion"), sep=" / ") %>%
                   mutate(region = recode(region,"United Kingdom" = "UK"), subregion = recode(subregion,"England" = "Great Britain"))

# Get world map
world_map <- map_data("world") %>% filter(region!="Antarctica")
# Find mean latitude and longitudes for each subregion
cleaned_world_map <- world_map %>% group_by(region, subregion) %>% dplyr::summarise(mean_lat = mean(lat), mean_long = mean(long))
# Do US and UK separately 
cleaned_world_map_sub <- cleaned_world_map %>% filter(region=="UK")
cleaned_us_map <- map_data("state") %>% group_by(region) %>% dplyr::summarise(mean_lat = mean(lat), mean_long = mean(long))
cleaned_us_map$region <- str_to_title(cleaned_us_map$region)
colnames(cleaned_us_map) <- c("subregion","mean_lat","mean_long")
cleaned_us_map$region <- "USA"
cleaned_world_map <- cleaned_world_map %>% filter(region!="USA" & region!="UK") %>% dplyr::summarise(mean_lat = mean(mean_lat), mean_long = mean(mean_long))
doublemuts_gisaid_uk <- doublemuts_gisaid %>% filter(region=="UK")
doublemuts_gisaid_us <- doublemuts_gisaid %>% filter(region=="USA")
double_muts_gisaid <- doublemuts_gisaid %>% filter(region!="USA" & region!="UK")
doublemuts_geo_uk <- inner_join(cleaned_world_map_sub, doublemuts_gisaid_uk, by=c("region","subregion"))
doublemuts_geo_us <- inner_join(cleaned_us_map, doublemuts_gisaid_us, by=c("region","subregion"))
doublemuts_geo <- inner_join(cleaned_world_map, doublemuts_gisaid, by=c("region"))
doublemuts_geo_final <- rbind(doublemuts_geo, doublemuts_geo_uk,doublemuts_geo_us)


# Make a color scale for clades
myColors <- c(brewer.pal(7, "Set2"),brewer.pal(8,"Accent"),brewer.pal(7, "Dark2"),brewer.pal(8,"Set3"))
names(myColors) <- levels(doublemuts_geo_final$`Pango lineage`)
colScale <- scale_color_manual(name = , values=myColors)

## A: Plot geographic locations of GISAID consensuses with D399N
doublemuts_geo_final_loc <- doublemuts_geo_final %>% group_by(region,subregion,mean_lat,mean_long,`Pango lineage`) %>% 
                            dplyr::summarise(num_samples = n_distinct(`Accession ID`))
worldplot <- ggplot(data = doublemuts_geo_final_loc, 
                    aes(x=mean_long, y=mean_lat, color=`Pango lineage`, size=num_samples)) + 
              #scale_size_continuous(limits=c(1,13),breaks=c(1,5,9,13)) + 
              scale_size_continuous(limits=c(1,30),breaks=c(1,10,20,30)) + 
              #scale_color_brewer(palette = "Dark2", direction=-1) + 
              #scale_color_brewer(palette = "Set2",direction=-1) + 
              colScale + 
            geom_polygon(data = world_map, aes(x=long, y= lat, group= group),fill="lightgray",colour="white", size=0.1) + 
            coord_fixed(1.3)  + 
            theme_void() + scale_x_continuous(limits=c(-170,175)) + geom_point(alpha=0.9) + 
            theme(legend.box="vertical",legend.position="bottom") + labs(color="Lineage",size="# Samples")
            
# Format collection date
doublemuts_geo_final$`Collection date` <- as.POSIXct.Date(as.Date(doublemuts_geo_final$`Collection date`))

#doublemuts_geo_final <- doublemuts_geo_final %>% filter(`Collection date` < as.POSIXct("2020-06-15 01:00:00", tz="UTC"))
doublemuts_geo_final$subregion[doublemuts_geo_final$region=="Luxembourg"] <- "Luxembourg"

## B: Plot GISAID consensuses with D399N over time
timemuts <- ggplot(data=doublemuts_geo_final, aes(x=`Collection date`,y=subregion, color=`Pango lineage`)) + #color=`Pango lineage` 
        #scale_color_brewer(palette="Dark2",direction=-1) +
        colScale + 
        geom_point() + facet_grid(rows = vars(region), scales="free",switch="y",space="free_y") +
          scale_x_datetime(date_breaks = "1 month") + theme_clean() + 
      theme(strip.placement="outside",strip.text.y.left = element_text(angle=0), strip.background.y = element_rect(fill="gray90"),
            axis.title.y = element_blank(), legend.position = "None", plot.background = element_blank(), 
      axis.text.x=element_text(angle=45,vjust=1,hjust=1), axis.title.x = element_text(vjust=-0.2)) + labs(x="Collection Date") 

## C: Find co-occurence of D399N with other N mutations
full_n_muts_list <- list()
doublemuts_upset <- data.frame()
# Go through and format GISAID N mutations into their own column
for(i in 1:nrow(doublemuts_geo_final)) {
  row <- doublemuts_geo_final[i,]
  muts_list <- unlist(strsplit(row$`AA Substitutions`,","))
  n_muts_list <- (sort(muts_list[grepl("^N_.*.",muts_list)]))
  for(n_mut in n_muts_list) {
    row[, n_mut] <- 1
  }
  doublemuts_upset <- rbind.fill(doublemuts_upset, row)
}
# Fill NAs with 0
doublemuts_upset[is.na(doublemuts_upset)] <- 0

# These mutations have <10 sequences, get thrown into other N mutations
doublemuts_upset <- doublemuts_upset %>% unite("Other_N_mutations", c("N_I292T","N_P13T","N_S193N","N_I130V","N_A220V","N_P365S","N_D3Y","N_S197L",
                                                                      "N_A376T","N_M234I","N_P80L","N_R385K","N_N345S","N_K388R","N_Q389K","N_Q389L",
                                                                      "N_Q406stop","N_H59R", "N_P364L","N_A152V","N_I94V","N_T391I","N_S183P"))
# Formatting for upset plot
doublemuts_upset <- doublemuts_upset %>% mutate(Other_N_mutations = ifelse(test = (Other_N_mutations == "0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0"),
                                                            yes = 0,
                                                            no = 1))
unique_n_muts <- unlist(colnames(doublemuts_upset))
unique_n_muts <- c(unique_n_muts[grepl("^N_.*.",colnames(doublemuts_upset))],"Other_N_mutations")
doublemuts_upset$Other_N_mutations <- as.numeric(as.character(doublemuts_upset$Other_N_mutations))

# saved separately and combined later
upsetplot <- upset(doublemuts_upset, sets = unique_n_muts, order.by = "freq", set_size.show = TRUE, set_size.scale_max = 250,
      mainbar.y.label = "Frequency of N Mutation Intersections", sets.x.label = "# Sequences with Mutation",
      text.scale = c(1.3,1,1.3,1,1.3,1), mb.ratio = c(0.6,0.4),
      queries = list(list(query = intersects, params = list("N_D399N"), color = "darkorange", active = T)))

combined <- plot_grid(worldplot, timemuts, labels=c("A","B"), ncol=1, rel_widths = c(1,0.5),rel_heights=c(0.6,0.7))

ggsave("doublemut_gisaid_locations_onlyD399N_0418.pdf",combined,units="in",width=8.5,height=16)
