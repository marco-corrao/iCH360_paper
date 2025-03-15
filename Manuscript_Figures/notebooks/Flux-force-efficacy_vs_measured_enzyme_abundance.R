library(plyr,include.only = "mapvalues")
library(dplyr)
library(ggplot2)
library(ggpubr)

pta_data=
  read.csv('../../Analysis/PTA/out/pta_fluxes.csv') 
reaction_enzyme_map=
  read.csv('../../Enzyme_constraints/EC_model_building_pipeline/reaction_enzyme_kcat_mapping/parsed/EC_table_w_manual_curation.csv')
enzyme_abundance_data=
  read.csv("../../Enzyme_Constraints/enzyme_abundance_estimation/schmidt_2016_copies_per_cell_abundance_mapped_NNLS.csv") 

glucose_uptake_flux=pta_data %>% 
  filter(id=="GLCptspp") %>% 
  .$v

RT=8.31446261815324e-3 * 310.15
ffe=function(drg){(exp(-drg/RT)-1)/(exp(-drg/RT)+1)}
data_processed=
  pta_data %>%
  mutate(dir=ifelse(v>=0,"fw","bw")) %>% 
  mutate(directional_id=paste(id,dir,sep='_')) %>% 
  mutate(enzyme=mapvalues(x=directional_id,
                          from=reaction_enzyme_map$directional_reaction_id,
                          to=reaction_enzyme_map$enzyme,
                          warn_missing=F)
         ) %>% 
  mutate(enzyme_abundance=mapvalues(x=enzyme,
                                    from=enzyme_abundance_data$enzyme,
                                    to=enzyme_abundance_data$glucose_copies_per_cell,
                                    warn_missing=F) %>% as.numeric()) %>% 
  mutate(relative_abs_fluxes=abs(v/glucose_uptake_flux))%>% 
  mutate(driving_force=-drg*sign(v)) %>% 
  mutate(flux_force_efficacy=ffe(-driving_force)) %>% 
  mutate(low_flux_force=(flux_force_efficacy<0.5)) %>% 
  mutate(group=if_else(low_flux_force,"<0.5",">0.5"))


      
data_plot=
  data_processed%>% 
  filter(relative_abs_fluxes>0.025) %>% 
  group_by(group) %>% 
  mutate(median=median(enzyme_abundance,na.rm=T))

number_of_observations=
  data_plot %>% 
  group_by(group) %>% 
  summarise(count=n(),
            label_y_pos=median(enzyme_abundance,na.rm=T)*1.5,
            )%>%
  mutate(text=paste("n= ",count))



#======================

summary_df=
    data_plot %>% 
  group_by(group) %>% 
  summarise(count=n(),
            median=median(enzyme_abundance,na.rm=T),
            label_x_pos=12000
  ) %>% 
  mutate(text=paste("median=",round(median,0)),
         label_y_pos=c(0.15,0.25))


  data_plot %>%
  ggplot(aes(x=enzyme_abundance,color=group))+
  stat_ecdf(geom="step",size=1)+
  scale_x_continuous(trans = "log10")+
  scale_color_manual(values=c("#b40000ff","#4C72B0"))+
  geom_segment(aes(x=data_plot %>% filter(group=="<0.5") %>% .$enzyme_abundance %>% median(na.rm=T),
                   y=0,
                   xend=data_plot %>% filter(group=="<0.5") %>% .$enzyme_abundance %>% median(na.rm=T),
                   yend=0.5),
               linetype="dotted",
               linewidth=1,
               color="#b40000ff"
  )+
  geom_segment(aes(x=data_plot %>% filter(group==">0.5") %>% .$enzyme_abundance %>% median(na.rm=T),
                   y=0,
                   xend=data_plot %>% filter(group==">0.5") %>% .$enzyme_abundance %>% median(na.rm=T),
                   yend=0.5),
               linetype="dotted",
               linewidth=1,
               color="#4C72B0"
               )+
  labs(x="Measured enzyme abundance (copies per cell)",
       y="Empirical CDF",
       color="Flux Force Efficacy")+
  geom_text(data=summary_df,
            aes(x=label_x_pos,y=label_y_pos,label=text),
            show.legend = F)+
  ylim(0.00,1.01)+
    coord_cartesian(expand=F)+
  theme_bw()->
  cdf_plot

#Run a p-value test
x=data_plot %>% filter(group=="<0.5") %>% .$enzyme_abundance
y=data_plot %>% filter(group==">0.5") %>% .$enzyme_abundance
test=wilcox.test(x,y,alternative="two.sided")

#Add pvalue onto plot
cdf_plot=
  cdf_plot+
  geom_text(aes(x=400,y=0.5,label="p< 0.001"),color="black")
cdf_plot
#ggsave("../figures/PTA_FFE_vs_enzyme_abundance_cdf_comparison.pdf",cdf_plot,width=8,height=6)


























