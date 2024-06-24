library(readxl)
library(dplyr)

# Open abundance quantification and read correct headers ============
data_BW25113=read.csv('raw_data/Schimdt_et_al_2016_Table_S6_protein_copies_per_cell.csv',check.names=F) %>% 
  filter(!grepl(pattern='Isoform',x=Description))

data_MG1655=read.csv('raw_data/Schmidt_et_al_2016_MG1655_protein_copies_per_cell.csv',check.names=F) %>% 
  filter(!grepl(pattern='Isoform',x=Description))

#select columns of interest for strain MG1655==========================
data_MG1655=
  data_MG1655 %>% 
  select('Uniprot Accession',
         'abundance_glucose_copies_per_cell') %>% 
  rename(glucose_MG1655_pp_per_cell=abundance_glucose_copies_per_cell)

#Select columns of interest for strain BW25113=========================
columns_of_interest=c('Uniprot Accession','Gene','Molecular weight (Da)')
condition_of_interest=c('Glucose','Acetate','Fumarate','Glycerol','Pyruvate','Fructose','Succinate','Xylose', 
                        'Chemostat µ=0.5','Chemostat µ=0.12','Chemostat µ=0.20','Chemostat µ=0.35'
                        )
data_BW25113=data_BW25113 %>% select(all_of(c(columns_of_interest,condition_of_interest)))

#Now we are going to merge the two tables, leaving a NaN whenever a protein was not detecte
all_data=
  merge(data_BW25113,
        data_MG1655,
        by='Uniprot Accession',
        all=T)
columns_renames=c(uniprot_id='Uniprot Accession',
                  gene='Gene',
                  MW_Da='Molecular weight (Da)',
                  acetate_pp_per_cell='Acetate',
                  fumarate_pp_per_cell='Fumarate',
                  succinate_pp_per_cell='Succinate',
                  glycerol_pp_per_cell='Glycerol',
                  pyruvate_pp_per_cell='Pyruvate',
                  fructose_pp_per_cell='Fructose',
                  xylose_pp_per_cell='Xylose',
                  glucose_chemostat_mu_0_5_pp_per_cell='Chemostat µ=0.5',
                  glucose_chemostat_mu_0_12_pp_per_cell='Chemostat µ=0.12',
                  glucose_chemostat_mu_0_20_pp_per_cell='Chemostat µ=0.20',
                  glucose_chemostat_mu_0_35_pp_per_cell='Chemostat µ=0.35',
                  glucose_pp_per_cell='Glucose')
all_data =
  all_data %>% 
  rename(all_of(columns_renames))



# Now, we estimate the abundances in gDW per cell
avogadro=6.022e23
gDW_per_cell=2.8e-13

for(condition in colnames(all_data)[4:16]){
  all_data[,condition %>% gsub('pp_per_cell','g_gDW',.)]=(all_data[,condition]*all_data[,'MW_Da'])/(avogadro*gDW_per_cell)
}

all_data %>% 
  write.csv('parsed_data/Schmidt_et_al_2016_parsed.csv')
#Now parse total measured proteome per cell=========================================================
total_measured_proteome=data.frame(row.names = colnames(all_data)[17:28])

for (condition in total_measured_proteome %>% row.names()){
  total_measured_proteome[condition,'total_measured_g_gDW']=sum(all_data[,condition],na.rm=T)
}
total_measured_proteome$condition=total_measured_proteome %>% row.names()
total_measured_proteome %>% write.csv('parsed_data/total_measured_proteome_g_gDW.csv',row.names = F)
