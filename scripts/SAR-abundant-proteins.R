library(dplyr)
library(tidyselect)
library(tidyr)

#Set WD
setwd("C:/Users/eKrys/Desktop/BIO708/Assignment_1")

#Reading in a subset of the data from Carella's phloem proteome rep #1
phloem_1 <- read.csv("phloem-proteome-1.csv", header = FALSE)
colnames(phloem_1) <-(phloem_1[3,])                #Setting names of headers
phloem_1 <- phloem_1[4:nrow(phloem_1), 1:20]  %>%  #Reformatting dataframe
  separate_rows(Accession, sep = ";")   #For now this is a quick and dirty fix which might confuse me later               
                                        #There are some rows with two accession numbers but its the same protein
                                        #For now I will just duplicate the rows with each of the protein accession #s
#Try something simple to begin with; let's figure out the proteins that are more
## abundant in the vir induced relative to mock induced plants. 

#Start by preparing a new dataframe with summarized data
#Begin by checking what my dataframe looks like, make sure vectors are correct type
str(phloem_1)

#Need to fix the abundance columns since they are being read as chr...
#I want to calculate the ratios and t-tests myself so I won't worry about those

phloem_1 <- phloem_1 %>%
  mutate(across(matches("-MockPEX\\d+$|-VirPEX\\d+$|-AvrPEX\\d+$"), as.numeric))

#This next part is making a summary dataframe; it calculates the mean FC, and 
# performs the t-tests on the data to look for more abundant proteins. 

#Note: the t-test here assumes unequal variance, and is a single-tailed t-test, 
# since I'm specifically interested inif protiens are more abundant rather than 
# just differentially abundant in the SAR vs. mock. If I want to also check for 
# less abundant I would need to adjust this to a two-tailed test because I think it 
# is an error to perform  multiple statistical tests on the same data... I think...
# Perhaps I need to consult someone who knows this stuff better. 

phloem_1_sum <- phloem_1 %>%
  rowwise() %>%
  mutate(
    mean_mock = mean(c_across(matches("-MockPEX\\d+$")), na.rm = TRUE), 
    mean_vir = mean(c_across(matches("-VirPEX\\d+$")), na.rm = TRUE), 
    mean_avr = mean(c_across(matches("-AvrPEX\\d+$")), na.rm = TRUE),
    FC_vir = mean(c_across(matches("-VirPEX\\d+$")))/mean(c_across(matches("-MockPEX\\d+$"))),
    FC_avr = mean(c_across(matches("-AvrPEX\\d+$")))/mean(c_across(matches("-MockPEX\\d+$"))),
    ttest_mock_vir = t.test(c_across(matches("-MockPEX\\d+$")), 
                            c_across(matches("-VirPEX\\d+$")), 
                            alternative = "less")$p.value, 
    ttest_mock_avr = t.test(c_across(matches("-MockPEX\\d+$")), 
                            c_across(matches("-AvrPEX\\d+$")), 
                            alternative = "less")$p.value) %>%
  ungroup() %>%
  select(Accession, Description, 
         mean_mock, mean_vir, mean_avr, 
         FC_vir, FC_avr, 
         ttest_mock_vir, ttest_mock_avr, 
         `Peptides used for quantitation`,`Confidence score`)

#Note: I used AI to help me figure out the correct regex to select the columns, and 
# also to troubleshoot an issue where I was struggling to get rowise() to work when 
# trying to calculate the mean of the abundance measurements for each protein

#Now let's try to extract a subset of this dataframe that is more abundant in
# virulent-induced plants relative to mock-induced plants. This should just be 
# a simple filter function

p1_inc_vir <- phloem_1_sum %>%
  filter(mean_vir > mean_mock, 
         ttest_mock_vir < 0.05, 
         `Peptides used for quantitation` >=2) #seems to be convention?

length(unique(p1_inc_vir$Description))   #Because I did that cheaty thing with duplicating rows
                                         # i'm looking at the number of unique descriptions

#Note to self:: Phil reported 42 for this dataset, he may have done some manual 
# filtering. I don't think he had access to R when he did this, and I believe he 
# did it all manually in excel rather than using a script or something. Maybe he 
# missed some things in his analysis? It's not clear how he went about filtering
# the data to get to the conclusion he did. Maybe this will reveal some interesting
# proteins...


#Let's do the Same thing for the avirulent induced plants

p1_inc_avr <- phloem_1_sum %>%
  filter(mean_avr > mean_mock,
         ttest_mock_avr < 0.05, 
         `Peptides used for quantitation` >= 2)

length(unique(p1_inc_avr$Description))

#This is reporting 67, but phil reported 71.. so what and where are those 4 missing?
# I probably need to put it into a simple spreadsheet and compare to quickly figure
# out where the differences are cause I don't want to manually ctrl + f each accession
# until I find one that is missing in his report. 

#Let's figure out which proteins are more abundant in both vir and avr

p1_inc_vir_avr <- phloem_1_sum %>%
  filter(mean_vir > mean_mock, 
         mean_avr > mean_mock, 
         ttest_mock_vir < 0.05, 
         ttest_mock_avr < 0.05, 
         `Peptides used for quantitation` >= 2)


length(p1_inc_vir_avr$Description)

#Phil only reported 30... Need to figure out how our analyses are different. 

#=========================CARELLA PHLOEM PROT. 2 ==============================#

#Read/format in a subset of the data from Carella's phloem protein rep #2
#For whatever reason the excel sheets are not in the same formats, so slightly
# different formatting process

phloem_2 <- read.csv("phloem-proteome-2.csv", header = FALSE)
colnames(phloem_2) <- (phloem_2[3,])               #Collecting names of actual headers
phloem_2 <- phloem_2[4:nrow(phloem_2),] %>%        #Reformatting dataframe 
  separate_rows(Accession, sep = ";")              #Seperating the columns with 
#multiple possible protein calls


#Prepare a new dataframe with the data summarized, begin by specifying all these
# values are numeric so the following code doesnt get confused

phloem_2 <- phloem_2 %>%
  mutate(
    across(matches("-Mock-PEX\\d+$|-Vir-PEX\\d+$|-Avr-PEX\\d+$"), as.numeric))


#Create the summary dataframe

phloem_2_sum <- phloem_2 %>%
  rowwise() %>%
  mutate(
    mean_mock = mean(c_across(matches("-Mock-PEX\\d+$")), na.rm = TRUE), 
    mean_vir = mean(c_across(matches("-Vir-PEX\\d+$")), na.rm = TRUE), 
    mean_avr = mean(c_across(matches("-Avr-PEX\\d+$")), na.rm = TRUE),
    FC_vir = mean(c_across(matches("-Vir-PEX\\d+$")))/mean(c_across(matches("-Mock-PEX\\d+$"))),
    FC_avr = mean(c_across(matches("-Avr-PEX\\d+$")))/mean(c_across(matches("-Mock-PEX\\d+$"))),
    ttest_mock_vir = t.test(c_across(matches("-Mock-PEX\\d+$")), 
                            c_across(matches("-Vir-PEX\\d+$")), 
                            alternative = "less")$p.value, 
    ttest_mock_avr = t.test(c_across(matches("-Mock-PEX\\d+$")), 
                            c_across(matches("-Avr-PEX\\d+$")), 
                            alternative = "less")$p.value) %>%
  ungroup() %>%
  select(Accession, Description, 
         mean_mock, mean_vir, mean_avr, 
         FC_vir, FC_avr, 
         ttest_mock_vir, ttest_mock_avr, 
         `Unique peptides`,`Confidence score`)

#I want to know where these proteins are usually, and what they do. 
# To figure out where they are and what they are doing, I'm going to be using 
# the TAIR functional descriptions and a dataset called SUBA like I did above

#Determine what proteins are more abundant in vir vs. mock
p2_inc_vir <- phloem_2_sum %>%
  filter(mean_vir > mean_mock, 
         ttest_mock_vir < 0.05, 
         `Unique peptides` >1)

length(unique(p2_inc_vir$Description))

#Phil reported 97 for this dataset;

#Determine what proteins are more abundant in the avr vs. mock
p2_inc_avr <- phloem_2_sum %>%
  filter(mean_avr > mean_mock,
         ttest_mock_avr < 0.05, 
         `Unique peptides` >= 2)

length(unique(p2_inc_avr$Description))

#Phil reported 162 for this dataset -- interesting how this time he found more...

#Figure out which proteins are more abundant in both vir and avr
p2_inc_vir_avr <- phloem_2_sum %>%
  filter(mean_vir > mean_mock, 
         mean_avr > mean_mock, 
         ttest_mock_vir < 0.05, 
         ttest_mock_avr < 0.05, 
         `Unique peptides` >1)


length(p2_inc_vir_avr$Description)

#Phil reported 94

#==============Time to compare Phloem 1 and Phloem 2==========================#

#Inner_join will compare the two df by the Accession row and keep rows that had
# matching values (i.e., keep proteins abundant in both proteomes)

SAR_up <- inner_join(p1_inc_vir_avr, 
                     p2_inc_vir_avr, 
                     by = "Accession")
SAR_up <- SAR_up %>%
  filter(!str_detect(Accession, "\\.2(\\b|;)")) %>% 
  select(Accession, Description.x,
         FC_vir.x, FC_vir.y,
         FC_avr.x, FC_avr.y, 
         Description.x)

length(unique(SAR_up$Description.x))

#Phil reported 18 herex