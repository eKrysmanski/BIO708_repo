library(dplyr)
library(tidyselect)
library(tidyr)

#Set WD to make sure I know where I am
setwd("C:/Users/eKrys/Desktop/GitRepo/BIO708_repo/data")

#Reading in a subset of the data from Carella's phloem proteome rep #1
phloem_1 <- read.csv("phloem-proteome-1.csv", header = FALSE)

#Formatting the dataframe
colnames(phloem_1) <-(phloem_1[3,])                       #Setting names of headers
phloem_1 <- phloem_1[4:nrow(phloem_1), 1:20]  %>%         #Reformatting dataframe
  separate_rows(Accession, sep = ";")   #For now this is a quick and dirty fix which might confuse me later               
                                        #There are some rows with two accession numbers but its the same protein
                                        #For now I will just duplicate the rows with each of the protein accession #s

#Working with the data:

#Try something simple to begin with; determine what  proteins that are more
# abundant in the vir induced relative to mock induced plants. To begin I think 
# it will be useful to summarize the data first so I have a clean data.frame to 
# work from. 

#Preparing a summarized dataframe
#Begin by checking what my dataframe looks like, make sure vectors are correct type
str(phloem_1)

#Need to fix the abundance columns since they are being read as chr...
#I want to calculate the ratios and t-tests myself so I won't worry about those
#Used AI to troubleshoot, was directed to use "Across". I think what was happening was
# mutate() expects mutate(new column = expression), and I'm giving it multiple
# columns?

phloem_1 <- phloem_1 %>% 
  mutate(across(starts_with("VLC"), as.numeric))

#Note: across() --> apply the same transformation to multiple columns, allowing you 
#      to use select() semantics inside in "data-masking" functions like summarise() 
#      and mutate().    [[  across(columns, function, names)]]

#Creating summary data.frame -- calculate the mean FC, and perform t-tests on 
# the data to look for more abundant proteins. 

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

#Note: I used AI to help me figure out the correct regex to select the columns, and it
# also directed me to use rowise(). Found statology page on how to do exactly what I wanted
# when googling (https://www.statology.org/dplyr-mean-for-multiple-columns/). c_across allows
# me to combine values from multiple columns. Basic syntax is like this: rowwise(function(c_across(columns)))
# to apply a function across multiple columns per row, like calculating the mean of 3 columns for each row
# in my case. 

#Now let's try to extract a subset of this dataframe that is more abundant in
# virulent-induced plants relative to mock-induced plants. This should just be 
# a simple filter function

p1_inc_vir <- phloem_1_sum %>%
  filter(mean_vir > mean_mock, 
         ttest_mock_vir < 0.05, 
         `Peptides used for quantitation` >=2) #seems to be convention to use >=2 unique peptides
                                               # for shotgun proteomics data. Maybe there is a way
                                               # to include high confidence single unique peptide 
                                               # proteins since I imagine this biases larger proteins

length(unique(p1_inc_vir$Description))   #Because I did that cheaty thing with duplicating rows
                                         # i'm using at the number of unique descriptions

#Note to self: Phil reported 42 for this dataset, he may have done some manual 
# filtering. I don't think he had access to R when he did this, and I believe he 
# did it all manually in excel rather than using a script or something. Maybe he 
# missed some things in his analysis? It's not clear how he went about filtering
# the data to get to the conclusion he did. Maybe this will reveal some interesting
# new proteins. Also note that it's not because I'm doing a single tailed t-test, 
# I checked with both methods. 


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
# values are numeric so the following code doesnt get confused. Regex tried to 
# mess with me here since DES also would include "DESCRIPTION", that was annoying

str(phloem_2)

phloem_2 <- phloem_2 %>%
  mutate(across(starts_with("DES7"), as.numeric))

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

#Combine the two inc_vir_avr dataframes by the accession number
#inner_join() checks if there is a matching value at the specified column ("by = ")
# then keeps it, if there is no match then it is removed. 

SAR_up <- inner_join(p1_inc_vir_avr, 
                     p2_inc_vir_avr, 
                     by = "Accession")

#Select the important columns with info I might want to check; descriptions will 
# be the same for both dataframes so keep only one, but the FC and other calculated
# values would be different. Keep FC for now so I can compare relative differences 
# in protien quantity, might want to include means for absolute as well, could mean 
# something if there is basically 0 in the mock plants and some quantity in the 
# induced plants. 

SAR_up <- SAR_up %>%
  select(Accession, Description.x,
         FC_vir.x, FC_vir.y,
         FC_avr.x, FC_avr.y, 
         Description.x)

length(unique(SAR_up$Description.x))

#Phil reported 18 here (yay?)

#Make a little table to summarize what I have so far for DAPs (diff. abundant proteins)
summ_DAPs <- data.frame(
  DAPs = c("Increased vir", "Increased avr", "Increased vir and avr"),
  Rep1 = c(length(unique(p1_inc_vir$Description)), 
           length(unique(p1_inc_avr$Description)), 
           length(unique(p1_inc_vir_avr$Description))), 
  Rep2 = c(length(unique(p2_inc_vir$Description)), 
           length(unique(p2_inc_avr$Description)), 
           length(unique(p2_inc_vir_avr$Description))), 
  Both = c("x", "x", length(unique(SAR_up$Description.x)))) %>% 
  print()
