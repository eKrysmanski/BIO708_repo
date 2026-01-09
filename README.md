New to this so not sure what the conventions are


#DIRECTORIES
data > Includes proteomics data from two independent experimental replicates, already normalized
       Data was collected by Dr. Phil Carella, formerly of the Cameron Lab (my lab)
       Published in Plant Physiology April 2016 (DOI: 10.1104/pp.16.00269)

scripts > Includes R scripts

Phloem_Proteome files:
   Accession -- Accession number for Arabidopsis protein
   Peptides used for quantitation -- Number of unique peptides used for quantitation
   Confidence -- Confidence score given by Progenesis QI software
   ANOVA (p) -- ANOVA test results between the groups (mock, vir, avr)
   Description -- I believe TAIR descriptions for the protein
   VLC.../DES... -- normalized abundance values for biological replicates 
         ends_with("MockPEX#") : protein abundance from plants that were infiltrated with mock inoculum (MgCl2)
         ends_with("VirPEX#") : protein abundance form plants that were infiltrated with Pst DC3000
         ends_with("AvrPEX#") : protein abundance from plants that were infiltrated with Pst expressing avrRpt2

#Context for the dataset:

The purpose of the proteomics was to identify proteins in the phloem exudates (i.e., phloem sap) of plants 
induced for a form of plant immunity  called Systemic Acquired Resistance (SAR). SAR is a immune response in 
plants where following infection of a leaf, long-distance signals are generated/activated to move systemically
via the phloem (vascular tissue) to distant leaves, where they signal the distant uninfected leaves to become
primed for defense, prior to interaction with any pathogens. This results in a more rapid and robust defense reponse
to subsequent pathogens. The phloem sap from mock-, vir-, or avr-induced plants was collected and subjected to liquid 
chromatography tandem mass-spectrometry to identify differentially abundant proteins in the phloem sap. 

#Biological Questions:

Initially I took this dataset as we had intended to perform additional LC MS/MS in our lab using different methods 
of SAR induction which did not rely on a live pathogen. We planned to induce SAR by expressing a protein that 
is recognized by plant immune receptors under and estrogen-inducible promoter, or infiltrating with a bacteria-derived 
peptide that is recognized by other plant immune receotors, to induce two distinct but overlapping immune pathways. We
then wanted to collect exudates from these plants for comparative proteomics. 

The question we wanted to ask was "Following SAR induction, what proteins enter the phloem in Arabidopsis"
More specifically we were asking "What proteins are more abundant in the phloem of plants induced for PAMP-triggered immunity" 
                  and..          "What proteins are more abundant in the phloem of plants induced for effector-triggered immunity"
Since both PAMP- and Effector triggered immunity are sufficient to induce SAR, and Pst avrRpt2 induces both immune pathways
simultaneously, we reasoned that common proteins found in the phloem proteomes of each induction method would comprise of 
important SAR long-distance signaling components. Moreover, since we are using a pathogen-free system, we could avoid 
any effects that live pathogens and their secreted effector proteins could have on the proteinaceous contents of the 
phloem during SAR signaling. 

Due to renovations and the building of the greenhouse, we could not perform the large-scale phloem sap collection
we intended since plants are incredibly sensitive, and stress can induce immune responses that would contaminate our 
results. Instead of abandoning the idea, I decided it might still be worth doing comparisons of our dataset to 
other publicly available SAR proteome datasets to see if any interesting proteins are revealed. I have found two
intercellular/apoplast SAR proteome datasets which I might like to compare our phloem proteome to. 

#Note: this might not be the dataset I continue with for the course, but it is one I could easily access and was familiar to me. 

Our lab also has an RNA-seq dataset which has not been well analyzed, and I may opt to use that to perform some analyses
once I have figured out a line of questioning pertinent to the data. I might want to compare the PTI dataset to published SAR
RNA seq data, RNA-seq data from plants induced in various ways. 
