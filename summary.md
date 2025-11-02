TO DO

1\. length clusters - histogram -> AA composition for short, mid, long

2\. elm -> common motives for certain group (s, m, l | c-,n-,inner)

3\. order/disorder analysis

4\. secondary structure prediction

4\. phosphorylation + glycosylation



\- remove outliers + polyHis in C-,N- terminals

\- Some non-domain regions have high disorder predicted -> are these linkers too?

\- How to define linkers really?

\- homologní proteiny -> odlišnost exonů?

\- disorder X secundary structure relationship

-----------------------------------------------------------------------------



Interpro uses two complementary domain annotations:

'Representative domains' are consensus sequence-based annotations derived from databases such as Pfam, SMART etc.

'TED domains' are automatically predicted structure-based segments from AlphaFold models.

Representative domains reflect curated functional units while TED domains reflect compact structural units that may merge, split or overlap those curated regions.



**Representative domains**



Representative domains are derived from overlapping models from InterPro’s member databases (e.g. Pfam, SMART, TIGRFAMs, CDD, ProSite). InterPro merges these when they describe the same biological function and keeps one as a “representative” entry.

Each underlying member database uses its own annotation method, such as:

&nbsp;	Pfam / SMART – profile HMMs trained from multiple sequence alignments (MSA + HMM)



&nbsp;	ATH-Gene3D – structure-based HMMs derived from 3D domain definitions



&nbsp;	DD – PSSMs (position-specific scoring matrices)



&nbsp;	ProSite – regular expression–style sequence motifs or profiles



**TED domains**



These are structure-based domains automatically predicted from AlphaFold models using CATH's 3D domain segmentation algorithms. The aim is to partition the continuous 3d model into independently folding units; domains. These domains are usually shorted than the "Representative domains" from InterPro.



For searching domains and their coordinates, I have created a python script that queries InterPro for a given accession and returns

\- representative domains if representative\_only==True

\- all domains if representative\_only==False

-----------------------------------------------------------------------------

**Phase 1 results**



We have analysed 100 proteins from e. coli, homo sapiens and box taurus. Firstly, We have gathered representative domains and their coordinates. Then we removed all proteins with no domains found leaving 95 remaining.

We took all regions that are not part of any domain and divided them into C-terminal, N-terminal and inner linkers. Then we ran python analysis on domains + each linker group to look for length distribution and AA composition. The regions we considered as linkers showed higher S, G, P and E frequency which are all residues that either promote flexibility, break secondary structure, or are common in disordered regions. C-terminus linkers showed the highest G, S % composition and inner linkers the highest E %. That is in accordance with the previously established theory.

We also found out that these "linker" regions have lower D, N, I and L frequency.
