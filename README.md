# LipidCruncher_Demo

The mission of Farese/Walther lab at Harvard/MSKCC is to understand cell metabolism: how our bodies manage the usage and storage of energy. 
They specifically focus on understanding the function of lipid molecules inside cells as lipids are where energy is stored. 
A typical study is as follows: the lab biologists come up with a hypothesis (driven by deep expertise and intuition built over many years), 
then they design/run an experiment to test the hypothsis. 
They usually have a control group where nothing is changed and one or more experimental groups where something is changed 
(i.g. a certain gene is knocked out or is over expressed). They prepare multiple samples (replicates) for each group and run them though mass spectrometer. 
One of the things that they are interested in is whether there is a significant difference in the abundance of lipid species between the control group 
and the experimental group(s).

In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. 
This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak reperesnts 
the relative abundance of the corresponding lipid species. There are two pieces of software that turn this spectrum into a lipidomics dataset: 
LipidSearch and LipidXplorer. 
I built LipidCruncher which is a web app that allows the user to perform lipidomics analysis on the LipidSearch and LipidXplorer datasets 
[(Link to LipidCruncher Youtube Demo)](https://www.youtube.com/watch?v=KC4eLuwYw3A).
