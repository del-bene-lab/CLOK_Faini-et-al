In this subfolder we report some MatlabScript and exemplary data to generate analysis and plots analogous to some figures of the paper. 

- Plot_RatioDistribution_F.m : Plot distribution and location of cell and their dye ratio values ( Supp. Figure 15) .

- ratio_PositionANDclustering_Whalocamp.m : Treat location and fluorescence values from dual dye whalocamp fluorescent to cluster early and late born popualtion (Fig. 6a and Supp.Fig.20)

- The file WhaloCamp_traceClusteringANDselection_essential.zip contains both the essential variablales .mat file and a full code to selct and merge properly late  and early cell ROI, cluster them and make some basic analysis on calcium WhaloCamp traces (normalization, detrendingn, smoothing and rough manual selection of active cells).

- The file Trace_Manual_Analysis.zip , starting from a summary structure containg traces and their properties of EALY and LATE born cells extracte automatically (output of  WhaloCamp_traceClusteringANDselection_essential.zip) , It also provide an additiona GUI to visualize and screen rapdialy through all traces, adjusting their properites (such ACTIVE status , baseline and peaks). Then a global overview of the calcium activity across group of cell is plot. 




