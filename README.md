# Simulation
 tree_surface.3<-OHA_Tree_phylo_NOOUTGROUPS
elevation_data.3<-read.csv("OHA_Elevation_2.csv", row.names = 1)
elevation_data.3<-elevation_data.3[!is.na(as.numeric(elevation_data.3$TmeanF,elevation_data.3$Prec, elevation_data.3$Elevation)), ]
 dat_surface.3<-log(elevation_data.3)
 tree_surface.3<-nameNodes(tree_surface.3)
 olist_surface.3<-convertTreeData(tree_surface.3,dat_surface.3)
 otree_surface.3<-olist_surface.3[[1]]
 odata_surface.3<-olist_surface.3[[2]]
 fwd.3<-surfaceForward(otree_surface.3, odata_surface.3, aic_threshold = 0, exclude = 0, verbose = FALSE, plotaic = FALSE)
 k.3<-length(fwd.3)
  bwd.3<-surfaceBackward(otree_surface.3, odata_surface.3, starting_model = fwd.3[[k.3]], aic_threshold = 0, only_best = FALSE, verbose = FALSE, plotaic = FALSE)
 bsum.3<-surfaceSummary(bwd.3)
 kk.3<-length(bwd.3)
 # To see the Surface Tree
 surfaceTreePlot(tree_surface.3, bwd.3[[kk.3]], labelshifts =T, show.tip=TRUE)
 # Code to run the models
 bm.3<-startingModel(otree_surface.3,odata_surface.3,brownian=TRUE)
 ou1.3<-startingModel(otree_surface.3,odata_surface.3)
 H12.3<-startingModel(otree_surface.3,odata_surface.3,shifts=c("26"="H1","13"="H1","5"="H2","19"="H2"))
 surfaceAICPlot(fwd.3, bwd.3, traitplot = "aic")
 abline(h=bm.3[[1]]$aic,lty="longdash")
 abline(h=H12.3[[1]]$aic,lty="longdash")
 text(c(6,6),c(bm.3[[1]]$aic, ou1.3[[1]]$aic, H12.3[[1]]$aic)-2,c("BM","OU1","H12"),cex=0.5)
 legend ("bottomleft", legend = c("Elevation", "TmeanF", "Prec"), title="Elevation Data", col = c("red", "blue","green"), lty=1, cex=0.5)
set.seed(100)
newsim.Hansen <-surfaceSimulate(tree_surface.3, type="hansen-fit", hansenfit=fwd.3[[k.3]]$fit,shifts=fwd.3[[k.3]]$savedshifts, sample_optima=TRUE, no_nested = TRUE, optima_type="even")
newsim.BM<-surfaceSimulate(tree_surface.3, type="BM",shifts=fwd.3[[k.3]]$savedshifts, sample_optima=TRUE, no_nested = TRUE, param = 0)
 par(mfrow=c(1,2),mai=c(0.8,0.8,0.2,0.2))
 surfaceTraitPlot(newsim.Hansen$data, newsim.Hansen, whattraits = c(1,2), convcol = TRUE)
surfaceTraitPlot(newsim.Hansen$data, newsim.Hansen, whattraits = c(3,2), convcol = TRUE)
surfaceTraitPlot(newsim.Hansen$data, newsim.Hansen, whattraits = c(3,1), convcol = TRUE)
#We can then run SURFACE on the simulated data set, here doing the entire analysis in one
#step using runSurface (this should take less than a minute). We can then use surfaceSummary
#to extract the results, and compare the number of regime shifts (k) and the extent of convergence
#(deltak or c) to what we saw in the `real' data set. In this case, one instance of convergence is
#recovered in the simulated data set where two regimes were relatively close to one another in trait
#space; as we know that the regimes were nonetheless distinct from one another in the generating
#model, this represents `incidental' convergence.
newout.Hansen<-runSurface(tree_surface.3, newsim.Hansen$dat)
newout.BM<-runSurface(tree_surface.3, newsim.BM$dat)
newsum.Hansen<-surfaceSummary(newout.Hansen$bwd)
newsum.BM<-surfaceSummary(newout.BM$bwd)
