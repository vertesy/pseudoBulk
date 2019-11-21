# ------------------------------------------------------------------------------------------
## Seurat3.Write.Out.CBCs.for.subset-bam.R
# ------------------------------------------------------------------------------------------
# for https://github.com/10XGenomics/subset-bam

"I did it in /Users/abel.vertesy/GitHub/OrgRepCard/ORC.elements/Write.barcodes.per.cluster.4.split-bam.R"

#  ------------------------------------------------------------
get.ls.of.CBCs <- function(scobj = combined.obj, ident = 'integrated_snn_res.0.3', 
                           plotit=T, trim.libName=F) { # 
  Idents(scobj) <- ident
  id_x = Idents(scobj)
  dsets = unique(stringr::str_split_fixed(names(id_x), pattern = '_', n = 2)[,2])
  iprint("Returns a list of lists (",l(dsets),"libraries [",dsets,"] /",
         l(levels(id_x)),"clusters [",levels(id_x),"] )")
  ls_CBCs = list.fromNames(levels(id_x))
  print("Cluster:")
  AllClusters = levels(id_x)
  cl=1
  for (cl in 1:l(AllClusters)) { 
    print(cl)
    cells_x = WhichCells(combined.obj, idents = id_x[cl])
    cells_perLib = stringr::str_split_fixed(cells_x, pattern = '_', n = 2)
    cells_perLib = cbind(cells_perLib,scobj$orig.ident[cells_x])
    ls_CBCs[[cl]] <- ls_cells_clX_perLib <- split(x = cells_perLib, f = cells_perLib[,3])
  }
  revlist = reverse.list.hierarchy(ls_CBCs)
  
  if (!isFALSE(trim.libName)) { # remove .WT and .TSC2
    names(revlist) = stringr::str_split_fixed(names(revlist), pattern = "\\.", n=2)[,1]
  }
  if (plotit) {
    for (i in 1:l(revlist)) {
      ClusterSizes = unlapply(revlist[[i]], l)
      wpie(ClusterSizes, savefile = F, plotname = paste('ClusterSizes, cl.', i))
    }
  }
  return(revlist)
}
# ls.of.CBCs = get.ls.of.CBCs(trim.libName = T); names(ls.of.CBCs)

#  ------------------------------------------------------------

write.out.CBCs.per.cl <- function(ls_CBCs = ls.of.CBCs, add.suffix="-1",
                                  ident = 'integrated_snn_res.0.3', 
                                  openOutDir=T, writeMeta=T) { # take the output of get.ls.of.CBCs() as input, and write out as csv
  (depth = l(ls_CBCs))
  (dsets = names(ls_CBCs))
  for (i in 1:l(dsets)) {
    outputDir = p0(OutDir,"CBCs/", dsets[i])
    dir.create(outputDir, recursive = T)
    inside.ls = ls_CBCs[[i]]
    for (j in 1:l(inside.ls)) {
      CBCs = inside.ls[[j]]
      if (!isFALSE(add.suffix)) CBCs = p0(CBCs,add.suffix)
      write.simple.vec(input_vec =  CBCs, ManualName =  p0(outputDir,"/Cl.", (j-1), ".csv"))
    }
  }
  if (writeMeta) {
    rez = str_split_fixed(ident,"snn_", n=2)[,2]
    (text_ = paste(names(ls_CBCs[[1]]), rez, idate()))
    write.simple.vec(input_vec =  text_, ManualName =
                       p0(OutDir,"CBCs/CBC.",rez,"__", idate(),".info"))
  }
  if (openOutDir) system(paste("open", outputDir))
}
# write.out.CBCs.per.cl()

#  ------------------------------------------------------------
write.out.config.file.template <- function(libname=names(ls.of.CBCs)[2], 
                                           bampath="/Volumes/abel/Data/bam.files.cellranger/",
                                           openOutDir=T) { # take the output of get.ls.of.CBCs() as input, and write out as csv
  outputDir2 = p0(OutDir,"CBCs/ConfigFiles/")
  
  ConfigFileTemplate = cbind(
    ppp("Cl", names(ls.of.CBCs[[1]]) ),
    rep(libname, length(ls.of.CBCs[[1]]) ),
    p0(bampath, names(ls.of.CBCs[[1]]) )
  )
  
  colnames(ConfigFileTemplate) = c("SAMPLEID","CONDITION", "SRC")
  write.csv(ConfigFileTemplate, file = p0(outputDir2,"/ConfigFileTemplate.",libname,".csv"),
            quote = F, row.names = F, col.names = NULL)
  if (openOutDir) system(paste("open", outputDir2))
}
write.out.config.file.template()

# Old deprecated, without names ------------------------------------------------------------

# old.get.ls.of.CBCs <- function(scobj = combined.obj, ident = 'integrated_snn_res.0.3', plotit=T) {
#   Idents(scobj) <- ident
#   id_x = Idents(scobj)
#   dsets = unique(stringr::str_split_fixed(names(id_x), pattern = '_', n = 2)[,2])
#   iprint("Retunts a list of lists (",l(dsets),"libraries [",dsets,"] /",
#          l(levels(id_x)),"clusters [",levels(id_x),"] )")
#   ls_CBCs = list.fromNames(levels(id_x))
#   print("Cluster:")
#   for (cl in 1:l(levels(id_x))) { 
#     print(cl)
#     cells_x = WhichCells(combined.obj, idents = id_x[cl])
#     cells_perLib = stringr::str_split_fixed(cells_x, pattern = '_', n = 2)
#     ls_CBCs[[cl]] <- ls_cells_clX_perLib <- split(x = cells_perLib, f = cells_perLib[,2])
#   }
#   revlist = reverse.list.hierarchy(ls_CBCs)
#   
#   if (plotit) {
#     for (i in 1:l(revlist)) {
#       ClusterSizes = unlapply(revlist[[i]], l)
#       wpie(ClusterSizes, savefile = F, plotname = paste('ClusterSizes, cl.', i))
#     }
#   }
#   return(revlist)
# }
