parBIOMOD_Modeling <- function(data, models = c("GLM", "GBM", "GAM", "CTA", "ANN", 
                                                "SRE", "FDA", "MARS", "RF", "MAXENT.Phillips", "MAXENT.Phillips.2"), 
                               models.options = NULL, NbRunEval = 1, DataSplit = 100, Yweights = NULL, 
                               Prevalence = NULL, VarImport = 0, models.eval.meth = c("KAPPA", 
                                                                                      "TSS", "ROC"), SaveObj = TRUE, rescal.all.models = FALSE, 
                               do.full.models = TRUE, modeling.id = as.character(format(Sys.time(), 
                                                                                        "%s")), ...) 
{
  .Models.dependencies(silent = TRUE, models.options = models.options)
  args <- .Models.check.args(data, models, models.options, 
                             NbRunEval, DataSplit, Yweights, VarImport, models.eval.meth, 
                             Prevalence, do.full.models, SaveObj, ...)
  models <- args$models
  models.options <- args$models.options
  NbRunEval <- args$NbRunEval
  DataSplit <- args$DataSplit
  Yweights <- args$Yweights
  VarImport <- args$VarImport
  models.eval.meth <- args$models.eval.meth
  Prevalence <- args$Prevalence
  do.full.models <- args$do.full.models
  DataSplitTable <- args$DataSplitTable
  SaveObj <- args$SaveObj
  compress.arg = TRUE
  rm(args)
  models.out <- new("BIOMOD.models.out", sp.name = data@sp.name, 
                    modeling.id = modeling.id, expl.var.names = colnames(data@data.env.var), 
                    has.evaluation.data = data@has.data.eval, rescal.all.models = rescal.all.models)
  .Models.prepare.workdir(data@sp.name, models.out@modeling.id)
  if (SaveObj) {
    save(data, file = file.path(models.out@sp.name, ".BIOMOD_DATA", 
                                models.out@modeling.id, "formated.input.data"), 
         compress = compress.arg)
    models.out@formated.input.data@inMemory <- FALSE
    models.out@formated.input.data@link <- file.path(models.out@sp.name, 
                                                     ".BIOMOD_DATA", models.out@modeling.id, "formated.input.data")
    save(models.options, file = file.path(models.out@sp.name, 
                                          ".BIOMOD_DATA", models.out@modeling.id, "models.options"), 
         compress = compress.arg)
    models.out@models.options@inMemory <- FALSE
    models.out@models.options@link <- file.path(models.out@sp.name, 
                                                ".BIOMOD_DATA", models.out@modeling.id, "models.options")
  }
  mod.prep.dat <- .Models.prepare.data(data, NbRunEval, DataSplit, 
                                       Yweights, Prevalence, do.full.models, DataSplitTable)
  rm(data)
  calib.lines <- mod.prep.dat[[1]]$calibLines
  if (length(mod.prep.dat) > 1) {
    for (pa in 2:length(mod.prep.dat)) {
      calib.lines <- abind(calib.lines, mod.prep.dat[[pa]]$calibLines, 
                           along = 3)
    }
  }
  save(calib.lines, file = file.path(models.out@sp.name, ".BIOMOD_DATA", 
                                     models.out@modeling.id, "calib.lines"), compress = compress.arg)
  models.out@calib.lines@inMemory <- FALSE
  models.out@calib.lines@link <- file.path(models.out@sp.name, 
                                           ".BIOMOD_DATA", models.out@modeling.id, "calib.lines")
  rm(calib.lines)
  .Models.print.modeling.summary(mod.prep.dat, models)
  
  require(snowfall)
  sfInit(parallel=T, cpu=20)
  modeling.out <- sfLapply(mod.prep.dat, .Biomod.Models.loop,
                           modeling.id = models.out@modeling.id, Model = models,
                           Options = models.options, VarImport = VarImport, mod.eval.method = models.eval.meth,
                           SavePred = SaveObj, scal.models = rescal.all.models)
  sfStop(nostop=F)
  
  # require(parallel)
  # elves <- makeCluster(20); clusterExport(elves, c('mod.prep.dat',
  #                                                  '.Biomod.Models.loop',
  #                                                  'models.out',
  #                                                  'models',
  #                                                  'models.options',
  #                                                  'VarImport',
  #                                                  'models.eval.meth',
  #                                                  'SaveObj',
  #                                                  'rescal.all.models'), envir = environment()); clusterEvalQ(elves, library(biomod2))
  # 
  # modeling.out <- parLapply(cluster = elves, mod.prep.dat, .Biomod.Models.loop, 
  #                          modeling.id = models.out@modeling.id, Model = models, 
  #                          Options = models.options, VarImport = VarImport, mod.eval.method = models.eval.meth, 
  #                          SavePred = SaveObj, scal.models = rescal.all.models)
  # closeAllConnections()
  
  models.out@models.computed <- .transform.outputs.list(modeling.out, 
                                                        out = "models.run")
  models.out@models.failed <- .transform.outputs.list(modeling.out, 
                                                      out = "calib.failure")
  if (SaveObj) {
    models.evaluation <- .transform.outputs.list(modeling.out, 
                                                 out = "evaluation")
    save(models.evaluation, file = file.path(models.out@sp.name, 
                                             ".BIOMOD_DATA", models.out@modeling.id, "models.evaluation"), 
         compress = compress.arg)
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- file.path(models.out@sp.name, 
                                                   ".BIOMOD_DATA", models.out@modeling.id, "models.evaluation")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)
    if (VarImport > 0) {
      variables.importances <- .transform.outputs.list(modeling.out, 
                                                       out = "var.import")
      save(variables.importances, file = file.path(models.out@sp.name, 
                                                   ".BIOMOD_DATA", models.out@modeling.id, "variables.importance"), 
           compress = compress.arg)
      models.out@variables.importances@inMemory <- TRUE
      models.out@variables.importances@link <- file.path(models.out@sp.name, 
                                                         ".BIOMOD_DATA", models.out@modeling.id, "variables.importance")
      models.out@variables.importances@val <- variables.importances
      rm(variables.importances)
    }
    models.prediction <- .transform.outputs.list(modeling.out, 
                                                 out = "prediction")
    save(models.prediction, file = file.path(models.out@sp.name, 
                                             ".BIOMOD_DATA", models.out@modeling.id, "models.prediction"), 
         compress = compress.arg)
    models.out@models.prediction@inMemory <- FALSE
    models.out@models.prediction@link <- file.path(models.out@sp.name, 
                                                   ".BIOMOD_DATA", models.out@modeling.id, "models.prediction")
    rm(models.prediction)
    models.prediction.eval <- .transform.outputs.list(modeling.out, 
                                                      out = "prediction.eval")
    save(models.prediction.eval, file = file.path(models.out@sp.name, 
                                                  ".BIOMOD_DATA", models.out@modeling.id, "models.prediction.eval"), 
         compress = compress.arg)
    models.out@models.prediction.eval@inMemory <- FALSE
    models.out@models.prediction.eval@link <- file.path(models.out@sp.name, 
                                                        ".BIOMOD_DATA", models.out@modeling.id, "models.prediction.eval")
    rm(models.prediction.eval)
  }
  rm(modeling.out)
  models.out@link <- file.path(models.out@sp.name, paste(models.out@sp.name, 
                                                         ".", models.out@modeling.id, ".models.out", sep = ""))
  assign(x = paste(models.out@sp.name, ".", models.out@modeling.id, 
                   ".models.out", sep = ""), value = models.out)
  save(list = paste(models.out@sp.name, ".", models.out@modeling.id, 
                    ".models.out", sep = ""), file = models.out@link)
  .bmCat("Done")
  return(models.out)
}