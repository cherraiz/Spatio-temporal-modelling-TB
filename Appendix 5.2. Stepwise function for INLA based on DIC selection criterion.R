INLAdf_stepBLogit<-function(fam1 = "binomial",
                            dataf,
                            invariant = "0 + Intercept",
                            direction = c("forwards","backwards"),
                            num_trials = NULL,
                            y = NULL,
                            y2 = NULL,
                            include = NULL,
                            n_threads = inla.getOption("num.threads"),
                            powerl = 1,
                            inter = 1,
                            thresh = 5,
                            ...) {
  
  # Basic checks
  if (is.null(nrow(dataf))) {
    stop("error in data frame")
  }
  if (nrow(dataf) == 0) {
    stop("no rows in data frame")
  }
  if (is.null(y)) {
    stop("no y variable")
  }
  if (!(class(dataf) == "data.frame" |
        class(dataf) == "SpatialPointsDataFrame")) {
    stop("data is not a data frame")
  }
  
  if (is.null(include))
    include <- (1:ncol(dataf))[!names(dataf) %in% c(y, y2, Ntrials)]
  
  z <- NULL
  
  facts <- sapply(dataf, is.factor)[include]
  explF<-names(dataf)[include]
  expl<-explF[!facts]
  expl <- expandExplanatoryVars(expl, explF, facts, powerl, inter)
  choice <- NULL
  chosen <- NULL
  new1 <- NULL
  dicloss <- 999
  dicold <- NULL
  
  while(length(expl) > 0){
    # If backwards.... ? 
    if (direction == "backwards") {
      runs <- c(1:length(expl), 9999999)
    } else{
      runs <- 1:length(expl)
    }
    
    for (ii in runs) {
      if (direction == "backwards") {
        if (ii == 9999999) {
          ii <- 1:length(expl)
        } else{
          ii <-
            {
              -1 * ii
            }}}
      
      if (is.null(chosen)) {
        if (length(expl[ii]) > 0) {
          formula2 <-
            formula(paste(y, "~", invariant, "+", paste(expl[ii], collapse = "+"), sep =
                            ""))
        } else {
          formula2 <- formula(paste(y, "~", invariant))
        }
      } else{
        if (length(expl[ii]) > 0) {
          formula2 <-
            formula(paste(y, "~", invariant, "+", chosen, " + ", expl[ii], sep = ""))
        } else {
          formula2 <- formula(paste(y, "~", invariant, "+", chosen))
        } }
      
      result2 <- INLA::inla(formula2,
                            data = dataf,
                            family = fam1,
                            Ntrials = dataf[,c(num_trials)],
                            control.family = list(control.link=list(model="logit")),
                            control.predictor = list(compute =TRUE, link = 1),
                            control.compute = list(dic = TRUE, waic =TRUE, config=TRUE, cpo = TRUE),
                            num.threads = n_threads
      )
      
      rmse <-
        sqrt(mean((
          dataf[, y2] - result2$summary.fitted.values$mean[1:nrow(dataf)]
        ) ^ 2, na.rm = TRUE))
      sumcpo <- sum(log(result2$cpo$cpo), na.rm = TRUE)
      if (length(ii) > 1) {
        var1 <- paste(expl[ii], collapse = "+")
      } else{
        var1 <- expl[abs(ii)]
      }
      if (is.null(choice)) {
        choice <-
          data.frame(
            var = var1,
            dic = result2$dic$dic,
            rmse,
            sumcpo,
            stringsAsFactors = FALSE
          )
      } else{
        choice <-
          rbind(
            choice,
            data.frame(
              var = var1,
              dic = result2$dic$dic,
              rmse,
              sumcpo,
              stringsAsFactors = FALSE
            )
          )
      }
    } 
    
    new1 <- choice[which.min(choice$dic), 1]
    
    if (!is.null(dicold)) {
      dicloss <- dicold - min(choice$dic, na.rm = TRUE)[1]
    }
    
    dicold <- choice[which.min(choice$dic), 2]
    
    if (is.null(z)) {
      progress <- choice[choice$var == new1, ]
      z <- 1
    } else {
      progress <- rbind(progress, choice[choice$var == new1, ])
    }
    
    message(paste(new1, " - ", min(choice$dic, na.rm = TRUE)), sep = "")
    choice <- NULL
    if (dicloss > thresh) {
      
      if (direction == "backwards") {
        expl <- expl[!expl == new1]
      }
      if (direction == "forwards") {
        if (is.null(chosen)) {
          chosen <- new1
          expl <- expl[!expl == new1]
        } else {
          chosen <- paste(chosen, " + ", new1, sep = "")
          expl <- expl[!expl == new1]
        }
      }
    } else {
      break
    }
    
  }
  
  if (direction == "backwards") {
    formulax <-
      formula(paste(y, "~", invariant, "+", paste(expl[ii], collapse = "+"), sep = ""))
  } else {
    if(!is.null(chosen)){
      formulax <- formula(paste(y, "~", invariant, "+", chosen, sep = ""))
    } else {
      formulax <- formula(paste(y, "~", invariant, sep = ""))
    }
  }
  
  resultx <- INLA::inla(formulax,
                        data = dataf,
                        family = fam1,
                        Ntrials = dataf[,c(num_trials)],
                        control.family = list(control.link=list(model="logit")),
                        control.predictor = list(compute =TRUE, link = 1),
                        control.compute = list(dic = TRUE, waic =TRUE, config=TRUE, cpo = TRUE),
                        num.threads = n_threads
  )
  
  
  
  output <- list(
    best_formula = formulax,
    dic = dicold,
    progress = progress,
    best_model = resultx
  )
  class(output) <- 'INLAstep'
  
  return(output)
  
}
expandExplanatoryVars <- function(expl, explF, facts, powerl, inter){
  
  if (length(expl) > 0) {
    if (powerl == 2) {
      expl2 <-
        paste("I(", expl, "^2)", sep = "")
      expl3 <- NULL
      expl4 <- NULL
      explX <- c(expl, expl2, expl3, expl4)
    }
    if (powerl == 3) {
      expl2 <-
        paste("I(", expl, "^2)", sep = "")
      expl3 <-
        paste("I(", expl, "^3)", sep = "")
      expl4 <- NULL
      explX <- c(expl, expl2, expl3, expl4)
    }
    if (powerl >= 4) {
      expl2 <-
        paste("I(", expl, "^2)", sep = "")
      expl3 <-
        paste("I(", expl, "^3)", sep = "")
      expl4 <-
        paste("I(", expl, "^4)", sep = "")
      explX <- c(expl, expl2, expl3, expl4)
    }
    if (powerl > 1) {
      expl <- explX
    }
    if (inter >= 2) {
      lvls <- data.frame(p1 = utils::combn(expl, 2)[1, ],
                         p2 = utils::combn(expl, 2)[2, ])
      lvls2 <- do.call(paste, c(lvls[names(lvls)], sep = ":"))
      expl2 <- c(expl, lvls2)
    }
    if (inter >= 3) {
      lvls <- data.frame(
        p1 = utils::combn(expl, 3)[1, ],
        p2 = utils::combn(expl, 3)[2, ],
        p3 = utils::combn(expl, 3)[3, ]
      )
      lvls3 <- do.call(paste, c(lvls[names(lvls)], sep = ":"))
      expl2 <- c(expl2, lvls3)
    }
    if (inter >= 4) {
      lvls <- data.frame(
        p1 = utils::combn(expl, 4)[1, ],
        p2 = utils::combn(expl, 4)[2, ],
        p3 = utils::combn(expl, 4)[3, ],
        p4 = utils::combn(expl, 4)[4, ]
      )
      lvls4 <- do.call(paste, c(lvls[names(lvls)], sep = ":"))
      expl2 <- c(expl2, lvls4)
    }
    if (inter > 1) {
      expl <- expl2
    }
  }
  if (length(explF[facts]) > 0) {
    expl <- c(expl, explF[facts])
  }
  return(expl)
}
