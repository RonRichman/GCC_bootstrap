#Function Script

#Function to create CC Input for each AY
ccInputs <- function(apvec, current_ap, exposure, latest, development, include){
  return(list(
    accident_period = current_ap,
    apvec = abs(current_ap - apvec),
    exposure = exposure,
    latest = latest,
    development =development, 
    include = include
  )
  )
}

#Function to calculate ccIELR
ccIELR <- function(apdiff, decay, accident_period, exposure, latest, development, include){
  decay_weight <- decay ^ apdiff
  cc_numerator <- sum(latest * include * decay_weight)
  cc_denominator <- sum(exposure * development* include * decay_weight)
  cc_IELR <- cc_numerator/cc_denominator
  return(list(
    accident_period = accident_period,
    decayfactor = decay,
    ccIELR = cc_IELR
  )
  )
}

# Helper functions for Chainladder package and dplyr interaction
triangleToDf <- function(triangle, replace_na = TRUE){
  df <- as.data.frame(triangle) %>% 
    tidyr::spread(dev, value) %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames(., var="origin") 
  if(replace_na == TRUE){
    df <- df %>% replace(is.na(.), 0) 
  }
  return(df)
}

# Helper functions for Chainladder package and dplyr interaction
dfToTriangle <- function(df){
  return(df %>%
           as.matrix() %>% 
           ChainLadder::as.triangle()
  )
}

#Function that estimate cumulative and incremental triangles based on cape cod inputs
ccEstimate <- function(ccDF, triange_names){
  #Get Exposure*IELR vector
  vec_exposure_IELR <- ccDF$Exposure*ccDF$ccIELR*ccDF$Include
  #Replicate Exposure*IELR vector into matrix dimensions
  mat_exposure_IELR <- matrix(rep(vec_exposure_IELR,length(vec_exposure_IELR)), ncol = length(vec_exposure_IELR), byrow = FALSE)
  vec_development <- rev(ccDF$Developed)
  mat_developed <-matrix(rep(vec_development,length(vec_development)), ncol = length(vec_development), byrow = TRUE)
  Cape_Cod_Cumulative_Estimation = as.data.frame(mat_developed * mat_exposure_IELR)
  dimnames(Cape_Cod_Cumulative_Estimation) <- list(origin = triange_names$origin, dev = triange_names$dev)
  Cape_Cod_Incremental_Estimation <- dfToTriangle(Cape_Cod_Cumulative_Estimation) %>% 
    ChainLadder::cum2incr() %>% 
    triangleToDf()
  return(list(
    Cape_Cod_Object = ccDF,
    Cape_Cod_Cumulative_Estimation = Cape_Cod_Cumulative_Estimation,
    Cape_Cod_Incremental_Estimation = Cape_Cod_Incremental_Estimation
  )
  )
}
#Function that estimates past triangle only up to final development point
ccEstimatePastTriangle <- function(ccobj, actual_triangle){
  inc_triangle_ind <- actual_triangle %>% mutate(across(.fns = ~ifelse(is.na(.), 0, 1)))
  estimate_past_triangle <- ccobj$Cape_Cod_Incremental_Estimation * inc_triangle_ind
  return (estimate_past_triangle)
}

#Function that creates pearson adjusted residual triangle to sample from & pearson scale paramater
ccAdjPearsonResidual <- function(ccobj, actual_triangle){
  estimate_past_triangle <- ccEstimatePastTriangle(ccobj, actual_triangle)
  unscaled_pearson_residual <- (actual_triangle - estimate_past_triangle)/(abs(estimate_past_triangle))^0.5
  #Cleaning Residual to replace NaN, Inf etc. with NA
  is.na(unscaled_pearson_residual)<-sapply(unscaled_pearson_residual, is.infinite)
  unscaled_pearson_residual[is.na(unscaled_pearson_residual)]<-NA

  n <- length(na.omit(unscaled_pearson_residual[,1]))
  nobs  <- apply(unscaled_pearson_residual, 2, function(x) length(which(!is.na(x)))) %>% sum()
  scale.factor <- (nobs - 2*n+1)
  pearson_scale <- sum(unscaled_pearson_residual^2,na.rm=TRUE)/scale.factor
  ## Adjust the Pearson residuals using
  adjusted_pearson_residual <- unscaled_pearson_residual * sqrt(nobs/scale.factor)
  adjusted_pearson_residual <- adjusted_pearson_residual %>% as.matrix()
  adjusted_pearson_residual <- adjusted_pearson_residual - mean(adjusted_pearson_residual, na.rm = TRUE)
  adjusted_pearson_residual <- adjusted_pearson_residual %>% as.data.frame()
  ccobj$adjusted_pearson_residual <- adjusted_pearson_residual
  ccobj$pearson_scale <- pearson_scale
  ccobj$Cape_Cod_Incremental_Past_Estimation <- estimate_past_triangle
  
  
  return(
    ccobj
  )
}

#Performs 1 bootstrap to be called multiple times via map
ccBootstrap <- function(decay, exposure, include, adjusted_pearson_residual, estimate_past_triangle, actual_inc_triangle, best_estimate_IBNR, best_estimate_Ult, sim_num, pearson_scale, process.distr, verbrose, LDF_threshold, DP_cutoff, re_reserve){
  is_stable_LDF <- FALSE
  count_while <- 0
  while (is_stable_LDF == FALSE & count_while < 15){
    count_while <- count_while + 1
    #While loop to continue simulating until stable LDFs
    pseudo_incremental_tri <- randomClaims(exp.clms = estimate_past_triangle,
                                           resids = adjusted_pearson_residual,
                                           n.sims = 1)
    pseudo_incremental_tri <- pseudo_incremental_tri %>% replace(is.na(.),0)
    pseudo_incremental_tri[row(pseudo_incremental_tri) + col(pseudo_incremental_tri) - 1 > length(pseudo_incremental_tri)] <- NA
    pseudo_incremental_tri <- dfToTriangle(pseudo_incremental_tri)
    pseudo_incremental_tri_cum <- ChainLadder::incr2cum(pseudo_incremental_tri)
    sim_linkratio <- c(attr(ata(pseudo_incremental_tri_cum), "vwtd"),1) %>% replace(is.na(.), 1)
    sim_dev <- 1 / cumprod(rev(sim_linkratio))

    #While suspension condition testing
    test_ata = data.table(LDF = sim_linkratio)
    test_ata[,DP := 1:.N]
    test_ata[, diff2 := c(0,0,diff(LDF, 2))]
    if( (test_ata[DP>DP_cutoff, sum(ifelse(abs(diff2) > LDF_threshold,1 ,0))] >= 1)){
      is_stable_LDF <- FALSE
    }else{
      is_stable_LDF <- TRUE
    }
      
  }
  names(sim_dev) <- colnames(pseudo_incremental_tri_cum)
  sim_latest <- pseudo_incremental_tri_cum %>% getLatestCumulative()
  
  #Build up ccDF to call function ccEstimate
  sim_cc_df <- ccBuildInputDF(sim_dev, exposure, sim_latest, include, decay) %>% bind_cols(., data_frame(Include = include))
  
  #Estimating future triangle
  sim_cc_obj <- ccEstimate(sim_cc_df, dimnames(pseudo_incremental_tri_cum))
  sim_IBNR <- sim_cc_obj$Cape_Cod_Incremental_Estimation
  sim_IBNR[row(sim_IBNR) + col(sim_IBNR) - 2 < length(sim_IBNR)] <- NA
  
  #Add Process Error
  tryCatch(
    # This is what I want to do...
    {
      if (process.distr == "od.pois"){
        sim_IBNR <- apply(sim_IBNR,c(1,2), function(x) ifelse(is.na(x), NA, sign(x)*rpois.od(1, abs(x), pearson_scale)))  
      }else{
        sim_IBNR[!is.na(sim_IBNR)] <- sign(sim_IBNR[!is.na(sim_IBNR)])*rgamma(length(sim_IBNR[!is.na(sim_IBNR)]), shape=abs(sim_IBNR[!is.na(sim_IBNR)]/pearson_scale), scale=pearson_scale)
      }
      
    },
    # ... but if an error occurs, tell me what happened: 
    error=function(error_message) {
      message("Error when adding process error")
    }
  )
  
  #Re reserving logic
  if (re_reserve){
    # Build list with future incremental triangles based on process error simulated above
    sim_inc_tri_re_reserving <- re_reserving_incremental_triangles(sim_IBNR, actual_inc_triangle%>%
                                                                              tail(ncol(sim_IBNR)) %>%
                                                                              select(1:ncol(sim_IBNR))
                                                                   )
    # Recalc LDFs and shift one on wards 
    rr_Developed <- list()
    for (i in 1:length(sim_inc_tri_re_reserving)){
      rr_Developed[[i]] <- sim_inc_tri_re_reserving[[i]] %>% incr2cum() %>% ata() %>% attr("vwtd") %>% replace(is.na(.), 1)
      rr_Developed[[i]] <- c(1,1/(rr_Developed[[i]] %>% rev %>% cumprod))%>% as.data.frame() %>% setNames("Developed") %>% tibble::remove_rownames()
    }
    #Build Data Frame of CC output for each future step
    re_reserving_cape_cod <- rereserving_cape_cod(sim_inc_tri_re_reserving, exposure, rr_Developed, decay)$re_reserving_cc_obj
  }

  
  #Building up outputs
  sim_IBNR_total <- sum(sim_IBNR, na.rm = TRUE)
  sim_ult <- sum(sim_latest) + sim_IBNR_total
  
  sim_output <- cbind(sim_num, sim_cc_df, sim_IBNR_total, sim_ult, process.distr, best_estimate_IBNR, best_estimate_Ult ) %>% 
                  dplyr::rename(
                    `Simulation Number` = sim_num,
                    Decay = decayfactor,
                    `Accident Period` = accident_period,
                    `Sim IBNR` = sim_IBNR_total,
                    `Sim Ult` = sim_ult,
                    `BE IBNR` = best_estimate_IBNR,
                    `BE Ult` = best_estimate_Ult, 
                    `Process Dist` = process.distr
                  )
  
  #Re reserving logic for output
  if (re_reserve){
    sim_output <-  sim_output %>% 
                    left_join(., re_reserving_cape_cod, by = c("Accident Period" = "CDR_n")) %>% 
                    mutate(CDR = rr_Ult - `BE Ult`)
  }
  
  if(verbrose){
    return(list(sim_output = list(sim_output),
                process_tri = sim_IBNR,
                pseudo_incremental_tri = pseudo_incremental_tri,
                sim_linkratio = sim_linkratio,
                disregarded_sims = count_while)
    )
  } else{
    return(list(sim_output))
  } 
  
}

#Extract future location from a triangle
future_triangle_position <- function(triangle_df){
  return(
    row(triangle_df) + col(triangle_df) - 1 > length(triangle_df)
  )
}

#Extract past location from a triangle
past_triangle_position <- function(triangle_df){
  return(
    row(triangle_df) + col(triangle_df) - 2 < length(triangle_df)
  )
}

#Bootstrap wrapper function
ccBootstrapWrapper <- function(n=2000, ccobj,actual_inc_triangle, best_estimate_IBNR, best_estimate_Ult, process.distr=c("gamma", "od.pois"), verbrose = FALSE, LDF_threshold, DP_cutoff, re_reserve = FALSE){
  if (verbrose){
    1:n %>% furrr::future_map(~ccBootstrap(decay = ccobj$Cape_Cod_Object$decayfactor %>% head(1),
                                           exposure = ccobj$Cape_Cod_Object$Exposure,
                                           include = rep(1,length(ccobj$Cape_Cod_Object$Exposure)),
                                           adjusted_pearson_residual = ccobj$adjusted_pearson_residual,
                                           estimate_past_triangle = ccobj$Cape_Cod_Incremental_Past_Estimation,
                                           actual_inc_triangle = actual_inc_triangle,
                                           best_estimate_IBNR = ccobj$Cape_Cod_Object$ccIBNR %>% sum(),
                                           best_estimate_Ult = ccobj$Cape_Cod_Object$ccUlt %>% sum(),
                                           sim_num = .x,
                                           pearson_scale = ccobj$pearson_scale,
                                           process.distr = process.distr,
                                           verbrose = verbrose, 
                                           LDF_threshold, 
                                           DP_cutoff,
                                           re_reserve
    ), .options = furrr_options(seed = TRUE), 
    .progress = TRUE
    )
  } else{
    1:n %>% furrr::future_map_dfr(~ccBootstrap(decay = ccobj$Cape_Cod_Object$decayfactor %>% head(1),
                                           exposure = ccobj$Cape_Cod_Object$Exposure,
                                           include = rep(1,length(ccobj$Cape_Cod_Object$Exposure)),
                                           adjusted_pearson_residual = ccobj$adjusted_pearson_residual,
                                           estimate_past_triangle = ccobj$Cape_Cod_Incremental_Past_Estimation,
                                           actual_inc_triangle = actual_inc_triangle,
                                           best_estimate_IBNR = ccobj$Cape_Cod_Object$ccIBNR %>% sum(),
                                           best_estimate_Ult = ccobj$Cape_Cod_Object$ccUlt %>% sum(),
                                           sim_num = .x,
                                           pearson_scale = ccobj$pearson_scale,
                                           process.distr = process.distr,
                                           verbrose = verbrose, 
                                           LDF_threshold, 
                                           DP_cutoff,
                                           re_reserve
    ), .options = furrr_options(seed = TRUE), 
    .progress = TRUE
    )
  }
}

#function to generate od.pois process error
rpois.od<-function (n, lambda,d=1) {
  if (d==1 | lambda<=0)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}

sampleResiduals <- function(resids, positions, n.sims){
  ## Author: Nigel de Silva
  ## Worry about excluding, zoning, etc. residuals later
  
  resids <- as.vector(resids) %>% unlist()
  resids <- resids[!is.na(resids)]
  
  .sample <- sample(resids, prod(dim(positions)[-3])*n.sims, replace=T)
  .sample <- array(.sample, dim=c(dim(positions)[-3], n.sims))
  .sample[is.na(positions)] <- NA
  
  return(.sample)
}

randomClaims <- function(exp.clms, resids, n.sims){
  ## Author: Nigel de Silva
  .residSample <- sampleResiduals(resids, exp.clms %>% replace(. == 0, NA), n.sims)
  out <- .residSample * sqrt(abs(exp.clms)) + exp.clms
  
  return(out)
}

re_reserving_incremental_triangles <- function(sim_IBNR, actual_inc_triangle){
  sim_inc_tri_rereserving <- list()
  
  #For loop that starts from Dev 2 up to development
  for (i in 2: nrow(sim_IBNR)){
    tmp_inc <- actual_inc_triangle%>% replace(is.na(.), 0) %>% as.matrix()
    tmp_cdr <- sim_IBNR
    tmp_cdr[row(sim_IBNR) + col(sim_IBNR) -1 - i >= nrow(sim_IBNR)] <- 0
    tmp_cdr <- tmp_cdr %>% replace(is.na(.), 0) %>%  as.matrix()
    tmp_inc <- tmp_inc + tmp_cdr
    tmp_inc <- tmp_inc %>% as.data.frame()
    tmp_inc[row(tmp_inc) + col(tmp_inc) -1 - i >= nrow(tmp_inc)] <- NA
    sim_inc_tri_rereserving[[paste0("CDR ", i-1)]] <- tmp_inc %>% dfToTriangle()
  }
  
  return(sim_inc_tri_rereserving)
}


rereserving_cape_cod <- function(sim_inc_tri_rereserving, exposure, rr_Developed_list, decay){
  re_reserving_cc_obj <- list()
  re_reserving_cc_DF <- list()
  for (i in 1 :length(sim_inc_tri_rereserving)){
    rr_developed <- lag(rr_Developed_list[[i]]$Developed, i, default = 1)
    #Add latest and do normal IBNR = ccIELR x (1-Dev) x Exposure calculation
    tmp_latest <- sim_inc_tri_rereserving[[i]] %>% 
      incr2cum() %>% 
      getLatestCumulative() %>% 
      as.data.frame() %>% 
      setNames("Latest") %>% tibble::remove_rownames()
    #Build up ccObj input for re-reserving
    sim_cc_df <- ccBuildInputDF(rr_developed, exposure, tmp_latest, include = rep(1, length(rr_developed)), decay)
    #Re-reserve Deterministically
    sim_cc_df <- bind_cols(sim_cc_df, tmp_latest)%>% 
                  mutate(
                    ccIBNR = `ccIELR` * `Exposure` * (1 - `Developed`),
                    ccUlt  = `ccIBNR` + `Latest`
                  ) %>% 
                  tibble::remove_rownames()
     
    
    sim_cc_df_total <- sim_cc_df %>% 
                            summarise(
                              rr_IBNR = sum(ccIBNR),
                              rr_Ult =sum(ccUlt),
                              CDR_n = i
                            )
    
    
    #Sum IBNR & Ult
    re_reserving_cc_obj[[i]] <- sim_cc_df_total
    re_reserving_cc_DF[[i]] <- sim_cc_df
    
  }
  
  re_reserving_cc_obj <- bind_rows(re_reserving_cc_obj)
  
  return(list (re_reserving_cc_obj = re_reserving_cc_obj, re_reserving_cc_DF = re_reserving_cc_DF))
}

ccBuildInputDF <- function (sim_dev, exposure, sim_latest, include, decay){
  
  #Set up input for ccIELR calculation
  sim_cc_input_list <-1:length(sim_dev) %>%
    purrr::map(~ ccInputs(1:length(sim_dev), 
                          .x, 
                          exposure, 
                          sim_latest, 
                          sim_dev, 
                          include)
    )
  #Calculate ccIELR for simulation 
  sim_cc_IELR <- sim_cc_input_list %>% 
    purrr::map_df(~ ccIELR(.x$apvec, 
                           decay, 
                           .x$accident_period,
                           .x$exposure, 
                           .x$latest, 
                           .x$development, 
                           .x$include)
    )
  #Build up ccDF to call function ccEstimate
  sim_cc_df <- data.frame(
    Exposure = exposure,
    sim_cc_IELR,
    Developed = sim_dev
  )
  
  return(sim_cc_df)
}

