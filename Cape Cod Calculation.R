#1. Loading packages

packages <- c("matrixStats",
              "dplyr",
              "rio",
              "purrr",
              "ChainLadder",
              "furrr",
              "tictoc",
              "plotly",
              "data.table",
              "lattice",
              "openxlsx",
              "tools")

install.packages(setdiff(packages, rownames(installed.packages())))


require(matrixStats)
require(rio)
require(purrr)
require(ChainLadder)
require(furrr)
require(tictoc)
require(plotly)
require(data.table)
require(lattice)
require(DBI)
require(RODBC)
require(dplyr)
require(openxlsx)
require(tools)
require(janitor)

source("./CapeCod Functions.R")

#2. Variables for Best Estimate, Bootstrap & Re-reserving

#2.1 set all of the parameters for the reserving exercise
### Best Estimate Inputs####
long_format <- FALSE #Triangle or long claims data
increments <- 0.1
decay <- c(0.000001,seq(from = 0.1, to = 1, by = increments)) #Decay factor grid

### Simulation Inputs####
sims = 3000
process.distr = "od.pois" #gamma or od.pois
verbrose = FALSE
LDF_threshold = 0.05
DP_cutoff = 8
re_reserve = FALSE

#2.2 read in premium & claims (triangle or long)
claims_data_list <-list.files(path = "./Claims/",
                                pattern = "*.xlsx",
                                full.names = TRUE) %>% 
                     map(~rio::import(.))

names(claims_data_list) <-    list.files(path = "./Claims/",
                                           pattern = "*.xlsx") %>% 
                              gsub(pattern = "\\.xlsx$", "", .)

premium_data_list <-list.files(path = "./Premium/",
                               pattern = "*.xlsx",
                               full.names = TRUE) %>% 
                    map(~rio::import(.))

names(premium_data_list) <- list.files(path = "./Premium/",
                                       pattern = "*.xlsx") %>% 
                            gsub(pattern = "\\.xlsx$", "", .)

#2.3 Run the reserving logic which will
    #Perform deterministic grid Cape Cod (best estimate)
    #Bootstrap Cape Cod
    #Re-reserve if required for each bootstrap

for (z in 1:length(claims_data_list)){
  
  if(long_format){
    #2. Loading required data
    dat = claims_data_list[[z]] %>% as.data.table()
    dat[, latest := max(dev), keyby = .(origin)]
    dat[, latest_ind := latest == dev]
    
    inputs = dat[latest_ind == T][, c(1,4,5)]
    inputs[, `AP Number` := as.integer(as.factor(origin))]
    inputs[, Include := 1]
    
    triang = dat %>% as.triangle(origin = "origin", dev = "dev", value = "claim")
    LDFs= triang %>% ata() %>% attr("vwtd")
    perc_dev = 1/(LDFs %>% rev %>% cumprod)
    
    inputs[, Developed := c(1,perc_dev)]
    inputs %>% setnames(names(inputs), c("Accident Period",
                                         "Latest",
                                         "Exposure",
                                         "AP Number",
                                         "Include",
                                         "Developed"))
    
    data = inputs
    cum_triangle_data = triang %>% triangleToDf(replace_na = FALSE)
    
    #Remove 1st Column as this is origin info
    cum_triangle <- cum_triangle_data
    # Convert to matrix
    cum_triangle <- as.matrix(cum_triangle)
    # Add dimension names
    dimnames(cum_triangle) <- list(origin=row.names(cum_triangle_data), dev=0:(ncol(cum_triangle)-1))
    #Convert to chainladder triangle object for calculating incurred
    cum_triangle <- ChainLadder::as.triangle(cum_triangle)
    input_inc_triangle <- ChainLadder::cum2incr(cum_triangle) %>% triangleToDf(. ,replace_na = FALSE)
  }else{
    #2. Loading required data
    data <- premium_data_list[[z]]
    #Read in Cum Paid/Incurred traingle
    cum_triangle_data <- claims_data_list[[z]]
    #Remove 1st Column as this is origin info
    cum_triangle <- cum_triangle_data[,-1]
    # Convert to matrix
    cum_triangle <- as.matrix(cum_triangle)
    # Add dimension names
    dimnames(cum_triangle) <- list(origin=cum_triangle_data[,1], dev=1:ncol(cum_triangle))
    #Convert to chainladder triangle object for calculating incurred
    cum_triangle <- ChainLadder::as.triangle(cum_triangle)
    input_inc_triangle <- ChainLadder::cum2incr(cum_triangle) %>% triangleToDf(. ,replace_na = FALSE)
    inc_triangle_ind <- input_inc_triangle %>% mutate(across(.fns = ~ifelse(is.na(.), 0, 1)))
    #Build up Accident Period Number, Latest, Dev% and Include Indicator
    Latest <- ChainLadder::as.triangle(cum_triangle) %>% getLatestCumulative() %>% as.data.frame() %>% setNames("Latest")
    Include <- data %>% select(Include)
    AP = 1:nrow(Latest) %>% as.data.frame() %>% setNames("AP Number")
    LDFs= cum_triangle %>% ata() %>% attr("vwtd")
    Developed = c(1,1/(LDFs %>% rev %>% cumprod))%>% as.data.frame() %>% setNames("Developed")
    
    
    data <- bind_cols(data, Latest, AP, Developed) %>% 
      setNames(c("Accident Period",
                 "Exposure",
                 "Include",
                 "Latest",
                 "AP Number",
                 "Developed")
      ) %>%
      select(
        `Accident Period`,
        `AP Number`,
        `Latest`,
        `Exposure`,
        `Developed`,
        `Include`
      ) %>% 
      tibble::remove_rownames()
  }
  
  #Calculating CC weights, CC IELR and creates a DF that stores CC best estimate results
  cc_input_list <- data$`AP Number` %>%
    purrr::map(~ ccInputs(data$`AP Number`, 
                          .x, 
                          data$Exposure, 
                          data$Latest, 
                          data$Developed, 
                          data$Include)
    )
  
  cross_data <- list(
    cc_inputs = cc_input_list,
    decay_factor = decay
  ) %>% 
    purrr::cross() %>% 
    purrr::map_df(~ ccIELR(.x$cc_inputs$apvec, 
                           .x$decay_factor, 
                           .x$cc_inputs$accident_period,
                           .x$cc_inputs$exposure, 
                           .x$cc_inputs$latest, 
                           .x$cc_inputs$development, 
                           .x$cc_inputs$include)
    ) %>% 
    left_join(. , data, by = c("accident_period" = "AP Number")) %>% 
    select(
      `Accident Period`,
      Exposure,
      Latest,
      Include,
      Developed,
      decayfactor,
      ccIELR
    ) %>% 
    group_by(decayfactor) %>% 
    mutate(
      ccIBNR = `ccIELR` * `Exposure` * (1 - `Developed`),
      ccUlt  = `ccIBNR` + `Latest`
    )
  #Summarising results by decay factor
  summary_data <- cross_data %>% 
    summarise(
      IBNR = sum(ccIBNR, na.rm = TRUE),
      Ult = sum(ccUlt, na.rm = TRUE),
    )
  
  cross_list <- cross_data %>%  group_split(.keep = TRUE)
  
  #Using list of CC Input data frames (with different decay factors) to estimate cumulative and incremental triangles
  if (length(dimnames(cum_triangle)[[1]]) != length(data$Exposure)){
    message("Exposure Inputs and Triangle dimentions not the same")
    stop()
  }else{
    cape_cod_estimation <- cross_list %>% purrr::map(~ccEstimate(.x, dimnames(cum_triangle)))
  }
  
  #Calculating adjusted pearson residuals, pearson scale parameter for each Cape Cod Object (where ccObj is based on different decay factors)
  
  ccobj <- list()
  for( i in 1:length(cape_cod_estimation)){
    ccobj[[i]] <-  ccAdjPearsonResidual(
      ccobj = cape_cod_estimation[[i]],
      actual_triangle = input_inc_triangle
    )
  }
  
  #Run simulation n times for each Cape Cod Object
  no_cores <- availableCores() - 1
  plan(multisession, workers = no_cores)
  
  bootstrap <- list()
  tic()
  writeLines(paste("\n #### Starting simulation for Class =",names(claims_data_list)[z], "###"))
  for( i in 1:length(ccobj)){
    writeLines(paste("\n Simulation for Decay factor=", decay[i]))
    bootstrap[[i]] <-  ccBootstrapWrapper(n = sims, 
                                          ccobj = ccobj[[i]],
                                          actual_inc_triangle = input_inc_triangle,
                                          best_estimate_IBNR = ccobj[[i]]$Cape_Cod_Object$ccIBNR %>% sum(),
                                          best_estimate_Ult = ccobj[[i]]$Cape_Cod_Object$ccUlt %>% sum(),
                                          process.distr = process.distr, 
                                          verbrose = verbrose, 
                                          LDF_threshold, 
                                          DP_cutoff,
                                          re_reserve)
    writeLines(paste("\nSimulation progress for class = ",names(claims_data_list)[z], round(i/length(decay)*100,2), "%"))
  }
  toc()
  if(verbrose){
    process_verbrose <- list()
    cc_obj_verbrose <- list()
    pseudo_incremental_tri_verbrose <- list()
    LDF_verbrose <- list()
    disregarded_sims_verbrose <- list()
    for (i in 1:length(bootstrap)){
      #loop for decay outputs
      for(j in 1:length(bootstrap[[i]]) ){
        #loop for simulation - builds list of (decay# ; sim#)
        process_verbrose[paste0(i,";",j)] <- list(bootstrap[[i]][[j]]$process_tri%>% data.frame())
        
        cc_obj_verbrose[paste0(i,";",j)] <- list(bootstrap[[i]][[j]]$sim_output %>% data.frame())
        
        pseudo_incremental_tri_verbrose[paste0(i,";",j)] <-  list(bootstrap[[i]][[j]]$pseudo_incremental_tri %>% triangleToDf())      
        
        LDF_verbrose[paste0(i,";",j)] <- list(bootstrap[[i]][[j]]$sim_linkratio)
        
        disregarded_sims_verbrose[paste0(i,";",j)] <- list(bootstrap[[i]][[j]]$disregarded_sims)
        
      }
    }
    bootstrap <- bind_rows(cc_obj_verbrose) %>% 
      dplyr::rename(
        `Simulation Number` = Simulation.Number,
        `Accident Period` = Accident.Period,
        `Process Dist` = Process.Dist,
        `Sim IBNR` = Sim.IBNR,
        `Sim Ult` = Sim.Ult,
        `BE IBNR` = BE.IBNR,
        `BE Ult` = BE.Ult
      )
    
  }else{
    bootstrap <- bind_rows(bootstrap)
  }
  
  
  bootstrap <- bootstrap %>% 
    left_join(., data %>% select(Latest, `AP Number`), by = c("Accident Period" = "AP Number")) %>% 
    mutate(
      `AY Ult` = Latest + Exposure * ccIELR * (1- Developed),
      `AY IBNR` = `AY Ult` - Latest,
      `AY ULR` = `AY Ult` / Exposure
    )
  
  
  #Stripping out AP info
  bootstrap_IBNR <- bootstrap %>% 
    group_by(Decay) %>%
    distinct_at(., .vars = vars(Latest, `Sim IBNR`, `Sim Ult`,`BE IBNR`, `BE Ult`, `Process Dist`,`Simulation Number`)) 
  
  #Stripping out re-reserving info
  if (re_reserve){
    bootstrap_re_reserving <- bootstrap %>% 
      rename(`CDR(n)` = `Accident Period`) %>% 
      group_by(Decay, `CDR(n)`) %>%
      select(
        Decay, 
        `CDR(n)`,
        CDR
      ) %>% 
      summarise(
        `Mean CDR` = mean(CDR, na.rm = TRUE),
        `SD CDR` = sd(CDR, na.rm = TRUE),
        `Median CDR` = median(CDR, na.rm = TRUE)
      ) %>% 
      replace(is.na(.), NA)

  }
  
  
  #Stripping out IELR info
  bootstrap_ccIELR <- bootstrap %>% 
    left_join(. , data %>% select(`Accident Period`,`AP Number`), by = c("Accident Period" = "AP Number")) %>% 
    rename(
      `Accident Number` = `Accident Period`,
      `Accident Period` = `Accident Period.y`
    ) %>% 
    group_by(Decay, 
             `Accident Period`,
             `Process Dist`) %>% 
    select(
      Decay,
      `Accident Period`,
      `Simulation Number`,
      ccIELR,
      `Process Dist`
    ) %>% 
    summarise(
      `Mean ccIELR` = mean(ccIELR),
      `SD ccIELR` = sd(ccIELR)) %>% 
    arrange(dplyr::across(c(`Accident Period`, Decay)), .by_group = TRUE) %>% 
    ungroup() %>% 
    group_by(Decay) %>% 
    mutate(
      `Accident Period Number` = 1:n()
    )
  
  #Summary Statistics
  IBNR_summary <- bootstrap_IBNR %>% 
    group_by(Decay, `Process Dist`) %>% 
    summarise(
      `Mean IBNR` = mean(`Sim IBNR`),
      `GCC BE IBNR` = mean(`BE IBNR`),
      `SD IBNR` = sd(`Sim IBNR`),
      `Mean Ult` = mean(`Sim Ult`),
      `GCC BE Ult` = mean(`BE Ult`),
      `50th IBNR Quantile` = quantile(`Sim IBNR`, c(0.5)),
      `75th IBNR Quantile` = quantile(`Sim IBNR`, c(0.75)),
      `90th IBNR Quantile` = quantile(`Sim IBNR`, c(0.90)),
      `99th IBNR Quantile` = quantile(`Sim IBNR`, c(0.99))
    ) %>% 
    mutate(
      `Coef of Var` = `SD IBNR`/`Mean IBNR`,
      `Risk Adjustment` = `75th IBNR Quantile` -  `Mean IBNR`
    ) %>% 
    mutate(
      Latest = sum(data$Latest)
    ) %>% 
    select(
      Decay,
      `Process Dist`,
      Latest,
      `Mean IBNR`,
      `GCC BE IBNR`,
      `SD IBNR`,
      `Mean Ult`,
      `GCC BE Ult`,
      `50th IBNR Quantile`,
      `75th IBNR Quantile`,
      `90th IBNR Quantile`,
      `99th IBNR Quantile`,
      `Risk Adjustment`,
      `Coef of Var`
    )
  
  if (re_reserve){
    CDR_ccIELR_IBNR_Summary <- bootstrap_re_reserving %>% 
      left_join(., bootstrap_ccIELR, by = c("CDR(n)" = "Accident Period Number", "Decay" = "Decay")) %>%
      left_join(., IBNR_summary %>% select(`SD IBNR`, `Mean IBNR`), by = c("Decay" = "Decay")) %>%
      left_join(., bootstrap_IBNR %>% 
                  group_by(Decay) %>% 
                  distinct( `BE IBNR`, `BE Ult`) , by = c("Decay" = "Decay")) %>%
      bind_cols(data %>% 
                  select(Exposure, Latest) %>% 
                  summarise(Exposure = sum(Exposure),
                            Latest = sum(Latest)
                  )
      )%>% 
      mutate(
        `Mean Ult` = Latest + `Mean IBNR`,
        ccULR = `Mean Ult` / Exposure,
        `BE ULR` = `BE Ult` / Exposure
      ) %>% 
      select(
        Decay,
        `CDR(n)`,
        `Mean CDR`,
        `SD CDR`,
        `SD IBNR`,
        `Accident Period`,
        Exposure,
        Latest,
        `Mean ccIELR`,
        `SD ccIELR`,
        `BE IBNR`,
        `Mean IBNR`,
        `BE Ult`,
        `Mean Ult`,
        `BE ULR`,
        ccULR,
        `Process Dist`
      ) 
    
    #Premium and Reserve Risk correlation and plots
    recent_ay = bootstrap %>%  filter(`Accident Period` == max(`Accident Period`)) %>% data.table()
    recent_ay[, res_risk := `AY ULR` - mean(`AY ULR`), by = Decay]
    recent_ay[, prem_risk := ccIELR - mean(ccIELR), by = Decay]
    recent_ay%>% ggplot(aes(x = prem_risk, y = res_risk)) + facet_wrap(~Decay)+geom_point()
    recent_ay[,cor(res_risk,prem_risk), keyby = .(Decay)]
    recent_ay%>% ggplot(aes(x = `Sim IBNR`)) + facet_wrap(~Decay)+geom_histogram()
  }
  
  #Mack, ODP & Gamma for back testing
  mack <- cum_triangle %>% replace(. ==0, NA) %>% MackChainLadder()
  od.pois <- cum_triangle %>% replace(. ==0, NA) %>% BootChainLadder(process.distr = "od.pois", R = 5000)
  gamma <- cum_triangle %>% replace(. ==0, NA) %>% BootChainLadder(process.distr = "gamma", R = 5000)
  
  chainladder_boot_IBNR_Summary <- data.frame(
    `Mean IBNR` = c(summary(od.pois)$Totals["Mean IBNR", ], 
                    summary(mack)$Totals["IBNR", ],
                    summary(gamma)$Totals["Mean IBNR", ]),
    `Mean Ultimate` = c(summary(od.pois)$Totals["Mean Ultimate", ], 
                        summary(mack)$Totals["Ultimate", ],
                        summary(gamma)$Totals["Mean Ultimate", ]),
    SD = c(summary(od.pois)$Totals["SD IBNR", ], 
           summary(mack)$Totals["Mack S.E.", ],
           summary(gamma)$Totals["SD IBNR", ]),
    CV = c(summary(od.pois)$Totals["SD IBNR", ] / summary(od.pois)$Totals["Mean IBNR", ],
           summary(mack)$Totals["Mack S.E.", ]/summary(mack)$Totals["IBNR", ],
           summary(gamma)$Totals["SD IBNR", ]/ summary(gamma)$Totals["Mean IBNR", ]) ,
    `75th` = c(summary(od.pois)$Totals["Total IBNR 75%", ], 
               NA_integer_,
               summary(gamma)$Totals["Total IBNR 75%", ]),
    RA = c(summary(od.pois)$Totals["Total IBNR 75%", ] - summary(od.pois)$Totals["Mean IBNR", ],
           NA_integer_,
           summary(gamma)$Totals["Total IBNR 75%", ] - summary(gamma)$Totals["Mean IBNR", ]),
    `Bootstrap Method` = c("OPD", 
                           "Mack",
                           "Gamma")
  )

  ##Writing Essential information - can tweak code for re-reserving and premium + reserve risk outputs
  dataset_names <- list()
  #Input
  dataset_names[['Input - Triangle']] = cum_triangle_data
  #outputs
  dataset_names[["Output - IBNR Summary"]] = IBNR_summary
  dataset_names[["Output - Back testing"]] = chainladder_boot_IBNR_Summary
  #Calculations
  for (i in 1:length(ccobj)){
    temp_key = paste0("Calc - CC Object(",decay[i],")")
    temp_df = ccobj[[i]]$Cape_Cod_Object
    dataset_names[[temp_key]] = temp_df
    
    temp_key = paste0("Calc - Adj Pearson Res (",decay[i],")")
    temp_df = ccobj[[i]]$adjusted_pearson_residual
    dataset_names[[temp_key]] = temp_df
    
    if (re_reserve){
      temp_key = paste0("Calc - CDR Results (",decay[i],")")
      temp_df = CDR_ccIELR_IBNR_Summary %>% filter(Decay == decay[i])
      dataset_names[[temp_key]] = temp_df
    }
  }
  openxlsx::write.xlsx(dataset_names, file = paste0("./Cape Cod Output/",names(claims_data_list)[z],".xlsx"))
  
  #Cleaning environment and saving workspace image for the run
  save.image(file = paste0("./Environments/",names(claims_data_list)[z]), safe = FALSE)
  rm(list=setdiff(ls(), c("claims_data_list",
                          "premium_data_list",
                          "sims",
                          "process.distr",
                          "verbrose",
                          "LDF_threshold",
                          "DP_cutoff",
                          "re_reserve",
                          "long_format",
                          "increments",
                          "decay",
                          "z"
                  )
                  )
     ) #will clear all objects except input variables.
  gc() #free up memrory and report the memory usage.
  source("./CapeCod Functions.R")
  require(matrixStats)
  require(rio)
  require(purrr)
  require(ChainLadder)
  require(furrr)
  require(tictoc)
  require(plotly)
  require(data.table)
  require(lattice)
  require(DBI)
  require(RODBC)
  require(dplyr)
  require(openxlsx)
  require(tools)
  require(janitor)
}

