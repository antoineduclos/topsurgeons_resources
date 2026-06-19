/************************************************
 * MACRO: funnel_plot

 * Parameters:
 *   data      = Input dataset with predictions (must contain: an_trim, eta, event)
 *   quarter   = Quarter to plot (format: "YYYY_Q", e.g., "2022_4")
 *   icc_value = ICC value (optional)
 *   max_rate  = Maximum rate for y-axis (optional, default = 0.5)
 *   smooth    = LOESS smoothing parameter (optional, default = 0.35)
 ************************************************/

%macro funnel_plot(data=, quarter=, icc_value=, max_rate=0.5, smooth=0.35);

  %put NOTE: Starting funnel_plot macro...;
  %put NOTE: Data = &data.;
  %put NOTE: Quarter = &quarter.;
  %put NOTE: ICC = &icc_value.;

  /*------------------------------------------------------------------*/
  /* Parameter validation                                              */
  /*------------------------------------------------------------------*/

  %if %sysevalf(%superq(data)=, boolean) %then %do;
    %put ERROR: The DATA parameter is required.;
    %return;
  %end;

  %if %sysevalf(%superq(quarter)=, boolean) %then %do;
    %put ERROR: The QUARTER parameter is required.;
    %return;
  %end;

  /* Check if dataset exists */
  %if %sysfunc(exist(&data.)) = 0 %then %do;
    %put ERROR: Dataset &data. does not exist.;
    %return;
  %end;

  /*------------------------------------------------------------------*/
  /* Step 1: Calculate hospital-specific rates                        */
  /*------------------------------------------------------------------*/

  /* Calculate sum of events by quarter and hospital */
  proc summary data=&data. nway;
    class an_trim eta;
    var event;
    output out=base_eta
      sum=event_sum;
  run;

  /* Calculate observed rate per hospital */
  data base_eta;
    set base_eta;
    n_patients = _FREQ_;
    tx_obs = event_sum / n_patients;
  run;

  /*------------------------------------------------------------------*/
  /* Step 2: Calculate overall (global) rate                          */
  /*------------------------------------------------------------------*/

  /* Calculate sum of events by quarter */
  proc summary data=&data. nway;
    class an_trim;
    var event;
    output out=base_global
      sum=event_sum_global;
  run;

  /* Calculate overall observed rate */
  data base_global;
    set base_global;
    n_global = _FREQ_;
    tx_global = event_sum_global / n_global;
  run;

  /*------------------------------------------------------------------*/
  /* Step 3: Merge hospital and global data                           */
  /*------------------------------------------------------------------*/

  /* Sort datasets by quarter for merging */
  proc sort data=base_eta; by an_trim; run;
  proc sort data=base_global; by an_trim; run;

  /*------------------------------------------------------------------*/
  /* Step 4: Calculate binomial control limits                        */
  /*------------------------------------------------------------------*/

  /* Merge hospital-specific and global data */
  data base_funnel;
    merge base_eta base_global;
    by an_trim;
    
    /* ICC retrieved from GEE model or provided as parameter */
    %if %sysevalf(%superq(icc_value)=, boolean) %then %do;
      /* If ICC not provided, use a default value */
      ICC = 0.001;
      %put WARNING: ICC value not provided. Using default ICC = 0.001;
    %end;
    %else %do;
      ICC = &icc_value.;
    %end;
    
    /* Inflation factor to account for clustering effect */
    IF_rate = 1 + (n_patients - 1) * ICC;
    IF_rate = min(IF_rate, 2);
    
    /* Define confidence levels */
    alpha_2SD = 0.05;    /* 2 standard deviations (95% confidence interval) */
    alpha_3SD = 0.0027;  /* 3 standard deviations (99.73% confidence interval) */
    
    /* Calculate quantiles from binomial distribution */
    /* Upper control limits (UCL) */
    q_upper_2SD = quantile('BINOMIAL', 1 - alpha_2SD/2, tx_global, n_patients);
    q_upper_3SD = quantile('BINOMIAL', 1 - alpha_3SD/2, tx_global, n_patients);
    
    /* Lower control limits (LCL) */
    q_lower_2SD = quantile('BINOMIAL', alpha_2SD/2, tx_global, n_patients);
    q_lower_3SD = quantile('BINOMIAL', alpha_3SD/2, tx_global, n_patients);
    
    /* Calculate deviations from global rate */
    diff_upper_2SD = (q_upper_2SD / n_patients) - tx_global;
    diff_upper_3SD = (q_upper_3SD / n_patients) - tx_global;
    
    diff_lower_2SD = tx_global - (q_lower_2SD / n_patients);
    diff_lower_3SD = tx_global - (q_lower_3SD / n_patients);
    
    /* Apply inflation factor (ICC adjustment) */
    UCL_2SD = min(1, tx_global + diff_upper_2SD * sqrt(IF_rate));
    UCL_3SD = min(1, tx_global + diff_upper_3SD * sqrt(IF_rate));
    
    LCL_2SD = max(0, tx_global - diff_lower_2SD * sqrt(IF_rate));
    LCL_3SD = max(0, tx_global - diff_lower_3SD * sqrt(IF_rate));
    
    /* Flag to identify real observations */
    real_obs = 1;
  run;

  /* Check if the specified quarter exists */
  proc sql noprint;
    select count(*) into :quarter_exists
    from base_funnel
    where an_trim = "&quarter.";
  quit;

  %if &quarter_exists. = 0 %then %do;
    %put ERROR: Quarter &quarter. does not exist in the dataset.;
    %put NOTE: Available quarters:;
    proc sql;
      select distinct an_trim as Available_Quarters
      from base_funnel;
    quit;
    %return;
  %end;

  /*------------------------------------------------------------------*/
  /* Step 5: Create complete dataset for smooth curves                */
  /*------------------------------------------------------------------*/

  /* Retrieve maximum number of patients and global rate for the quarter */
  proc sql noprint;
    select max(n_patients) into :max_n 
    from base_funnel 
    where an_trim="&quarter.";
    
    select tx_global into :tx_global_trim 
    from base_funnel 
    where an_trim="&quarter.";
  quit;

  /* Create sequence of values for the x-axis */
  data seq_n;
    do n_patients = 1 to &max_n.;
      output;
    end;
  run;

  /* Calculate control limits for all possible n_patients values */
  data base_funnel_complet;
    set seq_n;
    
    /* Global rate for the quarter */
    tx_global = &tx_global_trim.;
    
    /* ICC retrieved from GEE model or provided as parameter */
    %if %sysevalf(%superq(icc_value)=, boolean) %then %do;
      ICC = 0.001;
    %end;
    %else %do;
      ICC = &icc_value.;
    %end;
    
    /* Inflation factor to account for clustering effect */
    IF_rate = 1 + (n_patients - 1) * ICC;
    IF_rate = min(IF_rate, 2);
    
    /* Define confidence levels */
    alpha_2SD = 0.05;
    alpha_3SD = 0.0027;
    
    /* Calculate quantiles from binomial distribution */
    q_upper_2SD = quantile('BINOMIAL', 1 - alpha_2SD/2, tx_global, n_patients);
    q_upper_3SD = quantile('BINOMIAL', 1 - alpha_3SD/2, tx_global, n_patients);
    
    q_lower_2SD = quantile('BINOMIAL', alpha_2SD/2, tx_global, n_patients);
    q_lower_3SD = quantile('BINOMIAL', alpha_3SD/2, tx_global, n_patients);
    
    /* Calculate deviations from global rate */
    diff_upper_2SD = (q_upper_2SD / n_patients) - tx_global;
    diff_upper_3SD = (q_upper_3SD / n_patients) - tx_global;
    
    diff_lower_2SD = tx_global - (q_lower_2SD / n_patients);
    diff_lower_3SD = tx_global - (q_lower_3SD / n_patients);
    
    /* Apply inflation factor (ICC adjustment) */
    UCL_2SD = min(1, tx_global + diff_upper_2SD * sqrt(IF_rate));
    UCL_3SD = min(1, tx_global + diff_upper_3SD * sqrt(IF_rate));
    
    LCL_2SD = max(0, tx_global - diff_lower_2SD * sqrt(IF_rate));
    LCL_3SD = max(0, tx_global - diff_lower_3SD * sqrt(IF_rate));
    
    an_trim = "&quarter.";
    real_obs = 0;
  run;

  /*------------------------------------------------------------------*/
  /* Step 6: Filter data for the target quarter                       */
  /*------------------------------------------------------------------*/

  /* Keep only observations for the selected quarter */
  data base_funnel_v5; 
    set base_funnel; 
    where an_trim="&quarter."; 
  run;

  /*------------------------------------------------------------------*/
  /* Step 7: Merge datasets for plotting                              */
  /*------------------------------------------------------------------*/

  /* Combine observed data with complete smoothed curves dataset */
  data base_funnel_plot;
    set base_funnel_v5 base_funnel_complet;
  run;

  /*------------------------------------------------------------------*/
  /* Step 8: Generate funnel plot                                     */
  /*------------------------------------------------------------------*/

  /* Create funnel plot with LOESS-smoothed control limits */
  proc sgplot data=base_funnel_plot noborder noautolegend;
    /* Plot hospital points with labels */
    scatter x=n_patients y=tx_obs / 
      markerattrs=(size=8 color=black)
      datalabel=eta
      datalabelattrs=(size=8);
    
    /* Central line: global rate */
    series x=n_patients y=tx_global / 
      lineattrs=(color=blue thickness=1.5);
    
    /* 3 standard deviation limits with LOESS smoothing */
    loess x=n_patients y=UCL_3SD / 
      smooth=&smooth. 
      lineattrs=(color=red thickness=1.5) 
      nomarkers;
    
    loess x=n_patients y=LCL_3SD / 
      smooth=&smooth. 
      lineattrs=(color=green thickness=1.5) 
      nomarkers;
    
    /* 2 standard deviation limits with LOESS smoothing */
    loess x=n_patients y=UCL_2SD / 
      smooth=&smooth. 
      lineattrs=(color=lightred thickness=1.5 pattern=dash) 
      nomarkers;
    
    loess x=n_patients y=LCL_2SD / 
      smooth=&smooth. 
      lineattrs=(color=lightgreen thickness=1.5 pattern=dash) 
      nomarkers;
    
    /* Legend */
    keylegend "2DS" "3DS" / location=inside position=topright across=2;
    
    /* Axis formatting */
    yaxis label="Event rate per hospital (%)" 
          min=0 max=&max_rate. 
          valuesformat=percent9.0;
    xaxis label="Hospital volume" 
          min=0 max=&max_n.;
    
    /* Title */
    title "Funnel Plot - Quarter &quarter.";
    %if %sysevalf(%superq(icc_value)=, boolean) %then %do;
      *title2 "ICC = 0.001 (default)";
    %end;
    %else %do;
      *title2 "ICC = &icc_value.";
    %end;
  run;

  /*------------------------------------------------------------------*/
  /* Step 9: Clean up temporary datasets                              */
  /*------------------------------------------------------------------*/

  proc datasets lib=work nolist;
    delete base_eta base_global base_funnel base_funnel_complet 
           base_funnel_v5 base_funnel_plot seq_n;
  quit;

  %put NOTE: Funnel plot for quarter &quarter. completed successfully.;

%mend funnel_plot;
