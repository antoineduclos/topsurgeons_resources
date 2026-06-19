/************************************************
 * MACRO: Shewhart control_chart
 * Aim = To create control chart for a specific hospital
 *          with binomial control limits and ICC adjustment
 * 
 * Parameters:
 *   data            = Input dataset with predictions (must contain: an_trim, eta, event, pred, date)
 *   hospital        = Hospital identifier to plot (eta value)
 *   start_year      = Start year for the plot (e.g., 2020)
 *   end_year        = End year for the plot (e.g., 2024)
 *   icc_value       = ICC value from GEE model (optional)
 *   max_rate        = Maximum rate for y-axis (optional, default = 0.4)
 *   min_rate        = Minimum rate for y-axis (optional, default = 0)
 *   y_step          = Step for y-axis values (optional, default = 0.05)
 *   y_values        = Custom y-axis values (optional, overrides min_rate, max_rate, y_step)
 *   exclude_quarter = Quarter to exclude from plot (optional, e.g., "2025_1")
 *   auto_scale      = Automatically adjust y-axis limits (optional, default = 1)
 *   margin_top      = Top margin for auto-scaling (optional, default = 0.2)
 ************************************************/

%macro control_chart(data=, hospital=, start_year=, end_year=, 
                     icc_value=, max_rate=0.4, min_rate=0, y_step=0.05,
                     y_values=, exclude_quarter=, auto_scale=1, margin_top=0.2);

  %put NOTE: Starting control_chart macro...;
  %put NOTE: Data = &data.;
  %put NOTE: Hospital = &hospital.;
  %put NOTE: Start year = &start_year.;
  %put NOTE: End year = &end_year.;
  %put NOTE: ICC = &icc_value.;

  /*------------------------------------------------------------------*/
  /* Parameter validation                                              */
  /*------------------------------------------------------------------*/

  %if %sysevalf(%superq(data)=, boolean) %then %do;
    %put ERROR: The DATA parameter is required.;
    %return;
  %end;

  %if %sysevalf(%superq(hospital)=, boolean) %then %do;
    %put ERROR: The HOSPITAL parameter is required.;
    %return;
  %end;

  %if %sysevalf(%superq(start_year)=, boolean) %then %do;
    %put ERROR: The START_YEAR parameter is required.;
    %return;
  %end;

  %if %sysevalf(%superq(end_year)=, boolean) %then %do;
    %put ERROR: The END_YEAR parameter is required.;
    %return;
  %end;

  /* Check if dataset exists */
  %if %sysfunc(exist(&data.)) = 0 %then %do;
    %put ERROR: Dataset &data. does not exist.;
    %return;
  %end;

  /*------------------------------------------------------------------*/
  /* Step 1: Calculate overall period rates by hospital              */
  /*------------------------------------------------------------------*/

  /* Sort data by hospital */
  proc sort data=&data.; by eta; run;

  /* Calculate sum of observed and predicted events over the whole period */
  proc means data=&data. noprint; 
    output out=base_chir_carte_global 
           sum=obs_evenement pred_evenement 
           n=nb_pat_global;
    var event pred; 
    by eta; 
  run;

  /* Calculate observed and predicted global rates */
  data base_chir_carte_global; 
    set base_chir_carte_global;
    tx_obs_global = obs_evenement / nb_pat_global;
    tx_pred_global = pred_evenement / nb_pat_global;
    Calibr = tx_obs_global - tx_pred_global;
    keep eta tx_obs_global tx_pred_global Calibr nb_pat_global;
  run;

  /*------------------------------------------------------------------*/
  /* Step 2: Calculate quarterly rates by hospital                   */
  /*------------------------------------------------------------------*/

  /* Sort data by hospital and quarter */
  proc sort data=&data.; by eta an_trim; run;

  /* Calculate sum of observed and predicted events by quarter */
  proc means data=&data. noprint; 
    output out=base_chir_carte_trim0 
           sum=obs_evenement pred_evenement 
           n=nb_pat_trim;
    var event pred; 
    by eta an_trim; 
  run;

  /* Calculate observed and predicted quarterly rates */
  data base_chir_carte_trim1; 
    set base_chir_carte_trim0;
    tx_obs_trim = obs_evenement / nb_pat_trim;
    tx_pred_trim = pred_evenement / nb_pat_trim;
    drop _type_ _freq_;
  run;

  /*------------------------------------------------------------------*/
  /* Step 3: Merge global and quarterly data                          */
  /*------------------------------------------------------------------*/

  /* Sort datasets by hospital for merging */
  proc sort data=base_chir_carte_trim1; by eta; run;
  proc sort data=base_chir_carte_global; by eta; run;

  /* Merge global calibration with quarterly data */
  data base_chir_carte_v2; 
    merge base_chir_carte_trim1 base_chir_carte_global; 
    by eta;
    
    /* ICC from GEE model (provided as parameter or using default) */
    %if %sysevalf(%superq(icc_value)=, boolean) %then %do;
      ICC_evenement = 0.001;
      %put WARNING: ICC value not provided. Using default ICC = 0.001;
    %end;
    %else %do;
      ICC_evenement = &icc_value.;
    %end;
    
    /* Inflation factor to account for clustering */
    IF_rate = 1 + (nb_pat_trim - 1) * ICC_evenement;
    
    /* Cap inflation factor at 2 (sqrt(2)=1.41) when predicted rate < 20% */
    if (IF_rate > 2 and tx_pred_trim < 0.20) then IF_rate = 2;
    
    /* Define confidence levels */
    alpha3SD = 0.002699796;  /* 3SD = 99.7300204% */
    alpha2SD = 0.045500264;  /* 2SD = 95.4499736% */
    
    /* Calculate control limits with calibration */
    /* 3SD limits */
    UCL3ds = Calibr + tx_pred_trim + 
             abs((quantile('BINOMIAL', 1-alpha3SD, tx_pred_trim, nb_pat_trim)/nb_pat_trim) - tx_pred_trim) * sqrt(IF_rate);
    LCL3ds = Calibr + tx_pred_trim - 
             abs(tx_pred_trim - (quantile('BINOMIAL', alpha3SD, tx_pred_trim, nb_pat_trim)/nb_pat_trim)) * sqrt(IF_rate);
    
    /* 2SD limits */
    UCL2ds = Calibr + tx_pred_trim + 
             abs((quantile('BINOMIAL', 1-alpha2SD, tx_pred_trim, nb_pat_trim)/nb_pat_trim) - tx_pred_trim) * sqrt(IF_rate);
    LCL2ds = Calibr + tx_pred_trim - 
             abs(tx_pred_trim - (quantile('BINOMIAL', alpha2SD, tx_pred_trim, nb_pat_trim)/nb_pat_trim)) * sqrt(IF_rate);
    
    /* Calibrated central line */
    tx_pred_trim_calib = tx_pred_trim + Calibr;
    
    /* Cap limits to [0,1] range */
    if tx_pred_trim_calib > 1.0 then tx_pred_trim_calib = 1.0;
    if tx_pred_trim_calib < 0.0 then tx_pred_trim_calib = 0.0;
    
    if UCL3ds > 1.0 then UCL3ds = 1.0;
    if LCL3ds < 0.0 then LCL3ds = 0.0;
    if UCL2ds > 1.0 then UCL2ds = 1.0;
    if LCL2ds < 0.0 then LCL2ds = 0.0;
    
    /* Keep only useful variables */
    keep eta an_trim nb_pat_trim 
         tx_obs_trim tx_pred_trim tx_pred_trim_calib Calibr
         UCL3ds LCL3ds UCL2ds LCL2ds;
  run;

  /*------------------------------------------------------------------*/
  /* Step 4: Add year information from original data                 */
  /*------------------------------------------------------------------*/

  /* Extract year from date variable in the original dataset */
  proc sort data=&data.; by eta an_trim; run;

  /* Get unique year for each quarter */
  data base_chir_years;
    set &data.;
    year = year(date);
    keep eta an_trim year;
  run;

  /* Remove duplicates (one year per quarter) */
  proc sort data=base_chir_years nodupkey; by eta an_trim; run;

  /* Merge year information with control chart data */
  proc sort data=base_chir_carte_v2; by eta an_trim; run;
  proc sort data=base_chir_years; by eta an_trim; run;

  data base_chir_carte_v2_with_year;
    merge base_chir_carte_v2 base_chir_years;
    by eta an_trim;
  run;

  /*------------------------------------------------------------------*/
  /* Step 5: Filter data for the target hospital and year range      */
  /*------------------------------------------------------------------*/

  /* Keep only observations for the selected hospital and year range */
  data base_chir_carte_v2_5; 
    set base_chir_carte_v2_with_year; 
    where eta = "&hospital." 
      and year >= &start_year. 
      and year <= &end_year.;
    
    /* Exclude specific quarter if provided */
    %if %sysevalf(%superq(exclude_quarter) ne, boolean) %then %do;
      and an_trim ^= "&exclude_quarter.";
    %end;
  run;

  /* Check if there are data points after filtering */
  proc sql noprint;
    select count(*) into :n_obs
    from base_chir_carte_v2_5;
  quit;

  %if &n_obs. = 0 %then %do;
    %put ERROR: No data found for hospital &hospital. in years &start_year.-&end_year.;
    %return;
  %end;

  /*------------------------------------------------------------------*/
  /* Step 6: Calculate dynamic y-axis limits if auto_scale is enabled */
  /*------------------------------------------------------------------*/

  %if &auto_scale. = 1 and %sysevalf(%superq(y_values)=, boolean) %then %do;

    /* Calculate max values for dynamic scaling */
    proc sql noprint;
      select max(tx_obs_trim, UCL3ds, tx_pred_trim_calib) into :max_observed
      from base_chir_carte_v2_5;
      
      select min(LCL3ds, tx_obs_trim, tx_pred_trim_calib) into :min_observed
      from base_chir_carte_v2_5;
    quit;

    /* Calculate dynamic max with margin */
    %let dynamic_max = %sysevalf(&max_observed. * (1 + &margin_top.));
    
    /* Round up to the nearest multiple of y_step */
    %let dynamic_max = %sysevalf(ceil(&dynamic_max. / &y_step.) * &y_step.);
    
    /* Ensure min is at least 0 */
    %if &min_observed. < 0 %then %do;
      %let dynamic_min = %sysevalf(&min_observed. - 0.01);
      %let dynamic_min = %sysevalf(floor(&dynamic_min. / &y_step.) * &y_step.);
    %end;
    %else %do;
      %let dynamic_min = 0;
    %end;

    /* Override max_rate and min_rate with dynamic values */
    %let max_rate = &dynamic_max.;
    %let min_rate = &dynamic_min.;
    
    %put NOTE: Auto-scaling enabled;
    %put NOTE:   Max observed = &max_observed.;
    %put NOTE:   Min observed = &min_observed.;
    %put NOTE:   Dynamic max   = &max_rate.;
    %put NOTE:   Dynamic min   = &min_rate.;

  %end;

  /*------------------------------------------------------------------*/
  /* Step 7: Generate x-axis values with quotes for SGPLOT           */
  /*------------------------------------------------------------------*/

  /* Create list of quarters for the x-axis with quotes */
  proc sql noprint;
    select distinct quote(strip(an_trim)) into :quarter_values separated by ' '
    from base_chir_carte_v2_5
    order by an_trim;
    
    /* Also create display labels */
    select distinct an_trim,
           case when mod(input(scan(an_trim, 2, '_'), 2.), 4) = 1 
                then cats(scan(an_trim, 1, '_'), ' (Q1)')
                else ''
           end into :value_list separated by ' ', :label_list separated by ' '
    from base_chir_carte_v2_5
    order by an_trim;
  quit;

  %put NOTE: Quarter values with quotes: &quarter_values.;

  /*------------------------------------------------------------------*/
  /* Step 8: Generate control chart                                  */
  /*------------------------------------------------------------------*/

  /* Create the control chart using SGPLOT */
  proc sgplot data=base_chir_carte_v2_5 noautolegend noborder;
    /* Central line: calibrated predicted rate */
    series x=an_trim y=tx_pred_trim_calib / 
      lineattrs=(color=blue pattern=solid thickness=1.5)
      legendlabel="Calibrated predicted rate" name="center";
    
    /* 3SD control limits */
    series x=an_trim y=UCL3ds / 
      lineattrs=(color=red pattern=solid thickness=1.5)
      legendlabel="UCL 3SD" name="ucl3";
    series x=an_trim y=LCL3ds / 
      lineattrs=(color=green pattern=solid thickness=1.5)
      legendlabel="LCL 3SD" name="lcl3";
    
    /* 2SD control limits */
    series x=an_trim y=UCL2ds / 
      lineattrs=(color=crimson pattern=solid thickness=1)
      legendlabel="UCL 2SD" name="ucl2";
    series x=an_trim y=LCL2ds / 
      lineattrs=(color=darkgreen pattern=solid thickness=1)
      legendlabel="LCL 2SD" name="lcl2";
    
    /* Observed rates with markers and line */
    series x=an_trim y=tx_obs_trim / 
      MARKERS 
      markerattrs=(size=7 symbol=circleFilled color=black) 
      lineattrs=(THICKNESS=1.5 color=black)
      legendlabel="Observed rate" name="obs";
    
    /* X-axis formatting with dynamic values */
    xaxis label="Quarter over years" 
          offsetmax=0.05 offsetmin=0.05
          values=(&quarter_values.);
    
    /* Y-axis formatting with custom or dynamic values */
    %if %sysevalf(%superq(y_values)=, boolean) %then %do;
      /* Use min, max and step if y_values not provided */
      %if &auto_scale. = 1 %then %do;
        /* Auto-scaling: use dynamic min and max */
        yaxis label="Rate" 
              offsetmax=0.05 offsetmin=0.05
              min=&min_rate. max=&max_rate. 
              valuesformat=percent9.0
              values=(&min_rate. to &max_rate. by &y_step.);
      %end;
      %else %do;
        /* Manual scaling with user-provided values */
        yaxis label="Rate" 
              offsetmax=0.05 offsetmin=0.05
              min=&min_rate. max=&max_rate. 
              valuesformat=percent9.0
              values=(&min_rate. to &max_rate. by &y_step.);
      %end;
    %end;
    %else %do;
      /* Use custom y_values */
      yaxis label="Rate" 
            offsetmax=0.05 offsetmin=0.05
            valuesformat=percent9.0
            values=(&y_values.);
    %end;
    
    /* Legend */
    keylegend "center" "ucl3" "lcl3" "ucl2" "lcl2" "obs" / 
      location=inside position=topright across=2;
    
    /* Title */
    title "Control Chart - Hospital &hospital.";
    title2 "Years &start_year.-&end_year.";
    %if %sysevalf(%superq(icc_value)=, boolean) %then %do;
      *title3 "ICC = 0.001 (default)";
    %end;
    %else %do;
      *title3 "ICC = &icc_value.";
    %end;
  run;

  /*------------------------------------------------------------------*/
  /* Step 9: Clean up temporary datasets                              */
  /*------------------------------------------------------------------*/

  proc datasets lib=work nolist;
    delete base_chir_carte_global base_chir_carte_trim0 
           base_chir_carte_trim1 base_chir_carte_v2 
           base_chir_carte_v2_with_year base_chir_carte_v2_5
           base_chir_years;
  quit;

  %put NOTE: Control chart for hospital &hospital. completed successfully.;

%mend control_chart;
