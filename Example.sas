/************************************************
 * PROGRAM: Example with simulated data
 * AIM = Create funnel plot and control charts with ICC adjustment
 ************************************************/

/*------------------------------------------------------------------*/
/* Section 1: Library and Data Preparation                          */
/*------------------------------------------------------------------*/

/* Libname */
libname lib "D:/Git ressources";

/* Calculate quarter and year from date variable */
data base_chir;
    set lib.dataset;
    year = year(date);
    trimester = qtr(date);
    format date date9.;
    an_trim = catx("_", year, trimester);  /* Create year_quarter identifier */
run;

/*------------------------------------------------------------------*/
/* Section 2: Data Exploration                                      */
/*------------------------------------------------------------------*/

/* Display dataset structure and variable attributes */
proc contents data=base_chir; run;

/* Frequency distributions for key variables */
proc freq data=base_chir; 
    table age cha eta sex urg year trimester event; 
run;

/*------------------------------------------------------------------*/
/* Section 3: GEE Regression Model                                  */
/*------------------------------------------------------------------*/

/* Fit GEE model with exchangeable correlation structure */
proc genmod data=base_chir descending ;
    /* Define categorical variables with reference categories */
    class age (ref="18-59") 
          cha (ref="Charlson 0-1") 
          eta 
          sex (ref="F") 
          urg (ref="0") 
          event (ref="1") ;
    
    /* Specify the logistic regression model */
    model event = age cha sex urg / link=logit dist=binomial;
    
    /* Specify exchangeable correlation structure with hospital as subject */
    repeated subject=eta / type=exch;

    /* Estimate odds ratios for age groups */
    estimate "Age 60-74 vs 18-59" age 1 0 -1 / exp;
    estimate "Age 75+ vs 18-59" age 0 1 -1 / exp;
    
    /* Estimate odds ratios for sex and urgency */
    estimate "Female vs Male" sex 1 -1 / exp;
    estimate "Urgent vs non-urgent" urg 1 -1 / exp;
    
    /* Estimate odds ratios for Charlson comorbidity index */
    estimate "Charlson 2-3 vs 0-1" cha 1 0 -1 / exp;
    estimate "Charlson 4+ vs 0-1" cha 0 1 -1 / exp;

    /* Output correlation matrix and predictions */
    ods output GEEExchCorr = ICC_classique; 
    output out = base_chir_pred pred = pred;
    ods output estimates = Estimates_classique;
    store modeleGEE_classique;
run;

/*------------------------------------------------------------------*/
/* Section 4: Extract ICC for Funnel Plot                           */
/*------------------------------------------------------------------*/

/* First, look at the structure of ICC_classique dataset */
proc contents data=ICC_classique; run;
proc print data=ICC_classique; run;

/* Extract ICC - Version robuste */
data _null_;
    set ICC_classique;
    
    /* Afficher toutes les variables dans le log pour debug */
    put _all_;
    
    /* Try different possible variable names */
    if not missing(exchangeablecorr) then do;
        call symputx('ICC_value', exchangeablecorr);
        call symputx('ICC_exists', 'YES');
        put 'NOTE: ICC extracted from variable "exchangeablecorr" = ' exchangeablecorr;
    end;
    else if not missing(ExchCorr) then do;
        call symputx('ICC_value', ExchCorr);
        call symputx('ICC_exists', 'YES');
        put 'NOTE: ICC extracted from variable "ExchCorr" = ' ExchCorr;
    end;
    else if not missing(Corr) then do;
        call symputx('ICC_value', Corr);
        call symputx('ICC_exists', 'YES');
        put 'NOTE: ICC extracted from variable "Corr" = ' Corr;
    end;
    else do;
        /* Si aucune variable trouvée, utiliser une valeur par défaut */
        call symputx('ICC_value', 0.001);
        call symputx('ICC_exists', 'NO');
        put 'WARNING: No correlation variable found. Using default ICC = 0.001';
    end;
run;

/* Display ICC values in the log */
%put ===========================================;
%put ICC extraction results:;
%put   ICC_VALUE  = &ICC_value.;
%put   ICC_EXISTS = &ICC_exists.;

/* Si ICC_VALUE est vide ou égal à ., forcer une valeur par défaut */
%if &ICC_value. = . or &ICC_value. = %then %do;
    %let ICC_value = 0.001;
    %put WARNING: ICC_VALUE was empty. Forcing default ICC = 0.001;
%end;

/* Take absolute value for the funnel plot */
%let ICC_for_plot = %sysfunc(abs(&ICC_value.));
%put   ICC_for_plot = &ICC_for_plot.;
%put ===========================================;

/*====================================================================*/
/* INCLUDE MACROS                                                     */
/*====================================================================*/

%include "D:/Git ressources/funnel_plot_macro.sas";
%include "D:/Git ressources/control_chart_macro.sas";

/*====================================================================*/
/* GENERATE PLOTS                                                     */
/*====================================================================*/

/* Example 1: Funnel plot for 2022 T4 */
%funnel_plot(data=base_chir_pred, quarter=2022_4, icc_value=&ICC_for_plot.);

/* Example 2: Shewhart control chart for hospital 7 */
%control_chart(data=base_chir_pred, 
               hospital=etab07, 
               start_year=2020, 
               end_year=2024, 
               icc_value=&ICC_for_plot.,
               y_values=0 to 0.40 by 0.05);
