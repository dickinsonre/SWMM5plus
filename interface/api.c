#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "swmm5.h"
#include "headers.h"
#include "api_error.h"
#include "api.h"

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
#define REAL4 float

extern REAL4* NodeResults;             //  "
extern REAL4* LinkResults;             //  "
struct stat st = {0};

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------

Interface* api;
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input
char errmsg[1024];

//===============================================================================
// --- Simulation
//===============================================================================

// Initializes the EPA-SWMM simulation. It creates an Interface
// variable (details about Interface in interface.h). The
// function opens de SWMM input file and creates report and
// output files. Although the .inp is parsed within swmm_start,
// not every property is extracted, e.g., slopes of trapezoidal channels.
// In swmm.c
// to parse the .inp again using EPA-SWMM functionalities
// within api_load_vars
// The interface object is passed to many functions in interface.c
// but it is passed as a void pointer. This is because the object
// is also used by the Fortran engine. Treating Interface as void*
// facilitates interoperability

//===============================================================================
int DLLEXPORT api_initialize(
    char* f1, char* f2, char* f3, int run_routing)
//===============================================================================
{
    int error;
    api = (Interface*) malloc(sizeof(Interface));
    error = swmm_open(f1, f2, f3);
    if (error) {
        api->IsInitialized = FALSE;
        return 0;
    }
    error = swmm_start(0);
    if (error) {
        api->IsInitialized = FALSE;
        return 0;
    }
    api->IsInitialized = TRUE;
    error = api_load_vars();
    if (!run_routing) IgnoreRouting = TRUE;
    return 0;
}

//===============================================================================
int DLLEXPORT api_finalize()
//===============================================================================
    //
    //  Input: f_api is an Interface object passed as a void*
    //  Output: None
    //  Purpose: Closes the link with the EPA-SWMM library
{
    int ii;

     //printf("\n in api finalize \n");

    swmm_end();

    //printf(" ------------------\n");
    //printf("Fout.mode %d, %d %d \n",Fout.mode,SCRATCH_FILE,SAVE_FILE);
    //printf(" ------------------\n");

    // NOTE -- because we always initialize with
    // an *.out file name, the Fout.mode is SAVE_FILE
    // so the swmm_report() is NOT called at the end
    //if ( Fout.mode == SCRATCH_FILE ) swmm_report();

    swmm_close();

    // frees double variables in API
    for (ii = 0; ii < NUM_API_DOUBLE_VARS; ii++)
    {
        if (api->double_vars[ii] != NULL)
            free(api->double_vars[ii]);
    }

    // // frees integer variables in API --- these do not exist brh20211217
    // for (i = 0; i < NUM_API_INT_VARS; i++)
    // {
    //     if (api->int_vars[i] != NULL)
    //         free(api->int_vars[i]);
    // }

    free(api);

    return 0;
}

//===============================================================================
double DLLEXPORT api_run_step()
//===============================================================================
{
    swmm_step(&(api->elapsedTime));

    return api->elapsedTime;
}

//===============================================================================
// --- Property-extraction
// * During Simulation
//===============================================================================

//===============================================================================
int DLLEXPORT api_get_node_results(
    char* node_name, float* inflow, float* overflow, float* depth, float* volume)
//===============================================================================
    //
    //  Input:    f_api = Interface object passed as a void*
    //            node_name = string identifier of node
    //            inflow, overflow, depth, volume =
    //  Output: None
    //  Purpose: Closes the link with the SWMM C library
    //
{
    int jj, error;

    error = check_api_is_initialized("api_get_node_results");
    if (error != 0) return error;

    jj = project_findObject(NODE, node_name);

    *inflow   = Node[jj].inflow;
    *overflow = Node[jj].overflow;
    *depth    = Node[jj].newDepth;
    *volume   = Node[jj].newVolume;

    return 0;
}

//===============================================================================
int DLLEXPORT api_get_link_results(
    char* link_name, float* flow, float* depth, float* volume)
//===============================================================================
{
    int jj, error;

    error = check_api_is_initialized("api_get_link_results");
    if (error != 0) return error;

    jj = project_findObject(LINK, link_name);

    *flow   = Link[jj].newFlow;
    *depth  = Link[jj].newDepth;
    *volume = Link[jj].newVolume;

    return 0;
}

//===============================================================================
// --- Property-extraction
// * After Initialization
//===============================================================================

//===============================================================================
double DLLEXPORT api_get_start_datetime()
//===============================================================================
{
    return StartDateTime;
}

//===============================================================================
double DLLEXPORT api_get_end_datetime()
//===============================================================================
{
    return EndDateTime;
}

//===============================================================================
double DLLEXPORT api_get_TotalDuration()
//===============================================================================
// this is the SWMM-C Total Duration in milliseconds
{
    return TotalDuration;
}

//===============================================================================
int DLLEXPORT api_get_flowBC(
    int node_idx, double current_datetime, double* flowBC)
//===============================================================================
{
    int ptype, pcount, ii, jj, pp;
    int tyear, tmonth, tday;
    int thour, tmin, tsec, tweekday, attr;
    double val;
    double bline, sfactor;
    double total_factor = 1;
    double total_extinflow = 0;

    *flowBC = 0;
    datetime_decodeDate(current_datetime, &tyear, &tmonth, &tday);
    datetime_decodeTime(current_datetime, &thour, &tmin,   &tsec);
    tweekday = datetime_dayOfWeek(current_datetime);

    // handle dry weather inflows
    attr = nodef_has_dwfInflow;
    api_get_nodef_attribute(node_idx, attr, &val);
    if (val == 1) { // node_has_dwfInflow 
        for(ii=0; ii<4; ii++)
        {
            // jj is the pattern #
            jj = Node[node_idx].dwfInflow->patterns[ii];
            if (jj > 0)
            {
                ptype = Pattern[jj].type;
                // if (ptype == MONTHLY_PATTERN)
                //     total_factor *= Pattern[jj].factor[tmonth-1];
                // else if (ptype == DAILY_PATTERN)
                //     total_factor *= Pattern[jj].factor[tweekday-1];
                // else if (ptype == HOURLY_PATTERN)
                //     total_factor *= Pattern[jj].factor[thour];
                // else if (ptype == WEEKEND_PATTERN)
                // {
                //     if ((tweekday == 1) || (tweekday == 7))
                //         total_factor *= Pattern[jj].factor[thour];
                // }
                switch(ptype) {
                    case MONTHLY_PATTERN :
                        total_factor *= Pattern[jj].factor[tmonth-1];
                        break;
                    case DAILY_PATTERN :
                        total_factor *= Pattern[jj].factor[tweekday-1];
                        break;
                    case HOURLY_PATTERN :
                        total_factor *= Pattern[jj].factor[thour];
                        break;
                    case WEEKEND_PATTERN :
                        if ((tweekday == 1) || (tweekday == 7))
                            total_factor *= Pattern[jj].factor[thour];
                        break;
                    default :
                        printf(" unexpected default reached in api_get_flowBC at A -- needs debugging");
                        return -1;
                }
            }
        }
        *flowBC += total_factor * CFTOCM(Node[node_idx].dwfInflow->avgValue);
    }

    // handle external inflows
    attr = nodef_has_extInflow;
    api_get_nodef_attribute(node_idx, attr, &val);
    //printf("\n here node_idx %d \n",node_idx);
    //printf(" val = %f \n ",val);
    //printf(" MONTHLY %d \n ", MONTHLY_PATTERN);
    if (val == 1) { // node_has_extInflow
        // pp is the pattern #
        pp = Node[node_idx].extInflow->basePat; // pattern
        //printf(" base pat %d \n ",pp);
        bline = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline); // baseline value
        //printf(" bline %f \n ",bline);
        if (pp >= 0)
        {
            ptype = Pattern[pp].type;
            // if (ptype == MONTHLY_PATTERN)
            // {
            //     //total_extinflow += Pattern[j].factor[mm-1] * bline;  //brh20211221
            //     total_extinflow += Pattern[pp].factor[tmonth-1] * bline;  //brh20211221
            // } 
            // else if (ptype == DAILY_PATTERN)
            //     //total_extinflow += Pattern[j].factor[dow-1] * bline; //brh20211221
            //     total_extinflow += Pattern[pp].factor[tweekday-1] * bline;//brh20211221
            // else if (ptype == HOURLY_PATTERN)
            //     //total_extinflow += Pattern[j].factor[h] * bline;  /brh20211221
            //     total_extinflow += Pattern[pp].factor[thour] * bline;   //brh20211221
            // else if (ptype == WEEKEND_PATTERN)
            // {
            //     if ((tweekday == 1) || (tweekday == 7))
            //         //total_extinflow += Pattern[j].factor[h] * bline;
            //         total_extinflow += Pattern[pp].factor[thour] * bline; //brh20211221
            // }
            switch(ptype) {
                case MONTHLY_PATTERN :
                    total_extinflow += Pattern[pp].factor[tmonth-1] * bline;
                    break;
                case DAILY_PATTERN :
                    total_extinflow += Pattern[pp].factor[tweekday-1] * bline;
                    break;
                case HOURLY_PATTERN :
                    total_extinflow += Pattern[pp].factor[thour] * bline;
                    break;
                case WEEKEND_PATTERN :
                    if ((tweekday == 1) || (tweekday == 7))
                        total_extinflow += Pattern[pp].factor[thour] * bline;
                    break;
                default :
                    printf(" unexpected default reached in api_get_flowBC at B -- needs debugging");
                    return -1;
            }
        }
        else if (bline > 0)  // no pattern, but baseline inflow provided
        {
            total_extinflow += bline;
        }

        // external inflow from time series are added to baseline and pattern
        // jj is the time series
        jj = Node[node_idx].extInflow->tSeries; // tseries
        sfactor = Node[node_idx].extInflow->sFactor; // scale factor
        if (jj >= 0)
        {
            total_extinflow += table_tseriesLookup(&Tseries[jj], current_datetime, FALSE) * sfactor;
        }

        // add the external inflows to the dry weather flows stored in flowBC
        *flowBC += total_extinflow;
    }
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_headBC(
    int node_idx, double current_datetime, double* headBC)
//===============================================================================
{
    int ii = Node[node_idx].subIndex;

//     printf("  in api_get_headBC with outfall %d \n",Outfall[ii].type);
//     printf("     node_idx = %d \n",node_idx);
//     printf("     subIdx   = %d \n",ii);
//     printf("     FIXED_OUTFALL = %d \n",FIXED_OUTFALL);
//     printf("     NORMAL_OUTFALL = %d \n",NORMAL_OUTFALL);
//     printf("     FREE_OUTFALL = %d \n",FREE_OUTFALL);
//     printf("     TIDAL_OUTFALL = %d \n",TIDAL_OUTFALL);
//    printf("     TIMESERIES_OUTFALL = %d \n",TIMESERIES_OUTFALL);
    
    switch (Outfall[ii].type)
    {
        case FIXED_OUTFALL:
            *headBC = FTTOM(Outfall[ii].fixedStage);
            return 0;

        // case NORMAL_OUTFALL:
        //     *headBC = API_NULL_VALUE_I;
        //     sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[ii].type, Node[node_idx].ID);
        //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
        //     printf("Error for outfall type NORMAL_OUTFALL \n");
        //     return api_err_not_developed;

        // case FREE_OUTFALL:
        //     *headBC = API_NULL_VALUE_I;
        //     sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[ii].type, Node[node_idx].ID);
        //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
        //     printf("Error for outfall type FREE_OUTFALL \n");
        //     return api_err_not_developed;

        // case TIDAL_OUTFALL:
        //     *headBC = API_NULL_VALUE_I;
        //     sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[ii].type, Node[node_idx].ID);
        //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
        //     printf("Error for outfall type TIDAL_OUTFALL \n");
        //     return api_err_not_developed;

        // case TIMESERIES_OUTFALL:    
        //     *headBC = API_NULL_VALUE_I;
        //     sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[ii].type, Node[node_idx].ID);
        //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
        //     printf("Error for outfall type TIMESERIES_OUTFALL \n");
        //     return api_err_not_developed;
        
        case FREE_OUTFALL:
            *headBC = FTTOM(Node[node_idx].initDepth + Node[node_idx].invertElev);
            return 0;

        case NORMAL_OUTFALL:
            *headBC = FTTOM(Node[node_idx].initDepth + Node[node_idx].invertElev);
            return 0;

        default:
            *headBC = API_NULL_VALUE_I;
            sprintf(errmsg, "OUTFALL TYPE %d at NODE %s", Outfall[ii].type, Node[node_idx].ID);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            printf("Unexpected default case for Outfall[ii].type \n");
            return api_err_not_developed;
    }
}

//===============================================================================
int DLLEXPORT api_get_SWMM_controls(
    int*  flow_units,
    int*  route_model,
    int*  allow_ponding,
    int*  inertial_damping,
    int*  num_threads,
    int*  skip_steady_state,
    int*  force_main_eqn,
    int*  max_trials,
    int*  normal_flow_limiter,
    int*  surcharge_method,
    int*  tempdir_provided,
    double* variable_step,
    double* lengthening_step,
    double* route_step,
    double* min_route_step,
    double* min_surface_area,
    double* min_slope,
    double* head_tol,
    double* sys_flow_tol,
    double* lat_flow_tol)
//===============================================================================
// Note, at this time this only gets the SWMM controls that are important to hydraulics
{
    int error;

    error = check_api_is_initialized("api_get_SWMM_controls");
    if (error) return error;

    *flow_units = FlowUnits;

    *route_model = RouteModel;

    *allow_ponding = AllowPonding;

    *inertial_damping = InertDamping;

    *num_threads = NumThreads;

    *skip_steady_state = SkipSteadyState;

    *force_main_eqn = ForceMainEqn;

    *max_trials = MaxTrials;

    *normal_flow_limiter = NormalFlowLtd;

    *surcharge_method = SurchargeMethod;

    *tempdir_provided = 0;
    if (strlen(TempDir) >0)
        *tempdir_provided = 1;
    
    *variable_step = CourantFactor;

    *lengthening_step = LengtheningStep;

    *route_step = RouteStep;
    
    *min_route_step = MinRouteStep;

    *min_surface_area = MinSurfArea;

    *min_slope = MinSlope;

    *head_tol = HeadTol;

    *sys_flow_tol = SysFlowTol;

    *lat_flow_tol = LatFlowTol;


    //printf(" RouteModel = %d \n",RouteModel);

    //printf(" \n tempdir len is %ld \n ",strlen(TempDir));
    //printf(" courant factor %f \n", CourantFactor);
    //printf("\n testing variable %s \n",TempDir);


    return 0;
}
//===============================================================================
int DLLEXPORT api_get_SWMM_times(
    double* starttime_epoch,
    double* endtime_epoch,
    double* report_start_datetime, 
    int*    report_step, 
    int*    hydrology_step, 
    int*    hydrology_dry_step, 
    double* hydraulic_step,
    double* total_duration) 
//===============================================================================    
{
    int error;

    error = check_api_is_initialized("api_get_SWMM_times");
    if (error) return error;

    *starttime_epoch       = StartDateTime;
    *endtime_epoch         = EndDateTime;
    *report_start_datetime = ReportStart;
    *report_step           = ReportStep;
    *hydrology_step        = WetStep;
    *hydrology_dry_step    = DryStep;
    *hydraulic_step        = RouteStep;
    *total_duration        = TotalDuration / 1000.0;

    // printf(" report start datetime = %f \n",ReportStart);
    // printf(" report step = %d \n", ReportStep);
    // printf(" hydrology_step = %d \n",WetStep);
    // printf(" hydrology_dry_step = %d \n",DryStep);
    // printf(" hydraulic_step = %f \n",RouteStep);

    return 0;
}

//===============================================================================
double DLLEXPORT api_get_NewRunoffTime()
//===============================================================================
{
    return NewRunoffTime;
}

//===============================================================================
int DLLEXPORT api_get_nodef_attribute(
    int node_idx, int attr, double* value)
//===============================================================================
{
    int error, bpat, tseries_idx;

    //printf("==== %d \n",attr);

    error = check_api_is_initialized("api_get_nodef_attribute");
    if (error) return error;

    switch (attr) {

        case nodef_type :
            *value = Node[node_idx].type;
            break;

        case nodef_outfall_type  :
            switch (Node[node_idx].type) {
                case OUTFALL :
                    *value = Outfall[Node[node_idx].subIndex].type;
                    break;
                default :    
                    *value = API_NULL_VALUE_I;
                    sprintf(errmsg, "Extracting nodef_outfall_type for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                    api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                    return api_err_wrong_type;
            }
            break;   

        case nodef_invertElev  :
            *value = FTTOM(Node[node_idx].invertElev);
            break;

        case nodef_fullDepth  :
            *value = FTTOM(Node[node_idx].fullDepth);
            break;      

        case nodef_initDepth  :
            switch (Node[node_idx].type) {
                case OUTFALL :
                    error = api_get_headBC(node_idx, StartDateTime, value);
                    if (error) return error;
                    *value -= FTTOM(Node[node_idx].invertElev);
                    break;
                default :
                    *value = FTTOM(Node[node_idx].initDepth);
            }
            break;

        case nodef_StorageConstant  :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aConst;
                    break;
                default :
                    *value = -1;
            }
            break; 

        case nodef_StorageCoeff :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aCoeff;
                    break;
                default :
                    *value = -1;
            }
            break;

        case nodef_StorageExponent  :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aExpon;
                    break;
                default :
                    *value = -1;
            }
            break;   

        case nodef_StorageCurveID  :
            switch (Node[node_idx].type) {
                case STORAGE :
                    *value = Storage[Node[node_idx].subIndex].aCurve + 1;
                    break;
                default :
                    *value = -1;
            }
            break;

        case nodef_extInflow_tSeries :
            if (Node[node_idx].extInflow)
                *value = Node[node_idx].extInflow->tSeries;
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_extInflow_tSeries for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;   

        case nodef_extInflow_tSeries_x1  :
            tseries_idx = Node[node_idx].extInflow->tSeries;
            if (tseries_idx >= 0)
                *value = Tseries[tseries_idx].x1;
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_extInflow_tSeries_x1 for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_extInflow_tSeries_x2  :
            tseries_idx = Node[node_idx].extInflow->tSeries;
            if (tseries_idx >= 0)
                *value = Tseries[tseries_idx].x2;
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting tseries_idx for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;     

        case nodef_extInflow_basePat_idx  :
            if (Node[node_idx].extInflow)
            {
                *value = Node[node_idx].extInflow->basePat;
                //*value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->basePat); // BRH20211221 HACK THIS IS WRONG
                //printf("%g \n",Node[node_idx].extInflow->cFactor);
                //printf("%d \n",Node[node_idx].extInflow->basePat);
            }    
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_extInflow_basePat for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_extInflow_basePat_type  :
            bpat = Node[node_idx].extInflow->basePat;
            //printf(" bpat %d",bpat);
            if (bpat >= 0) // baseline pattern exists
                *value = Pattern[bpat].type;
            else
            {
                *value = bpat;  // brh changed to bpat (-1) because API_NULL_VALUE_I does not have scope for where its needed
                //*value = API_NULL_VALUE_I;
                //printf("  bpat = %d \n", bpat);
                //printf("  location 3098705 problem with basePatType \n");
                // brh20211207s  commenting this so that it moves through with null result
                //sprintf(errmsg, "Extracting node_extInflow_basePat_type for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                //api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                //return api_err_wrong_type;
                // brh20211207e
            }
            break;   

        case nodef_extInflow_baseline :
            if (Node[node_idx].extInflow)
                *value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline);
            else
                *value = 0;
            break;

        case nodef_extInflow_sFactor  :
            if (Node[node_idx].extInflow)
                *value = Node[node_idx].extInflow->sFactor;
            else
                *value = 1;
            break;      

        case nodef_has_extInflow :
            if (Node[node_idx].extInflow)
                *value = 1;
            else
                *value = 0;
            break;   

        case nodef_dwfInflow_monthly_pattern  :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[0];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_monthly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_dwfInflow_daily_pattern :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[1];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_daily_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;    

        case nodef_dwfInflow_hourly_pattern  :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[2];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_hourly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;

        case nodef_dwfInflow_weekend_pattern  :
            if (Node[node_idx].dwfInflow)
                *value = Node[node_idx].dwfInflow->patterns[3];
            else
            {
                *value = API_NULL_VALUE_I;
                sprintf(errmsg, "Extracting node_dwfInflow_weekend_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
                api_report_writeErrorMsg(api_err_wrong_type, errmsg);
                return api_err_wrong_type;
            }
            break;   

        case nodef_dwfInflow_avgvalue  :
            if (Node[node_idx].dwfInflow)
                *value = CFTOCM(Node[node_idx].dwfInflow->avgValue);
            else
                *value = 0;
            break;

        case nodef_has_dwfInflow  :
            if (Node[node_idx].dwfInflow)
                *value = 1;
            else
                *value = 0;
            break;      

        case nodef_newDepth  :
            *value = FTTOM(Node[node_idx].newDepth);
            break;

        case nodef_inflow  :
            *value = CFTOCM(Node[node_idx].inflow);
            break;   

        case nodef_volume  :
            *value = CFTOCM(Node[node_idx].newVolume);
            break;

        case nodef_overflow  :
            *value = CFTOCM(Node[node_idx].overflow);
            break;     

        case nodef_rptFlag  :
            if (Node[node_idx].rptFlag)
                *value = 1;
            else
                *value = 0;
            break;

        default :
            printf(" ****** api_get_nodef_attribute called without supported attr at 3979874 %d ",attr);
            *value = API_NULL_VALUE_I;
    }

    // if (attr == nodef_type)
    // {
    //     // printf("    nodef_type attr = 2 \n");
    //     *value = Node[node_idx].type;
    // }
    // else if (attr == nodef_outfall_type)
    // {
    //     // printf("    nodef_outfall_type attr = 3 \n");
    //     if (Node[node_idx].type == OUTFALL)
    //     {
    //         *value = Outfall[Node[node_idx].subIndex].type;
    //     }
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting nodef_outfall_type for NODE %s, which is not an outfall [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_invertElev)
    // {
    //     // printf("    nodef_invertElev attr = 4 \n");
    //     *value = FTTOM(Node[node_idx].invertElev);
    // }
    // else if (attr == nodef_fullDepth)
    // {
    //     // printf("   nodef_fullDepth attr = 25 \n ");
    //     *value = FTTOM(Node[node_idx].fullDepth);
    // }
    // else if (attr == nodef_initDepth)
    // {
    //     // printf("    nodef_initDepth attr attr = 5 \n");
    //     if (Node[node_idx].type == OUTFALL)
    //     {
    //         // printf("   *** call api_get_headBC \n");
    //         error = api_get_headBC(node_idx, StartDateTime, value);
    //         if (error) return error;
    //         *value -= FTTOM(Node[node_idx].invertElev);
    //     }
    //     else
    //         *value = FTTOM(Node[node_idx].initDepth);
    // }
    // else if (attr == nodef_StorageConstant)
    // {
    //     // printf("   nodef_StorageConstant attr = 6 \n");
    //     if (Node[node_idx].type == STORAGE)
    //         *value = Storage[Node[node_idx].subIndex].aConst;
    //     else
    //         *value = -1;
    // }
    // else if (attr == nodef_StorageCoeff)
    // {
    //     // printf("   nodef_storage_Coeff attr = 7 \n");
    //     if (Node[node_idx].type == STORAGE)
    //         *value = Storage[Node[node_idx].subIndex].aCoeff;
    //     else
    //         *value = -1;
    // }
    // else if (attr == nodef_StorageExponent)
    // {
    //     // printf("   nodef_StorageExponent attr = 8 \n");
    //     if (Node[node_idx].type == STORAGE)
    //         *value = Storage[Node[node_idx].subIndex].aExpon;
    //     else
    //         *value = -1;
    // }
    // else if (attr == nodef_StorageCurveID)
    // {
    //     // printf("   nodef_StorageCurveID attr = 9 \n");
    //     if (Node[node_idx].type == STORAGE)
    //         *value = Storage[Node[node_idx].subIndex].aCurve + 1;
    //     else
    //         *value = -1;
    // }
    // else if (attr == nodef_extInflow_tSeries)
    // {
    //     //printf("   nodef_extInflow_tSeries attr = 10 \n");
    //     //printf("    %d \n ",Node[node_idx].extInflow->tSeries);
    //     if (Node[node_idx].extInflow)
    //         *value = Node[node_idx].extInflow->tSeries;
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_extInflow_tSeries for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_extInflow_tSeries_x1)
    // {
    //     // printf("   nodef_extInflow_tSeries_x1 attr = 11 \n");
    //     tseries_idx = Node[node_idx].extInflow->tSeries;
    //     if (tseries_idx >= 0)
    //         *value = Tseries[tseries_idx].x1;
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_extInflow_tSeries_x1 for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_extInflow_tSeries_x2)
    // {
    //     // printf("   nodef_extInflow_tSeries_x2 attr = 12 \n");
    //     tseries_idx = Node[node_idx].extInflow->tSeries;
    //     if (tseries_idx >= 0)
    //         *value = Tseries[tseries_idx].x2;
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting tseries_idx for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_extInflow_basePat)
    // {
    //     //printf("   nodef_extInflow_basePat attr = 13 \n");
    //     if (Node[node_idx].extInflow)
    //     {
    //         *value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->basePat);
    //         printf("%g \n",Node[node_idx].extInflow->cFactor);
    //         printf("%d \n",Node[node_idx].extInflow->basePat);
    //     }    
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_extInflow_basePat for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_extInflow_basePat_type)
    // {
    //     //printf("   nodef_extInflow_basePat_type attr = 14 \n");
    //     bpat = Node[node_idx].extInflow->basePat;
        
    //     //printf(" bpat %d",bpat);
    //     if (bpat >= 0) // baseline pattern exists
    //         *value = Pattern[bpat].type;
    //     else
    //     {
    //         *value = bpat;  // brh changed to bpat (-1) because API_NULL_VALUE_I does not have scope for where its needed
    //         //*value = API_NULL_VALUE_I;
    //         //printf("  bpat = %d \n", bpat);
    //         //printf("  location 3098705 problem with basePatType \n");
    //         // brh20211207s  commenting this so that it moves through with null result
    //         //sprintf(errmsg, "Extracting node_extInflow_basePat_type for NODE %s, which doesn't have an extInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         //api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         //return api_err_wrong_type;
    //         // brh20211207e
    //     }
    // }
    // else if (attr == nodef_extInflow_baseline)
    // {
    //     // printf("   nodef_extInflow_baseline attr = 15 \n");
    //     if (Node[node_idx].extInflow)
    //     {
    //         //printf(" in api baseline %f \n",Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline);
    //         *value = CFTOCM(Node[node_idx].extInflow->cFactor * Node[node_idx].extInflow->baseline);
    //     }
    //     else
    //         *value = 0;
    // }
    // else if (attr == nodef_extInflow_sFactor)
    // {
    //     // printf("   nodef_extInflow_sFactor attr = 16 \n");
    //     if (Node[node_idx].extInflow)
    //         *value = Node[node_idx].extInflow->sFactor;
    //     else
    //         *value = 1;
    // }
    // else if (attr == nodef_has_extInflow)
    // {
    //     // printf("   nodef_has_extInflow attr = 17 \n");
    //     if (Node[node_idx].extInflow)
    //     {
    //         *value = 1;
    //     }
    //     else
    //         *value = 0;
    // }
    // else if (attr == nodef_dwfInflow_monthly_pattern)
    // {
    //     // printf("   nodef_dwInflow_monthly_pattern attr = 18 \n");
    //     if (Node[node_idx].dwfInflow)
    //         *value = Node[node_idx].dwfInflow->patterns[0];
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_dwfInflow_monthly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_dwfInflow_daily_pattern)
    // {
    //     // printf("   nodef_dwdInflow_daily_pattern attr = 19 \n");
    //     if (Node[node_idx].dwfInflow)
    //     {
    //         *value = Node[node_idx].dwfInflow->patterns[1];
    //     }
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_dwfInflow_daily_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_dwfInflow_hourly_pattern)
    // {
    //     // printf("   nodef_dwInflow_hourly_pattern attr = 20 \n");
    //     if (Node[node_idx].dwfInflow)
    //         *value = Node[node_idx].dwfInflow->patterns[2];
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_dwfInflow_hourly_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_dwfInflow_weekend_pattern)
    // {
    //     // printf("   nodef_dwInflow_weekend_pattern attr = 21 \n");
    //     if (Node[node_idx].dwfInflow)
    //         *value = Node[node_idx].dwfInflow->patterns[3];
    //     else
    //     {
    //         *value = API_NULL_VALUE_I;
    //         sprintf(errmsg, "Extracting node_dwfInflow_weekend_pattern for NODE %s, which doesn't have a dwfInflow [api.c -> api_get_nodef_attribute]", Node[node_idx].ID);
    //         api_report_writeErrorMsg(api_err_wrong_type, errmsg);
    //         return api_err_wrong_type;
    //     }
    // }
    // else if (attr == nodef_dwfInflow_avgvalue)
    // {
    //     // printf("   nodef_dwInflow_avgvalue attr = 22 \n");
    //     if (Node[node_idx].dwfInflow)
    //         *value = CFTOCM(Node[node_idx].dwfInflow->avgValue);
    //     else
    //         *value = 0;
    // }
    // else if (attr == nodef_has_dwfInflow)
    // {
    //     // printf("   nodef_has_dwfInflow attr = 23 \n");
    //     if (Node[node_idx].dwfInflow)
    //         *value = 1;
    //     else
    //         *value = 0;
    // }
    // // brh20211207s
    // //else if (attr == node_depth)
    // else if (attr == nodef_newDepth)

    // {
    //     // printf("   node_depth = 24 \n");
    //     // printf("   nodef_newDepth attr = 24 \n");
    // // brh20211207e
    //     *value = FTTOM(Node[node_idx].newDepth);
    // }         
    // else if (attr == nodef_inflow)
    // {
    //     // printf("   nodef_inflow attr = 26 \n");
    //     *value = CFTOCM(Node[node_idx].inflow);
    // }
    // else if (attr == nodef_volume)
    // {
    //     // printf("   nodef_volume attr = 27 \n");
    //     *value = CFTOCM(Node[node_idx].newVolume);
    // }    
    // else if (attr == nodef_overflow)
    // {
    //     // printf("   nodef_overflow attr = 28 \n");
    //     *value = CFTOCM(Node[node_idx].overflow);
    // }
    // // brh20211207s
    // else if (attr = nodef_rptFlag)
    // {
    //     // printf("    nodef_rptFlag attr = 29 \n");
    //     if (Node[node_idx].rptFlag)
    //         *value = 1;
    //     else
    //         *value = 0;
    // }
    // // brh20211207e
    // else
    // {
    //     printf(" ****** api_get_nodef_attribute called without supported attr at 3979874 %d ",attr);
    //     *value = API_NULL_VALUE_I;
    // } 
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_linkf_attribute(
    int link_idx, int attr, double* value)
//===============================================================================
{
    int error;

    error = check_api_is_initialized("api_get_linkf_attribute");
    if (error) return error;

    switch (attr) {

        case linkf_subIndex :
            *value = Link[link_idx].subIndex;
            break;

        case linkf_type : 
            *value = Link[link_idx].type;
            break;

        case linkf_node1 :
            *value = Link[link_idx].node1;
            break;

        case linkf_node2 : 
            *value = Link[link_idx].node2;
            break;

        case linkf_offset1 :
            *value = FTTOM(Link[link_idx].offset1);
            break;

        case linkf_offset2 : 
            *value = FTTOM(Link[link_idx].offset2);
            break;

        case linkf_xsect_type :
            *value = Link[link_idx].xsect.type;
            break;

        case linkf_xsect_wMax :
            *value = FTTOM(Link[link_idx].xsect.wMax); 
            break;

        case linkf_xsect_yBot :
            *value = FTTOM(Link[link_idx].xsect.yBot);
            break;

        case linkf_xsect_yFull : 
            *value = FTTOM(Link[link_idx].xsect.yFull);
            break;

        case linkf_q0 :
            *value = CFTOCM(Link[link_idx].q0);
            break;

        case linkf_pump_type : 
            *value =  Pump[Link[link_idx].subIndex].type;
            break;
        case linkf_orifice_type  :
            *value = Orifice[Link[link_idx].subIndex].type;
            break;

        case linkf_outlet_type : 
            *value = Outlet[Link[link_idx].subIndex].curveType;
            break;

        case linkf_weir_type : 
            *value = Weir[Link[link_idx].subIndex].type;
            break;

        case linkf_conduit_roughness :
            switch (Link[link_idx].type) {
                case CONDUIT :
                    *value = Conduit[Link[link_idx].subIndex].roughness;
                    break;
                default :
                    *value = 0;
            }    
            break;

        case linkf_conduit_length : 
            switch (Link[link_idx].type) {
                case CONDUIT :
                    *value = FTTOM(Conduit[Link[link_idx].subIndex].length);
                    break;
                case ORIFICE :
                    *value = 0.01;
                    break;
                case WEIR :
                    *value = 0.01;
                case OUTLET :
                    *value = 0.01;
                    break;
                case PUMP :
                    *value = 0.01;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_weir_end_contractions :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].endCon;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_curveid :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].cdCurve+1;
                    break;
                case PUMP :
                    *value = Pump[Link[link_idx].subIndex].pumpCurve+1;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].qCurve+1;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_discharge_coeff1 :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].cDisch1;
                    break;
                case ORIFICE :
                    *value = Orifice[Link[link_idx].subIndex].cDisch;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].qCoeff;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_discharge_coeff2 :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].cDisch2;
                    break;
                case OUTLET :
                    *value = Outlet[Link[link_idx].subIndex].qExpon;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_weir_side_slope :
            switch (Link[link_idx].type) {
                case WEIR :
                    *value = Weir[Link[link_idx].subIndex].slope;
                    break;
                default :
                    *value = 0;
            }
            break;

        case linkf_flow :
            *value = CFTOCM(Link[link_idx].newFlow);
            break;

        case linkf_depth :
            *value = FTTOM(Link[link_idx].newDepth);
            break;

        case linkf_volume :
            *value = CFTOCM(Link[link_idx].newVolume);
            break;

        case linkf_froude :
            *value = Link[link_idx].froude;
            break;

        case linkf_setting :
            *value = Link[link_idx].setting;
            break;

        case linkf_left_slope :
            *value = api->double_vars[api_left_slope][link_idx];
            break;

        case linkf_right_slope :
            *value = api->double_vars[api_right_slope][link_idx];
            break;

        case linkf_geometry :
            printf(" ****** api_get_linkf_attribute called for unsupported attr = linkf_geometry at 2875 %d ",attr);
            break;
            case linkf_rptFlag :
                if (Link[link_idx].rptFlag)
                    *value = 1;
                else
                    *value = 0; 
            break;        
        default :             
            printf(" ****** api_get_linke_attribute called without supported attr at 837954 %d ",attr);
            *value = API_NULL_VALUE_I;               
    }

    // if (attr == linkf_subIndex)
    //     *value = Link[link_idx].subIndex;
    // else if (attr == linkf_type)
    //     *value = Link[link_idx].type;
    // else if (attr == linkf_node1)
    //     *value = Link[link_idx].node1;
    // else if (attr == linkf_node2)
    //     *value = Link[link_idx].node2;
    // else if (attr == linkf_offset1)
    //     *value = FTTOM(Link[link_idx].offset1);
    // else if (attr == linkf_offset2)
    //     *value = FTTOM(Link[link_idx].offset2);
    // else if (attr == linkf_xsect_type)
    //     *value = Link[link_idx].xsect.type;
    // else if (attr == linkf_xsect_wMax)
    //     *value = FTTOM(Link[link_idx].xsect.wMax);
    // else if (attr == linkf_xsect_yBot)
    //     *value = FTTOM(Link[link_idx].xsect.yBot);
    // else if (attr == linkf_xsect_yFull)
    //     *value = FTTOM(Link[link_idx].xsect.yFull);
    // else if (attr == linkf_q0)
    //     *value = CFTOCM(Link[link_idx].q0);
    // // brh20211207s duplicate ov above    
    // // rm else if (attr == linkf_type)
    // // rm    *value =  Link[link_idx].type;
    // // brh20211207e
    // else if (attr == linkf_pump_type)
    //     *value =  Pump[Link[link_idx].subIndex].type;
    // else if (attr == linkf_orifice_type)
    //     *value = Orifice[Link[link_idx].subIndex].type;
    // else if (attr == linkf_outlet_type)
    //     *value = Outlet[Link[link_idx].subIndex].curveType;
    // else if (attr == linkf_weir_type)
    //     *value = Weir[Link[link_idx].subIndex].type;
    // else if (attr == linkf_conduit_roughness)
    // {
    //     if (Link[link_idx].type == CONDUIT)
    //         *value = Conduit[Link[link_idx].subIndex].roughness;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_conduit_length)
    // {
    //     if (Link[link_idx].type == CONDUIT)
    //         *value = FTTOM(Conduit[Link[link_idx].subIndex].length);

    //     else if (Link[link_idx].type == ORIFICE)
    //         *value = 0.01;
    //     else if (Link[link_idx].type == WEIR)
    //         *value = 0.01;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = 0.01;
    //     else if (Link[link_idx].type == PUMP)
    //         *value = 0.01;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_weir_end_contractions)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].endCon;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_curveid)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].cdCurve+1;
    //     else if (Link[link_idx].type == PUMP)
    //         *value = Pump[Link[link_idx].subIndex].pumpCurve+1;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = Outlet[Link[link_idx].subIndex].qCurve+1;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_discharge_coeff1)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].cDisch1;
    //     else if (Link[link_idx].type == ORIFICE)
    //         *value = Orifice[Link[link_idx].subIndex].cDisch;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = Outlet[Link[link_idx].subIndex].qCoeff;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_discharge_coeff2)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].cDisch2;
    //     else if (Link[link_idx].type == OUTLET)
    //         *value = Outlet[Link[link_idx].subIndex].qExpon;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_weir_side_slope)
    // {
    //     if (Link[link_idx].type == WEIR)
    //         *value = Weir[Link[link_idx].subIndex].slope;
    //     else
    //         *value = 0;
    // }
    // else if (attr == linkf_flow)
    //     *value = CFTOCM(Link[link_idx].newFlow);
    // else if (attr == linkf_depth)
    //     *value = FTTOM(Link[link_idx].newDepth);
    // else if (attr == linkf_volume)
    //     *value = CFTOCM(Link[link_idx].newVolume);
    // else if (attr == linkf_froude)
    //     *value = Link[link_idx].froude;
    // else if (attr == linkf_setting)
    //     *value = Link[link_idx].setting;
    // else if (attr == linkf_left_slope)
    //     *value = api->double_vars[api_left_slope][link_idx];
    // else if (attr == linkf_right_slope)
    //     *value = api->double_vars[api_right_slope][link_idx];
    // // brh 20211207s
    // else if (attr == linkf_geometry)  
    // {
    //     printf(" ****** api_get_linkf_attribute called for unsupported attr = linkf_geometry at 2875 %d ",attr); 
    // }    
    // else if (attr == linkf_rptFlag)
    // {
    //     if (Link[link_idx].rptFlag)
    //         *value = 1;
    //     else
    //         *value = 0;   
    // }        
    // // brh 20211207e    
    // else
    // {
    //     printf(" ****** api_get_linke_attribute called without supported attr at 837954 %d ",attr);
    //     *value = API_NULL_VALUE_I;
    // }    
    return 0;
}
//===============================================================================
int DLLEXPORT api_get_num_objects(
    int object_type)
//===============================================================================
{
    int error;
    error = check_api_is_initialized("api_get_num_objects");
    if (error) return error;
    // if (object_type > API_START_INDEX) // Objects for API purposes
    //     return api->num_objects[object_type - API_START_INDEX];
    return Nobjects[object_type];
}
//===============================================================================
int DLLEXPORT api_get_object_name(
    int object_idx, char* object_name, int object_type)
//===============================================================================
{
    int error, ii;
    int obj_len = -1;

    error = check_api_is_initialized("api_get_object_name");
    if (error) return error;
    error = api_get_object_name_len(object_idx, object_type, &obj_len);
    if (error) return error;

    // switch (object_type) {
    //     case NODE :
    //         for(ii=0; ii<obj_len; ii++)
    //         {
    //             object_name[ii] = Node[object_idx].ID[ii];
    //         }
    //         break;
    //     case LINK :
    //         for(ii=0; ii<obj_len; ii++)
    //         {
    //             object_name[ii] = Link[object_idx].ID[ii];
    //         }
    //     default :
    //         strcpy(object_name, "");
    //         sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name]", object_type);
    //         api_report_writeErrorMsg(api_err_not_developed, errmsg);
    //         return api_err_not_developed;
    // }

    if (object_type == NODE)
    {
        for(ii=0; ii<obj_len; ii++)
        {
            object_name[ii] = Node[object_idx].ID[ii];
        }
    }
    else if (object_type == LINK)
    {
        for(ii=0; ii<obj_len; ii++)
        {
            object_name[ii] = Link[object_idx].ID[ii];
        }
    }
    else
    {
        strcpy(object_name, "");
        sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name]", object_type);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}
//===============================================================================
int DLLEXPORT api_get_object_name_len(
    int object_idx, int object_type, int* len)
//===============================================================================
{
    int error;

    error = check_api_is_initialized("api_get_object_name_len");
    if (error) {
        *len = API_NULL_VALUE_I;
        return error;
    }

    switch (object_type) {
        case NODE :
            *len = strlen(Node[object_idx].ID);
            return 0;
            break;
        case LINK :
            *len = strlen(Link[object_idx].ID);
            return 0;
            break;
        default :
            *len = API_NULL_VALUE_I;
            sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name_len]", object_type);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            return api_err_not_developed;
    }

    // if (object_type == NODE)
    // {
    //     *len = strlen(Node[object_idx].ID);
    //     return 0;
    // }
    // else if (object_type == LINK)
    // {
    //     *len = strlen(Link[object_idx].ID);
    //     return 0;
    // }
    // else
    // {
    //     *len = API_NULL_VALUE_I;
    //     sprintf(errmsg, "OBJECT_TYPE %d [api.c -> api_get_object_name_len]", object_type);
    //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
    //     return api_err_not_developed;
    // }
}
//===============================================================================
int DLLEXPORT api_get_num_table_entries(
    int table_idx, int table_type, int* num_entries)
//===============================================================================
{
    double xx, yy;
    int success;

    *num_entries = 0;
    // printf("1 Number of entries in Curve %d\n", *num_entries);

    switch (table_type) {
        case CURVE :
            // ERROR handling
            if (table_idx >= Nobjects[CURVE] || Nobjects[CURVE] == 0) return -1;
            success = table_getFirstEntry(&Curve[table_idx], &xx, &yy); // first values in the table
            (*num_entries)++;
            while (success)
            {
                success = table_getNextEntry(&(Curve[table_idx]), &xx, &yy);
                if (success) (*num_entries)++;
                // printf("0 Number of entries in Curve %d\n", *num_entries);
            }    
            break;
        default :
            return -1;
    }

    // if (table_type == CURVE)
    // {
    //     // ERROR handling
    //     if (table_idx >= Nobjects[CURVE] || Nobjects[CURVE] == 0) return -1;
    //     success = table_getFirstEntry(&Curve[table_idx], &x, &y); // first values in the table
    //     (*num_entries)++;
    //     while (success)
    //     {
    //         success = table_getNextEntry(&(Curve[table_idx]), &x, &y);
    //         if (success) (*num_entries)++;
    //         // printf("0 Number of entries in Curve %d\n", *num_entries);
    //     }
    // }
    // else
    // {
    //     return -1;
    // }
    // // printf("Number of entries in Curve %d\n", *num_entries);
    // // printf("SUCCESS %d\n", success);

    return 0;
}
//===============================================================================
int DLLEXPORT api_get_table_attribute(
    int table_idx, int attr, double* value)
//===============================================================================
{
    int error;

    error = check_api_is_initialized("api_get_table_attribute");
    if (error) return error;

    switch (attr) {
        case table_ID :
            *value = table_idx;
            break;
        case table_type :
            *value = Curve[table_idx].curveType;
            break;
        case table_refers_to :
            *value = Curve[table_idx].refersTo;
            break;
        default :
            *value = API_NULL_VALUE_I;
            sprintf(errmsg, "attr %d [api.c -> api_get_table_attribute]", attr);
            api_report_writeErrorMsg(api_err_not_developed, errmsg);
            return api_err_not_developed;
    }

    // if (attr == table_ID)
    //     *value = table_idx;
    // else if (attr == table_type)
    //     *value = Curve[table_idx].curveType;
    // else if (attr == table_refers_to)
    //     *value = Curve[table_idx].refersTo;
    // else
    // {
    //     *value = API_NULL_VALUE_I;
    //     sprintf(errmsg, "attr %d [api.c -> api_get_table_attribute]", attr);
    //     api_report_writeErrorMsg(api_err_not_developed, errmsg);
    //     return api_err_not_developed;
    // }
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_first_entry_table(
    int table_idx, int table_type, double* xx, double* yy)
//===============================================================================
{
    int success;

    switch (table_type) {
        case CURVE :
            success = table_getFirstEntry(&(Curve[table_idx]), xx, yy);
            // printf("...success, %d \n",success);
            // printf("...curveType, %d \n",Curve[table_idx].curveType);
            // unit conversion depending on the type of curve
            switch (Curve[table_idx].curveType) {
                case STORAGE_CURVE:
                    *xx = FTTOM(*xx);
                    *yy = FT2TOM2(*yy);
                    break;
                case DIVERSION_CURVE:
                    *xx = CFTOCM(*xx);
                    *yy = CFTOCM(*yy);
                    break;
                case TIDAL_CURVE:
                    *yy = FTTOM(*yy);
                    break;
                case RATING_CURVE:
                    *xx = FTTOM(*xx);
                    *yy = CFTOCM(*yy);
                    break;
                case SHAPE_CURVE:
                    *xx = FTTOM(*xx);
                    *yy = FTTOM(*yy);
                    break;
                case CONTROL_CURVE:
                    break;
                case WEIR_CURVE:
                    *xx = FTTOM(*xx);
                    break;
                case PUMP1_CURVE:
                    *xx = CFTOCM(*xx);
                    *yy = CFTOCM(*yy);
                    break;
                case PUMP2_CURVE:
                    *xx = FTTOM(*xx);
                    *yy = CFTOCM(*yy);
                    break;
                case PUMP3_CURVE:
                    *xx = FTTOM(*xx);
                    *yy = CFTOCM(*yy);
                    break;
                case PUMP4_CURVE:
                    *xx = FTTOM(*xx);
                    *yy = CFTOCM(*yy);
                    break;
                default:
                    return 0;
            }
            break;
        case TSERIES :
            success = table_getFirstEntry(&(Tseries[table_idx]), xx, yy);
            break;          
        default :
            return 0;
    }
    // if (table_type == CURVE)
    //     success = table_getFirstEntry(&(Curve[table_idx]), x, y);
    // else if (table_type == TSERIES)
    //     success = table_getFirstEntry(&(Tseries[table_idx]), x, y);
    // else
    //     return 0;

    return success;
}

//===============================================================================
int DLLEXPORT api_get_next_entry_table(
    int table_idx, int table_type, double* xx, double* yy)
//===============================================================================
{
    int success;

    switch (table_type) {
        case TSERIES :
            success = table_getNextEntry(&(Tseries[table_idx]), &(Tseries[table_idx].x2), &(Tseries[table_idx].y2));
            if (success)
            {
                *xx = Tseries[table_idx].x2;
                *yy = Tseries[table_idx].y2;
            }
            break;
        case CURVE :
            success = table_getNextEntry(&(Curve[table_idx]), &(Curve[table_idx].x2), &(Curve[table_idx].y2));
            if (success)
            {
                    // unit conversion depending on the type of curve
            switch (Curve[table_idx].curveType) {
                case STORAGE_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = FT2TOM2(Curve[table_idx].y2);
                    // printf("...xx, %f \n",*xx);
                    // printf("...yy, %f \n",*yy);
                    break;
                case DIVERSION_CURVE:
                    *xx = CFTOCM(Curve[table_idx].x2);
                    *yy = CFTOCM(Curve[table_idx].y2);
                    break;
                case TIDAL_CURVE:
                    *xx = Curve[table_idx].x2;
                    *yy = FTTOM(Curve[table_idx].y2);
                    break;
                case RATING_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = CFTOCM(Curve[table_idx].y2);
                    break;
                case CONTROL_CURVE:
                    *xx = Curve[table_idx].x2;
                    *yy = Curve[table_idx].y2;
                    break;
                case SHAPE_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = FTTOM(Curve[table_idx].y2);
                    break;
                case WEIR_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = Curve[table_idx].y2;
                    break;
                case PUMP1_CURVE:
                    *xx = CFTOCM(Curve[table_idx].x2);
                    *yy = CFTOCM(Curve[table_idx].y2);
                    break;
                case PUMP2_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = CFTOCM(Curve[table_idx].y2);
                    break;break;
                case PUMP3_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = CFTOCM(Curve[table_idx].y2);
                    break;
                case PUMP4_CURVE:
                    *xx = FTTOM(Curve[table_idx].x2);
                    *yy = CFTOCM(Curve[table_idx].y2);
                    break;
                default:
                    *xx = Curve[table_idx].x2;
                    *yy = Curve[table_idx].y2;
                    break;
                }
            }
            break;
        default :
            return 0;
    }
    // if (table_type == TSERIES)
    // {
    //     success = table_getNextEntry(&(Tseries[table_idx]), &(Tseries[table_idx].x2), &(Tseries[table_idx].y2));
    //     if (success)
    //     {
    //         *x = Tseries[table_idx].x2;
    //         *y = Tseries[table_idx].y2;
    //     }
    // }
    // else if (table_type == CURVE)
    // {
    //     success = table_getNextEntry(&(Curve[table_idx]), &(Curve[table_idx].x2), &(Curve[table_idx].y2));
    //     if (success)
    //     {
    //         *x = Curve[table_idx].x2;
    //         *y = Curve[table_idx].y2;
    //     }
    // }

    return success;
}

//===============================================================================
int DLLEXPORT api_get_next_entry_tseries(
    int tseries_idx)
//===============================================================================
{
    int success;
    double x2, y2;

    x2 = Tseries[tseries_idx].x2;
    y2 = Tseries[tseries_idx].y2;
    success = table_getNextEntry(&(Tseries[tseries_idx]), &(Tseries[tseries_idx].x2), &(Tseries[tseries_idx].y2));
    if (success == TRUE)
    {
        Tseries[tseries_idx].x1 = x2;
        Tseries[tseries_idx].y1 = y2;
    }
    return success;
}

//===============================================================================
// --- Output Writing (Post Processing)
// * The follwing functions should only be executed after finishing
//   and writing SWMM5+ report files. The following functions are
//   meant to be called from Fortran in order to export .rpt and
//   .out files according to the SWMM 5.13 standard. Fortran-generated
//   report files are not manipulated here, the manipulation of
//   SWMM5+ report files is kept within the Fortran code to ensure
//   compatibility with future updates of the SWMM5+ standard
//===============================================================================

//===============================================================================
int DLLEXPORT api_write_output_line(
    double t)
//===============================================================================
// t: elapsed time in seconds
{

    // --- check that simulation can proceed
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // Update routing times to skip interpolation when saving results.
    OldRoutingTime = 0; NewRoutingTime = t*1000; // times in msec
    output_saveResults(t*1000);
    return 0;
}

//===============================================================================
int DLLEXPORT api_update_nodeResult(
    int node_idx, int resultType, double newNodeResult)
//===============================================================================
{
    // --- check that simulation can proceed
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    if (resultType == output_node_depth)
        Node[node_idx].newDepth = newNodeResult;
    else if (resultType == output_node_volume)
        Node[node_idx].newVolume = newNodeResult;
    else if (resultType == output_node_latflow)
        Node[node_idx].newLatFlow = newNodeResult;
    else if (resultType == output_node_inflow)
        Node[node_idx].inflow = newNodeResult;
    else
    {
        sprintf(errmsg, "resultType %d [api.c -> api_update_nodeResult]", resultType);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}

//===============================================================================
int DLLEXPORT api_update_linkResult(
    int link_idx, int resultType, double newLinkResult)
//===============================================================================
{
    // --- check that simulation can proceed
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    if (resultType == output_link_depth)
        Link[link_idx].newDepth = newLinkResult;
    else if (resultType == output_link_flow)
        Link[link_idx].newFlow = newLinkResult;
    else if (resultType == output_link_volume)
        Link[link_idx].newVolume = newLinkResult;
    else if (resultType == output_link_direction)
        Link[link_idx].direction = newLinkResult;
    else
    {
        sprintf(errmsg, "resultType %d [api.c -> api_update_linkResult]", resultType);
        api_report_writeErrorMsg(api_err_not_developed, errmsg);
        return api_err_not_developed;
    }
    return 0;
}

//===============================================================================
// --- Print-out
//===============================================================================

//===============================================================================
int DLLEXPORT api_export_linknode_properties(
    int units)
//===============================================================================
{
    //  link
    int li_idx[Nobjects[LINK]];
    int li_link_type[Nobjects[LINK]];
    int li_geometry[Nobjects[LINK]];
    int li_Mnode_u[Nobjects[LINK]];
    int li_Mnode_d[Nobjects[LINK]];
    float lr_Length[Nobjects[LINK]];
    float lr_Slope[Nobjects[LINK]];
    float lr_Roughness[Nobjects[LINK]];
    float lr_InitialFlowrate[Nobjects[LINK]];
    float lr_InitialUpstreamDepth[Nobjects[LINK]];
    float lr_InitialDnstreamDepth[Nobjects[LINK]];
    int li_InitialDepthType[Nobjects[LINK]]; //
    float lr_BreadthScale[Nobjects[LINK]]; //
    float lr_InitialDepth[Nobjects[LINK]]; //

    //  node
    int ni_idx[Nobjects[NODE]];
    int ni_node_type[Nobjects[NODE]];
    int ni_N_link_u[Nobjects[NODE]];
    int ni_N_link_d[Nobjects[NODE]];
    int ni_Mlink_u1[Nobjects[NODE]];
    int ni_Mlink_u2[Nobjects[NODE]];
    int ni_Mlink_u3[Nobjects[NODE]];
    int ni_Mlink_d1[Nobjects[NODE]];
    int ni_Mlink_d2[Nobjects[NODE]];
    int ni_Mlink_d3[Nobjects[NODE]];
    float nr_Zbottom[Nobjects[NODE]]; //

    int NNodes = Nobjects[NODE];
    int NLinks = Nobjects[LINK];

    float length_units;
    float manning_units;
    float flow_units;

    FILE *f_nodes;
    FILE *f_links;

    int i;
    int error;

    error = check_api_is_initialized("api_export_linknode_properties");
    if (error) return error;

    // Initialization
    for (i = 0; i<Nobjects[NODE]; i++) {
        ni_N_link_u[i] = 0;
        ni_N_link_d[i] = 0;
        ni_Mlink_u1[i] = API_NULL_VALUE_I;
        ni_Mlink_u2[i] = API_NULL_VALUE_I;
        ni_Mlink_u3[i] = API_NULL_VALUE_I;
        ni_Mlink_d1[i] = API_NULL_VALUE_I;
        ni_Mlink_d2[i] = API_NULL_VALUE_I;
        ni_Mlink_d3[i] = API_NULL_VALUE_I;
    }

    // Choosing unit system
    if (units == US)
    {
        flow_units = 1;
        manning_units = 1;
        length_units = 1;
    }
    else if (units == SI)
    {
        flow_units = M3perFT3;
        manning_units = pow(1/MperFT, 1/3);
        length_units = MperFT;
    }
    else
    {
        sprintf(errmsg, "Incorrect type of units %d [api.c -> api_export_linknode_properties]", units);
        api_report_writeErrorMsg(api_err_wrong_type, errmsg);
        return api_err_wrong_type;
    }

    // Links
    for (i=0; i<Nobjects[LINK]; i++) {
        int li_sub_idx;
        float h;

        li_idx[i] = i;
        li_link_type[i] = Link[i].type;
        li_geometry[i] = Link[i].xsect.type;

        li_Mnode_u[i] = Link[i].node1;
        error = add_link(i, li_Mnode_u[i], DOWNSTREAM, 
            ni_N_link_u, 
                ni_Mlink_u1, 
                ni_Mlink_u2, 
                ni_Mlink_u3, 
            ni_N_link_d, 
                ni_Mlink_d1, 
                ni_Mlink_d2, 
                ni_Mlink_d3);
        if (error) return error;

        li_Mnode_d[i] = Link[i].node2;
        error = add_link(i, li_Mnode_d[i], UPSTREAM, 
            ni_N_link_u, 
                ni_Mlink_u1, 
                ni_Mlink_u2, 
                ni_Mlink_u3, 
            ni_N_link_d, 
                ni_Mlink_d1, 
                ni_Mlink_d2, 
                ni_Mlink_d3);
        if (error) return error;

        li_sub_idx = Link[i].subIndex;
        // [li_InitialDepthType] This condition is associated to nodes in SWMM
        if (li_link_type[i] == CONDUIT) {
            lr_Length[i] = Conduit[li_sub_idx].length * length_units;
            lr_Roughness[i] = Conduit[li_sub_idx].roughness * manning_units;
            h = (Node[li_Mnode_u[i]].invertElev - Node[li_Mnode_d[i]].invertElev) * length_units;
            lr_Slope[i] = -SSIGN(h)*lr_Length[i]/(pow(lr_Length[i],2) - pow(fabsf(h),2));
        } else {
            lr_Length[i] = 0;
            lr_Roughness[i] = 0;
            lr_Slope[i] = 0;
        }

        lr_InitialFlowrate[i] = Link[i].q0 * flow_units;
        lr_InitialUpstreamDepth[i] = Node[li_Mnode_u[i]].initDepth * length_units;
        lr_InitialDnstreamDepth[i] = Node[li_Mnode_d[i]].initDepth * length_units;
    }

    // Nodes
    for (i=0; i<Nobjects[NODE]; i++) {
        ni_idx[i] = i;
        ni_node_type[i] = Node[i].type;
    }

    f_nodes = fopen("debug_input/node/nodes_info.csv", "w");
    f_links = fopen("debug_input/link/links_info.csv", "w");

    fprintf(f_nodes,
        "n_left,node_id,ni_idx,ni_node_type,ni_N_link_u,ni_N_link_d,ni_Mlink_u1,ni_Mlink_u2,ni_Mlink_u3,ni_Mlink_d1,ni_Mlink_d2,ni_Mlink_d3\n");
    for (i=0; i<NNodes; i++) {
        fprintf(f_nodes, "%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
            NNodes-i,
            Node[i].ID,
            ni_idx[i],
            ni_node_type[i],
            ni_N_link_u[i],
            ni_N_link_d[i],
            ni_Mlink_u1[i],
            ni_Mlink_u2[i],
            ni_Mlink_u3[i],
            ni_Mlink_d1[i],
            ni_Mlink_d2[i],
            ni_Mlink_d3[i]);
    }
    fclose(f_nodes);

    fprintf(f_links,
        "l_left,link_id,li_idx,li_link_type,li_geometry,li_Mnode_u,li_Mnode_d,lr_Length,lr_Slope,lr_Roughness,lr_InitialFlowrate,lr_InitialUpstreamDepth,lr_InitialDnstreamDepth\n");
    for (i=0; i<NLinks; i++) {
        fprintf(f_links, "%d,%s,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
            NLinks-i,
            Link[i].ID,
            li_idx[i],
            li_link_type[i],
            li_geometry[i],
            li_Mnode_u[i],
            li_Mnode_d[i],
            lr_Length[i],
            lr_Slope[i],
            lr_Roughness[i],
            lr_InitialFlowrate[i],
            lr_InitialUpstreamDepth[i],
            lr_InitialDnstreamDepth[i]);
    }
    fclose(f_links);

    return 0;
}

//===============================================================================
int DLLEXPORT api_export_link_results(
    int link_idx)
//===============================================================================
{
	FILE* tmp;
    DateTime days;
    int period;
    char theTime[20];
    char theDate[20];
	char path[50];
    int error;

    error = check_api_is_initialized("api_export_link_results");
    if (error) return error;

    /* File path writing */
    strcpy(path, "debug_output/swmm5/link/");
    strcat(path, Link[link_idx].ID); strcat(path, ".csv");
    tmp = fopen(path, "w");
    fprintf(tmp, "date,time,flow,velocity,depth,volume,capacity\n");

    for ( period = 1; period <= Nperiods; period++ )
    {
        output_readDateTime(period, &days);
        datetime_dateToStr(days, theDate);
        datetime_timeToStr(days, theTime);
        output_readLinkResults(period, link_idx);
        fprintf(tmp, "%10s,%8s,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            theDate,
            theTime,
            LinkResults[LINK_FLOW],
            LinkResults[LINK_VELOCITY],
            LinkResults[LINK_DEPTH],
            LinkResults[LINK_VOLUME],
            LinkResults[LINK_CAPACITY]);
    }
    fclose(tmp);

    return 0;
}

//===============================================================================
int DLLEXPORT api_export_node_results(
    int node_idx)
//===============================================================================
{
	FILE* tmp;
    DateTime days;
    int period;
    char theTime[20];
    char theDate[20];
	char path[50];
    int error;

    error = check_api_is_initialized("api_export_node_results");
    if (error) return error;

    if (stat("NodeResults", &st) == -1) {
        mkdir("NodeResults", 0700);
    }

    /* File path writing */
    strcpy(path, "NodeResults/");
    strcat(path, Node[node_idx].ID);
    strcat(path, ".csv");
    tmp = fopen(path, "w");
    fprintf(tmp, "date,time,inflow,overflow,depth,volume\n");

    for ( period = 1; period <= Nperiods; period++ ) {
        output_readDateTime(period, &days);
        datetime_dateToStr(days, theDate);
        datetime_timeToStr(days, theTime);
        output_readNodeResults(period, node_idx);
        fprintf(tmp, "%10s,%8s,%.4f,%.4f,%.4f,%.4f\n",
            theDate,
            theTime,
            NodeResults[NODE_INFLOW],
            NodeResults[NODE_OVERFLOW],
            NodeResults[NODE_DEPTH],
            NodeResults[NODE_VOLUME]);
    }
    fclose(tmp);
    return 0;
}

//===============================================================================
// --- Utils
//===============================================================================

//===============================================================================
int DLLEXPORT api_find_object(
    int object_type, char *id)
//===============================================================================
{
    return project_findObject(object_type, id);
}


// -------------------------------------------------------------------------
// |
// |  Hydrology
// v
// -------------------------------------------------------------------------
//===============================================================================
int DLLEXPORT api_call_runoff_execute()
//===============================================================================
// calls the runoff_execute() procedure in SWMM-C
{
    // printf(" in api_call_runoff_execute");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    runoff_execute();
    
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_subcatch_runoff(
    int sc_idx, double *runoff)
//===============================================================================
{
    // printf(" in api_get_subcatch_runoff \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // Get runoff and convert to cubic meters per second
    *runoff = CFTOCM(Subcatch[sc_idx].newRunoff);
    // printf("... sc_idx, newRunoff CMS %d , %f \n",sc_idx,Subcatch[sc_idx].newRunoff);
    
    return 0;
}

//===============================================================================
int DLLEXPORT api_get_subcatch_runoff_nodeIdx(
    int sc_idx, int *node_idx)
//===============================================================================
{
    // printf(" in api_get_subcatch_runoff_nodeIdx \n");

    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( ! api->IsInitialized )
    {
        report_writeErrorMsg(ERR_NOT_OPEN, "");
        return error_getCode(ErrorCode);
    }

    // Get node index
    *node_idx = Subcatch[sc_idx].outNode;

    //printf("... sc_idx, node_idx %d , %d \n",sc_idx,*node_idx);
    
    return 0;
}

// //===============================================================================
// int DLLEXPORT api_get_subcatch_runoff_nodeName(
//     int sc_idx, int *nodeIdx)
// //===============================================================================
// {
//     printf(" in api_get_subcatch_runoff_nodeName");

//     if ( ErrorCode ) return error_getCode(ErrorCode);
//     if ( ! api->IsInitialized )
//     {
//         report_writeErrorMsg(ERR_NOT_OPEN, "");
//         return error_getCode(ErrorCode);
//     }

    
//     *runoff = CFTOCM(Subcatch[sc_idx].newRunoff);
//     printf("... sc_idx, newRunoff CMS %d , %f \n",sc_idx,Subcatch[sc_idx].newRunoff);

//     api_get_object_name(
//     int sc_indx, char* object_name, int object_type)
    
//     return 0;
// }

// -------------------------------------------------------------------------
// |
// |  Private functionalities
// v
// -------------------------------------------------------------------------
//===============================================================================
int api_load_vars()
//===============================================================================
{
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    int sect, ii, jj, kk, error;
    int found = 0;
    double xx[4];

    error = check_api_is_initialized("api_load_vars");
    if (error) return error;

    for (ii = 0; ii < NUM_API_DOUBLE_VARS; ii++)
    {
        api->double_vars[ii] = (double*) calloc(Nobjects[LINK], sizeof(double));
    }

    rewind(Finp.file);
    while ( fgets(line, MAXLINE, Finp.file) != NULL )
    {
        // --- make copy of line and scan for tokens
        strcpy(wLine, line);
        Ntokens = getTokens(wLine);

        // --- skip blank lines and comments
        if ( Ntokens == 0 ) continue;
        if ( *Tok[0] == ';' ) continue;

        if (*Tok[0] == '[')
        {
            if (found) break;
            sect = findmatch(Tok[0], SectWords);
        }
        else
        {
            if (sect == s_XSECTION)
            {
                found = 1;
                jj = project_findObject(LINK, Tok[0]);
                kk = findmatch(Tok[1], XsectTypeWords);
                if ( kk == TRAPEZOIDAL )
                {
                    // --- parse and save geometric parameters
                    for (ii = 2; ii <= 5; ii++)
                        getDouble(Tok[ii], &xx[ii-2]);

                    // --- extract left and right slopes for trapezoidal channel
                    api->double_vars[api_left_slope] [jj] = xx[2];
                    api->double_vars[api_right_slope][jj] = xx[3];
                }
            }
        }
        continue;
    }
    return 0;
}
//===============================================================================
// int add_link_alt(
//     int li_idx,
//     int ni_idx,
//     int maxUp,
//     int maxDn,
//     int direction,
//     int* ni_N_link_u,
//     int* ni_N_link_d,
//     int* MlinkUp[3],
//     int*,MlinkDn[3])
// //===============================================================================
// {
//     if (direction == UPSTREAM) {
//         ni_N_link_u[ni_idx] ++;
//         if (ni_N_link_up[ni_idx] <= maxUp){
//             MlinkUp[ni_N_link_up[ni_idx]] = li_idx;
//         } else {
//             sprintf(errmsg, "incoming links for NODE %s > Max allowed [api.c -> add_link_alt]", Node[ni_idx].ID);
//             api_report_writeErrorMsg(api_err_model_junctions, errmsg);
//             return api_err_model_junctions;
//         }
//         return 0;
//     } else {
//         ni_N_link_d[ni_idx] ++;
//         if (ni_N_link_dn[ni_idx] <= maxDn){
//             MlinkDn[ni_N_link_dn[ni_idx]] = li_idx;
//         } else {
//             sprintf(errmsg, "outgoing links for NODE %s > 3 [api.c -> add_link_alt]", Node[ni_idx].ID);
//             api_report_writeErrorMsg(api_err_model_junctions, errmsg);
//             return api_err_model_junctions;
//         }
//         return 0;    
//     }
//     api_report_writeErrorMsg(api_err_internal, "[api.c -> add_link_alt]");
//     return api_err_internal;
// }

//===============================================================================
int add_link(
    int li_idx,
    int ni_idx,
    int direction,
    int* ni_N_link_u,
    int* ni_Mlink_u1,
    int* ni_Mlink_u2,
    int* ni_Mlink_u3,
    int* ni_N_link_d,
    int* ni_Mlink_d1,
    int* ni_Mlink_d2,
    int* ni_Mlink_d3)
//===============================================================================
{
    if (direction == UPSTREAM) {
        ni_N_link_u[ni_idx] ++;
        if (ni_N_link_u[ni_idx] <= 3) {
            if (ni_N_link_u[ni_idx] == 1) {
                ni_Mlink_u1[ni_idx] = li_idx;
            } else if (ni_N_link_u[ni_idx] == 2) {
                ni_Mlink_u2[ni_idx] = li_idx;
            } else if (ni_N_link_u[ni_idx] == 3) {
                ni_Mlink_u3[ni_idx] = li_idx;
            } else {
                sprintf(errmsg, "ni_Mlink_u3 == 0 at NODE %s [api.c -> add_link]", Node[ni_idx].ID);
                api_report_writeErrorMsg(api_err_internal, "errmsg");
                return api_err_internal;
            }
            return 0;
        } else {
            sprintf(errmsg, "incoming links for NODE %s > 3 [api.c -> add_link]", Node[ni_idx].ID);
            api_report_writeErrorMsg(api_err_model_junctions, errmsg);
            return api_err_model_junctions;
        }
    } else {
        ni_N_link_d[ni_idx] ++;
        if (ni_N_link_d[ni_idx] <= 3) {
            if (ni_N_link_d[ni_idx] == 1) {
                ni_Mlink_d1[ni_idx] = li_idx;
            } else if (ni_N_link_d[ni_idx] == 2) {
                ni_Mlink_d2[ni_idx] = li_idx;
            } else if (ni_N_link_d[ni_idx] == 3) {
                ni_Mlink_d3[ni_idx] = li_idx;
            } else {
                sprintf(errmsg, "ni_Mlink_d3 == 0 at NODE %s [api.c -> add_link]", Node[ni_idx].ID);
                api_report_writeErrorMsg(api_err_internal, "errmsg");
                return api_err_internal;
            }
            return 0;
        } else {
            sprintf(errmsg, "outgoing links for NODE %s > 3 [api.c -> add_link]", Node[ni_idx].ID);
            api_report_writeErrorMsg(api_err_model_junctions, errmsg);
            return api_err_model_junctions;
        }
    }
    api_report_writeErrorMsg(api_err_internal, "[api.c -> add_link]");
    return api_err_internal;
}
//===============================================================================
int check_api_is_initialized(
    char * function_name)
//===============================================================================
    //  provides error if api has not been initialized
{
    if ( ErrorCode ) return error_getCode(ErrorCode);
    if ( !api->IsInitialized )
    {
        sprintf(errmsg, "[api.c -> %s -> check_api_is_initialized]", function_name);
        api_report_writeErrorMsg(api_err_not_initialized, errmsg);
        return api_err_not_initialized;
    }
    return 0;
}
//===============================================================================
int getTokens(
    char *ss)
//===============================================================================
// Copy pasted getTokens from src/input.c to ensure independence
// from the original EPA-SWMM code. In the original code
// getTokens is not defined as an external API function
//
//  Input:   s = a character string
//  Output:  returns number of tokens found in s
//  Purpose: scans a string for tokens, saving pointers to them
//           in shared variable Tok[].
//
//  Notes:   Tokens can be separated by the characters listed in SEPSTR
//           (spaces, tabs, newline, carriage return) which is defined
//           in CONSTS.H. Text between quotes is treated as a single token.
//
{
    int  len, mm, nn;
    char *cc;

    // --- begin with no tokens
    for (nn = 0; nn < MAXTOKS; nn++) Tok[nn] = NULL;
    nn = 0;

    // --- truncate s at start of comment
    cc = strchr(ss,';');
    if (cc) *cc = '\0';
    len = strlen(ss);

    // --- scan s for tokens until nothing left
    while (len > 0 && nn < MAXTOKS)
    {
        mm = strcspn(ss,SEPSTR);              // find token length
        if (mm == 0) ss++;                    // no token found
        else
        {
            if (*ss == '"')                  // token begins with quote
            {
                ss++;                        // start token after quote
                len--;                      // reduce length of s
                mm = strcspn(ss,"\"\n");      // find end quote or new line
            }
            ss[mm] = '\0';                    // null-terminate the token
            Tok[nn] = ss;                     // save pointer to token
            nn++;                            // update token count
            ss += mm+1;                       // begin next token
        }
        len -= mm+1;                         // update length of s
    }
    return(nn);
}
//===============================================================================
// EOF
//===============================================================================