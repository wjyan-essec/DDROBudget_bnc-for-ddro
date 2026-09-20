cdef extern from "cplex.h":

    #########
    # Types #
    #########

    ctypedef void* CPXENVptr
    ctypedef void* CPXCENVptr
    ctypedef void* CPXLPptr
    ctypedef void* CPXCLPptr
    ctypedef int (*LazyCallbackCFunc)(CPXENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
    ctypedef int (*UserCutCallbackCFunc)(CPXENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
    ctypedef int (*BranchCallbackCFunc)(CPXCENVptr xenv, void *cbdata, int wherefrom, void *cbhandle, int brtype, int brset, int nodecnt, int bdcnt, const int *nodebeg, const int *xindex, const char *lu, const double *bd, const double *nodeest, int *useraction_p)
    ctypedef int (*NodeCallbackCFunc)(CPXCENVptr xenv, void *cbdata, int wherefrom, void *cbhandle, int *nodeindex, int *useraction)
    ctypedef int (*IncumbentCallbackCFunc)(CPXCENVptr xenv, void *cbdata, int wherefrom, void *cbhandle, double objval, double *x, int *isfeas_p, int *useraction_p)

    #############
    # Functions #
    #############

    # Create environment
    CPXENVptr CPXopenCPLEX(int *status_p)
    int CPXcloseCPLEX(CPXENVptr *env_p)
    
    # Create model
    CPXLPptr CPXcreateprob(CPXENVptr env, int *status_p, char *name)
    int CPXfreeprob(CPXENVptr env, CPXLPptr * lp_p)

    # Read model
    int CPXgetnumcols(CPXCENVptr env, CPXLPptr lp)
    int CPXgetnumrows(CPXCENVptr env, CPXLPptr lp)
    int CPXgetx (CPXCENVptr env, CPXCLPptr lp, double *x, int begin, int end)
    int CPXgetsense (CPXCENVptr env, CPXCLPptr lp, char *sense, int begin, int end)
    int CPXgetlb (CPXCENVptr env, CPXCLPptr lp, double *lb, int begin, int end)
    int CPXgetub (CPXCENVptr env, CPXCLPptr lp, double *ub, int begin, int end)
    int CPXgetbase (CPXCENVptr env, CPXCLPptr lp, int *cstat, int *rstat)
    int CPXgetbhead (CPXCENVptr env, CPXCLPptr lp, int *head, double *x)
    int CPXbinvacol (CPXCENVptr env, CPXCLPptr lp, int j, double *x)
    int CPXbinvarow (CPXCENVptr env, CPXCLPptr lp, int j, double *x)
    int CPXbinvcol (CPXCENVptr env, CPXCLPptr lp, int j, double *x)
    int CPXbinvrow (CPXCENVptr env, CPXCLPptr lp, int j, double *x)
    int CPXgetrhs( CPXCENVptr env, CPXCLPptr lp, double * rhs, int begin, int end )
    int CPXgetrows( CPXCENVptr env, CPXCLPptr lp, int * nzcnt_p, int * rmatbeg, int * rmatind, double * rmatval, int rmatspace, int * surplus_p, int begin, int end )

    # Statistics
    int CPXgettime ( CPXCENVptr env , double * timestamp_p )
    int CPXgetnumcuts( CPXCENVptr env, CPXCLPptr lp, int cuttype, int * num_p )
    int CPXgetmiprelgap( CPXCENVptr env, CPXCLPptr lp, double * gap_p )
    int CPXgetnodecnt( CPXCENVptr env, CPXCLPptr lp )
    int CPXgetsolnpoolnumsolns( CPXCENVptr env, CPXCLPptr lp )
    int CPXgetbestobjval( CPXCENVptr env, CPXCLPptr lp, double * objval_p )
    int CPXgetobjval( CPXCENVptr env, CPXCLPptr lp, double * objval_p )
    int CPXgetstat( CPXCENVptr env, CPXCLPptr lp )
    
    # Build model
    int CPXchgobjsen (CPXCENVptr env, CPXLPptr lp, int maxormin)
    int CPXnewcols(CPXENVptr env, CPXLPptr lp, int ccnt, const double * obj, const double * lb, const double * ub, const char * xctype, char ** colname)
    int CPXaddrows(CPXENVptr env, CPXLPptr lp, int ccnt, int rcnt, int nzcnt, const double * rhs, const char * sense, const int * rmatbeg, const int * rmatind, const double * rmatval, char ** colname, char ** rowname)
    int CPXcopyorder( CPXCENVptr env, CPXLPptr lp, int cnt, const int * indices, const int * priority, const int * direction )
    int CPXcopybase(CPXCENVptr env, CPXLPptr lp, const int * cstat, const int * rstat )
    
    # I/O functions
    int CPXwriteprob(CPXENVptr env, CPXLPptr lp, const char *filename_str, const char *filetype)
    int CPXordwrite(CPXCENVptr env, CPXCLPptr lp, const char * filename_str)

    # Parameters
    int CPXsetintparam(CPXENVptr env, int whichparam, int newvalue)
    int CPXsetdblparam(CPXENVptr env, int whichparam, double newvalue)
    double CPXgetdblparam(CPXENVptr env, int whichparam, double* value)

    # Solve
    int CPXmipopt(CPXENVptr env, CPXLPptr lp)

    # Lazy Constraint Callbacks
    int CPXsetlazyconstraintcallbackfunc(CPXENVptr env, LazyCallbackCFunc lazyconcallback, void* cbhandle)
    int CPXgetcallbacknodelp (CPXCENVptr env, void *cbdata, int wherefrom, CPXLPptr *nodelp_p)
    int CPXcutcallbackadd( CPXCENVptr env, void * cbdata, int wherefrom, int nzcnt, double rhs, int sense, const int * cutind, const double* cutval, int purgeable )
    int CPXcutcallbackaddlocal( CPXCENVptr env, void * cbdata, int wherefrom, int nzcnt, double rhs, int sense, const int * cutind, const double * cutval )
    
    # User Cut Callbacks
    int CPXsetusercutcallbackfunc(CPXENVptr env, UserCutCallbackCFunc cutcb, void* cbhandle)
    
    # Branch Callbacks
    int CPXsetbranchcallbackfunc(CPXENVptr env, BranchCallbackCFunc branchcallback, void * cbhandle)

    # Node Callbacks
    int CPXsetnodecallbackfunc(CPXENVptr env, NodeCallbackCFunc nodecallback, void * cbhandle)

    # Incumbent Callbacks
    int CPXsetincumbentcallbackfunc(CPXENVptr env, IncumbentCallbackCFunc incumbentcallback, void * cbhandle)

    #############
    # Constants #
    #############

    double CPX_INFBOUND

    # Parameter Names
    int CPX_PARAM_SCRIND
    int CPX_PARAM_STARTALG
    int CPX_PARAM_SUBALG
    int CPX_PARAM_PREIND
    int CPX_PARAM_TILIM
    int CPX_PARAM_THREADS
    int CPX_PARAM_EPINT
    int CPX_PARAM_MIPCBREDLP
    int CPX_PARAM_HEURFREQ
    int CPX_PARAM_CUTPASS
    int CPX_PARAM_MIPDISPLAY
    int CPX_PARAM_MIPSEARCH
    int CPX_PARAM_PRELINEAR
    int CPX_PARAM_MIPORDIND
    int CPX_PARAM_PROBE
    int CPX_PARAM_REDUCE
    int CPX_PARAM_VARSEL
    int CPX_PARAM_NUMERICALEMPHASIS
    int CPX_PARAM_TUNINGTILIM
    int CPX_PARAM_CLOCKTYPE
    int CPX_PARAM_EPGAP

    # Parameter Values
    int CPX_ALG_DUAL
    int CPX_ALG_PRIMAL
    int CPX_ON
    int CPX_OFF
    int CPX_MIPSEARCH_TRADITIONAL
    int CPX_PREREDUCE_NOPRIMALORDUAL
    int CPX_VARSEL_MININFEAS
    int CPX_VARSEL_MAXINFEAS
    int CPX_VARSEL_DEFAULT
    int CPX_VARSEL_PSEUDO
    int CPX_VARSEL_STRONG
    int CPX_VARSEL_PSEUDOREDUCED

    # Objective Senses
    int CPX_MIN
    int CPX_MAX

    # Basis status
    int CPX_AT_LOWER
    int CPX_BASIC
    int CPX_AT_UPPER
    int CPX_FREE_SUPER
    
    # Cut Purgeability
    # https://www.ibm.com/docs/en/cofz/12.9.0?topic=cpxxcutcallbackadd-cpxcutcallbackadd
    int CPX_USECUT_FORCE
    int CPX_USECUT_PURGE
    int CPX_USECUT_FILTER

    # Callback Return codes
    int CPX_CALLBACK_DEFAULT
    int CPX_CALLBACK_FAIL
    int CPX_CALLBACK_SET
    int CPX_CALLBACK_ABORT_CUT_LOOP

    # Cut Types
    # https://www.ibm.com/docs/en/cofz/12.9.0?topic=cpxxgetnumcuts-cpxgetnumcuts
    int CPX_CUT_COVER
    int CPX_CUT_GUBCOVER
    int CPX_CUT_FLOWCOVER
    int CPX_CUT_CLIQUE
    int CPX_CUT_FRAC
    int CPX_CUT_MIR
    int CPX_CUT_FLOWPATH
    int CPX_CUT_DISJ
    int CPX_CUT_IMPLBD
    int CPX_CUT_ZEROHALF
    int CPX_CUT_MCF
    int CPX_CUT_LANDP
    int CPX_CUT_USER
    int CPX_CUT_TABLE
    int CPX_CUT_SOLNPOOL
    int CPX_CUT_LOCALIMPLBD
    int CPX_CUT_BQP
    int CPX_CUT_RLT
    int CPX_CUT_BENDERS

    # MIP Status
    # https://www.ibm.com/docs/en/icos/22.1.1?topic=micclcarm-solution-status-symbols-specific-mip-in-cplex-callable-library-c-api
    int CPXMIP_ABORT_FEAS
    int CPXMIP_OPTIMAL
    int CPXMIP_ABORT_INFEAS
    int CPXMIP_ABORT_RELAXATION_UNBOUNDED
    int CPXMIP_ABORT_RELAXED
    int CPXMIP_DETTIME_LIM_FEAS
    int CPXMIP_DETTIME_LIM_INFEAS
    int CPXMIP_FAIL_FEAS
    int CPXMIP_FAIL_FEAS_NO_TREE
    int CPXMIP_FAIL_INFEAS
    int CPXMIP_FAIL_INFEAS_NO_TREE
    int CPXMIP_FEASIBLE
    int CPXMIP_FEASIBLE_RELAXED_INF
    int CPXMIP_FEASIBLE_RELAXED_QUAD
    int CPXMIP_FEASIBLE_RELAXED_SUM
    int CPXMIP_INFEASIBLE
    int CPXMIP_INForUNBD
    int CPXMIP_MEM_LIM_FEAS
    int CPXMIP_MEM_LIM_INFEAS
    int CPXMIP_NODE_LIM_FEAS
    int CPXMIP_NODE_LIM_INFEAS
    int CPXMIP_OPTIMAL_INFEAS
    int CPXMIP_OPTIMAL_POPULATED
    int CPXMIP_OPTIMAL_POPULATED_TOL
    int CPXMIP_OPTIMAL_RELAXED_INF
    int CPXMIP_OPTIMAL_RELAXED_QUAD
    int CPXMIP_OPTIMAL_RELAXED_SUM
    int CPXMIP_OPTIMAL_TOL
    int CPXMIP_POPULATESOL_LIM
    int CPXMIP_SOL_LIM
    int CPXMIP_TIME_LIM_FEAS
    int CPXMIP_TIME_LIM_INFEAS
    int CPXMIP_UNBOUNDED

    # Error codes
    int CPXERR_NEGATIVE_SURPLUS
