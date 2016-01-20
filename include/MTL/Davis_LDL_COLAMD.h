//
// Math Template Library
//
// Copyright (c) 2015: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
//   1. Redistributions of source code must retain the above copyright notice, this list of
//      conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright notice, this list of
//      conditions and the following disclaimer in the documentation and/or other materials
//      provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
// WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

using MTL::I32;
using MTL::Min;
using MTL::Max;

//
// From LDL library by Timothy A. Davis.
//
static void LDLt_Symbolic
(int n,          /* A and L are n-by-n, where n >= 0 */
 const I32* Ap,  /* input of size n+1, not modified */
 const I32* Ai,  /* input of size nz=Ap[n], not modified */
 I32* Lp,        /* output of size n+1, not defined on input */
 I32* Parent,    /* output of size n, not defined on input */
 I32* Lnz,       /* output of size n, not defined on input */
 I32* Flag,      /* workspace of size n, not defn. on input or output */
 const I32* P,   /* optional input of size n */
 const I32* Pinv /* optional output of size n (used if P is not NULL) */
)
{
  int i, k, kk, p, p2;

  for (k = 0 ; k < n ; k++)
  {
    /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
    Parent[k] = -1;         /* parent of k is not yet known */
    Flag[k] = k;            /* mark node k as visited */
    Lnz[k] = 0;             /* count of nonzeros in column k of L */
    kk = P[k];
    p2 = Ap [kk+1] ;
    for (p = Ap[kk]; p < p2; p++)
    {
      /* A (i,k) is nonzero */
      i = Pinv[Ai[p]];
      if (i < k)
      {
        /* follow path from i to root of etree, stop at flagged node */
        for (; Flag[i] != k; i = Parent[i])
        {
          /* find parent of i if not yet determined */
          if (Parent[i] == -1)
            Parent[i] = k;

          Lnz[i]++;                         /* L (k,i) is nonzero */
          Flag[i] = k;                      /* mark i as visited */
        }
      }
    }
  }

  /* construct Lp index array from Lnz column counts */
  Lp[0] = 0;
  for (k = 0 ; k < n ; k++)
    Lp[k+1] = Lp[k] + Lnz[k];
}

template<class T>
static void LDLt_Numeric
(int n,              /* A and L are n-by-n, where n >= 0 */
 const I32* Ap,      /* input of size n+1, not modified */
 const I32* Ai,      /* input of size nz=Ap[n], not modified */
 const T* Ax,        /* input of size nz=Ap[n], not modified */
 const I32* Lp,      /* input of size n+1, not modified */
 const I32* Parent,  /* input of size n, not modified */
 I32* Lnz,           /* output of size n, not defn. on input */
 I32* Li,            /* output of size lnz=Lp[n], not defined on input */
 T* Lx,              /* output of size lnz=Lp[n], not defined on input */
 T* D,               /* output of size n, not defined on input */
 T* Y,               /* workspace of size n, not defn. on input or output */
 I32* Pattern,       /* workspace of size n, not defn. on input or output */
 I32* Flag ,         /* workspace of size n, not defn. on input or output */
 const I32* P,             /* optional input of size n */
 const I32* Pinv           /* optional output of size n (used if P is not NULL) */
)
{
  T yi, l_ki;
  int i, k, kk, p, p2, len, top ;
  for (k = 0 ; k < n ; k++)
  {
    /* compute nonzero Pattern of kth row of L, in topological order */
    Y[k] = 0.0;               /* Y(0:k) is now all zero */
    top = n;                  /* stack for pattern is empty */
    Flag[k] = k;              /* mark node k as visited */
    Lnz[k] = 0;               /* count of nonzeros in column k of L */
    kk = P[k];
    p2 = Ap[kk+1];
    for (p = Ap[kk] ; p < p2; p++)
    {
      i = Pinv[Ai[p]];   /* get A(i,k) */
      if (i <= k)
      {
        Y[i] += Ax[p] ;  /* scatter A(i,k) into Y (sum duplicates) */
        for (len = 0; Flag[i] != k; i = Parent[i])
        {
          Pattern[len++] = i;   /* L(k,i) is nonzero */
          Flag[i] = k;          /* mark i as visited */
        }

        while (len > 0)
          Pattern[--top] = Pattern[--len];
      }
    }

    /* compute numerical values kth row of L (a sparse triangular solve) */
    D[k] = Y[k];            /* get D(k,k) and clear Y(k) */
    Y[k] = 0.0;
    for (; top < n; top++)
    {
      i = Pattern[top];     /* Pattern [top:n-1] is pattern of L(:,k) */
      yi = Y[i];            /* get and clear Y(i) */
      Y[i] = 0.0;
      p2 = Lp[i] + Lnz[i];
      for (p = Lp[i]; p < p2; p++)
        Y [Li[p]] -= Lx[p] * yi;

      l_ki = yi / D[i];     /* the nonzero entry L(k,i) */
      D[k] -= l_ki * yi;
      Li[p] = k;            /* store L(k,i) in column form of L */
      Lx[p] = l_ki;
      Lnz[i]++;             /* increment count of nonzeros in col i */
    }
  }
}

//
// From COLAMD library by Stefan I. Larimore and Timothy A. Davis.
//
#define Int long
#define COLAMD_EMPTY   (-1)

/* Row and column status */
#define COLAMD_ALIVE   (0)
#define COLAMD_DEAD    (-1)

/* Column status */
#define COLAMD_DEAD_PRINCIPAL          (-1)
#define COLAMD_DEAD_NON_PRINCIPAL      (-2)

/* Macros for row and column status update and checking. */
#define ROW_IS_DEAD(r)                  ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)    (row_mark < COLAMD_ALIVE)
#define ROW_IS_ALIVE(r)                 (Row [r].shared2.mark >= COLAMD_ALIVE)
#define COL_IS_DEAD(c)                  (Col [c].start < COLAMD_ALIVE)
#define COL_IS_ALIVE(c)                 (Col [c].start >= COLAMD_ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)        (Col [c].start == COLAMD_DEAD_PRINCIPAL)
#define KILL_ROW(r)                     { Row [r].shared2.mark = COLAMD_DEAD ; }
#define KILL_PRINCIPAL_COL(c)           { Col [c].start = COLAMD_DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)       { Col [c].start = COLAMD_DEAD_NON_PRINCIPAL ; }

#define DENSE_DEGREE(alpha,n) \
  ((Int) Max (16.0, (alpha) * sqrt ((double) (n))))

#define ONES_COMPLEMENT(r) (-(r)-1)

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
#define COLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..6] are currently used. */
#define COLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define COLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define COLAMD_DENSE_COL 1

/* knobs [2]: aggressive absorption */
#define COLAMD_AGGRESSIVE 2

/* stats [2]: memory defragmentation count output statistic */
#define COLAMD_DEFRAG_COUNT 2

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define COLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */ 
#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6

/* error codes returned in stats [3]: */
#define COLAMD_OK                               (0)
#define COLAMD_OK_BUT_JUMBLED                   (1)
#define COLAMD_ERROR_A_not_present              (-1)
#define COLAMD_ERROR_p_not_present              (-2)
#define COLAMD_ERROR_nrow_negative              (-3)
#define COLAMD_ERROR_ncol_negative              (-4)
#define COLAMD_ERROR_nnz_negative               (-5)
#define COLAMD_ERROR_p0_nonzero                 (-6)
#define COLAMD_ERROR_A_too_small                (-7)
#define COLAMD_ERROR_col_length_negative        (-8)
#define COLAMD_ERROR_row_index_out_of_bounds    (-9)
#define COLAMD_ERROR_out_of_memory              (-10)
#define COLAMD_ERROR_internal_error             (-999)

typedef struct Colamd_Col_struct
{
    Int start ;         /* index for A of first row in this column, or DEAD */
                        /* if column is dead */
    Int length ;        /* number of rows in this column */
    union
    {
        Int thickness ; /* number of original columns represented by this */
                        /* col, if the column is alive */
        Int parent ;    /* parent in parent tree super-column structure, if */
                        /* the column is dead */
    } shared1 ;
    union
    {
        Int score ;     /* the score used to maintain heap, if col is alive */
        Int order ;     /* pivot ordering of this column, if col is dead */
    } shared2 ;
    union
    {
        Int headhash ;  /* head of a hash bucket, if col is at the head of */
                        /* a degree list */
        Int hash ;      /* hash value, if col is not in a degree list */
        Int prev ;      /* previous column in degree list, if col is in a */
                        /* degree list (but not at the head of a degree list) */
    } shared3 ;
    union
    {
        Int degree_next ;       /* next column, if col is in a degree list */
        Int hash_next ;         /* next column, if col is in a hash list */
    } shared4 ;

} Colamd_Col ;

typedef struct Colamd_Row_struct
{
    Int start ;         /* index for A of first col in this row */
    Int length ;        /* number of principal columns in this row */
    union
    {
        Int degree ;    /* number of principal & non-principal columns in row */
        Int p ;         /* used as a row pointer in init_rows_cols () */
    } shared1 ;
    union
    {
        Int mark ;      /* for computing set differences and marking dead rows*/
        Int first_column ;/* first column in row (used in garbage collection) */
    } shared2 ;

} Colamd_Row ;

/* add two values of type size_t, and check for integer overflow */
static size_t t_add (size_t a, size_t b, int *ok)
{
  (*ok) = (*ok) && ((a + b) >= Max (a,b)) ;
    return ((*ok) ? (a + b) : 0) ;
}

/* compute a*k where k is a small integer, and check for integer overflow */
static size_t t_mult (size_t a, size_t k, int *ok)
{
    size_t i, s = 0 ;
    for (i = 0 ; i < k ; i++)
    {
        s = t_add (s, a, ok) ;
    }
    return (s) ;
}

/* size of the Col and Row structures */
#define COLAMD_C(n_col,ok) \
    ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (Int)))

#define COLAMD_R(n_row,ok) \
    ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (Int)))

static void COLAMD_set_defaults
(
    /* === Parameters ======================================================= */
    double knobs [COLAMD_KNOBS]         /* knob array */
)
{
    /* === Local variables ================================================== */
    Int i ;

    if (!knobs)
    {
        return ;                        /* no knobs to initialize */
    }
    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
        knobs [i] = 0 ;
    }
    knobs [COLAMD_DENSE_ROW] = 10 ;
    knobs [COLAMD_DENSE_COL] = 10 ;
    knobs [COLAMD_AGGRESSIVE] = 1 ;     /* default: do aggressive absorption*/
}

static size_t COLAMD_recommended        /* returns recommended value of Alen. */
(
    /* === Parameters ======================================================= */
    Int nnz,                    /* number of nonzeros in A */
    Int n_row,                  /* number of rows in A */
    Int n_col                   /* number of columns in A */
)
{
    size_t s, c, r ;
    int ok = 1 ;
    if (nnz < 0 || n_row < 0 || n_col < 0)
    {
        return (0) ;
    }
    s = t_mult (nnz, 2, &ok) ;      /* 2*nnz */
    c = COLAMD_C (n_col, &ok) ;     /* size of column structures */
    r = COLAMD_R (n_row, &ok) ;     /* size of row structures */
    s = t_add (s, c, &ok) ;
    s = t_add (s, r, &ok) ;
    s = t_add (s, n_col, &ok) ;     /* elbow room */
    s = t_add (s, nnz/5, &ok) ;     /* elbow room */
    ok = ok && (s < INT_MAX) ;
    return (ok ? s : 0) ;
}

static Int init_rows_cols
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int p [],
    Int stats [COLAMD_STATS]
) ;

static void init_scoring
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,
    Int *p_n_col2,
    Int *p_max_deg
) ;

static Int find_ordering
(
    Int n_row,
    Int n_col,
    Int Alen,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int n_col2,
    Int max_deg,
    Int pfree,
    Int aggressive
) ;

static void order_children
(
    Int n_col,
    Colamd_Col Col [],
    Int p []
) ;

/* ========================================================================== */
/* === colamd =============================================================== */
/* ========================================================================== */

/*
    The colamd routine computes a column ordering Q of a sparse matrix
    A such that the LU factorization P(AQ) = LU remains sparse, where P is
    selected via partial pivoting.   The routine can also be viewed as
    providing a permutation Q such that the Cholesky factorization
    (AQ)'(AQ) = LL' remains sparse.
*/

static Int colamd               /* returns TRUE if successful, FALSE otherwise*/
(
    /* === Parameters ======================================================= */

    Int n_row,                  /* number of rows in A */
    Int n_col,                  /* number of columns in A */
    Int Alen,                   /* length of A */
    Int A [],                   /* row indices of A */
    Int p [],                   /* pointers to columns in A */
    double knobs [COLAMD_KNOBS],/* parameters (uses defaults if NULL) */
    Int stats [COLAMD_STATS]    /* output statistics and error codes */
)
{
    /* === Local variables ================================================== */

    Int i ;                     /* loop index */
    Int nnz ;                   /* nonzeros in A */
    size_t Row_size ;           /* size of Row [], in integers */
    size_t Col_size ;           /* size of Col [], in integers */
    size_t need ;               /* minimum required length of A */
    Colamd_Row *Row ;           /* pointer into A of Row [0..n_row] array */
    Colamd_Col *Col ;           /* pointer into A of Col [0..n_col] array */
    Int n_col2 ;                /* number of non-dense, non-empty columns */
    Int n_row2 ;                /* number of non-dense, non-empty rows */
    Int ngarbage ;              /* number of garbage collections performed */
    Int max_deg ;               /* maximum row degree */
    double default_knobs [COLAMD_KNOBS] ;       /* default knobs array */
    Int aggressive ;            /* do aggressive absorption */
    int ok ;

    /* === Check the input arguments ======================================== */

    assert(stats);
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
      stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    assert(A);
    assert(p);

    if (n_row < 0)      /* n_row must be >= 0 */
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
      stats [COLAMD_INFO1] = n_row ;
      return 0;
    }

    if (n_col < 0)      /* n_col must be >= 0 */
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
      stats [COLAMD_INFO1] = n_col ;
      return 0;
    }

    nnz = p [n_col] ;
    if (nnz < 0)        /* nnz must be >= 0 */
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
      stats [COLAMD_INFO1] = nnz ;
      return 0;
    }

    if (p [0] != 0)
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero ;
      stats [COLAMD_INFO1] = p [0] ;
      return 0;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
        COLAMD_set_defaults (default_knobs) ;
        knobs = default_knobs ;
    }

    aggressive = (knobs [COLAMD_AGGRESSIVE] != 0) ;

    /* === Allocate the Row and Col arrays from array A ===================== */

    ok = 1 ;
    Col_size = COLAMD_C (n_col, &ok) ;      /* size of Col array of structs */
    Row_size = COLAMD_R (n_row, &ok) ;      /* size of Row array of structs */

    /* need = 2*nnz + n_col + Col_size + Row_size ; */
    need = t_mult (nnz, 2, &ok) ;
    need = t_add (need, n_col, &ok) ;
    need = t_add (need, Col_size, &ok) ;
    need = t_add (need, Row_size, &ok) ;

    if (!ok || need > (size_t) Alen || need > INT_MAX)
    {
        /* not enough space in array A to perform the ordering */
        stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
        stats [COLAMD_INFO1] = (Int)need ;
        stats [COLAMD_INFO2] = Alen ;
        return 0;
    }

    Alen -= Int(Col_size + Row_size);
    Col = (Colamd_Col *) &A [Alen] ;
    Row = (Colamd_Row *) &A [Alen + Col_size] ;

    /* === Construct the row and column data structures ===================== */

    if (!init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
    {
        /* input matrix is invalid */
        return 0;
    }

    /* === Initialize scores, kill dense rows/columns ======================= */

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
        &n_row2, &n_col2, &max_deg) ;

    /* === Order the supercolumns =========================================== */

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
        n_col2, max_deg, 2*nnz, aggressive) ;

    /* === Order the non-principal columns ================================== */

    order_children (n_col, Col, p) ;

    /* === Return statistics in stats ======================================= */

    stats [COLAMD_DENSE_ROW] = n_row - n_row2 ;
    stats [COLAMD_DENSE_COL] = n_col - n_col2 ;
    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
    return 1;
}

/* ========================================================================== */
/* === symamd =============================================================== */
/* ========================================================================== */
static Int symamd                       /* return TRUE if OK, FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n,                              /* number of rows and columns of A */
    const Int A [],                     /* row indices of A */
    const Int p [],                     /* column pointers of A */
    Int perm [],                        /* output permutation, size n+1 */
    double knobs [COLAMD_KNOBS],        /* parameters (uses defaults if NULL) */
    Int stats [COLAMD_STATS],           /* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
                                        /* pointer to calloc (ANSI C) or */
                                        /* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
                                        /* pointer to free (ANSI C) or */
                                        /* mxFree (for MATLAB mexFunction) */
)
{
    /* === Local variables ================================================== */

    Int *count ;                /* length of each column of M, and col pointer*/
    Int *mark ;                 /* mark array for finding duplicate entries */
    Int *M ;                    /* row indices of matrix M */
    size_t Mlen ;               /* length of M */
    Int n_row ;                 /* number of rows in M */
    Int nnz ;                   /* number of entries in A */
    Int i ;                     /* row index of A */
    Int j ;                     /* column index of A */
    Int k ;                     /* row index of M */ 
    Int mnz ;                   /* number of nonzeros in M */
    Int pp ;                    /* index into a column of A */
    Int last_row ;              /* last row seen in the current column */
    Int length ;                /* number of nonzeros in a column */

    double cknobs [COLAMD_KNOBS] ;              /* knobs for colamd */
    double default_knobs [COLAMD_KNOBS] ;       /* default knobs for colamd */

    /* === Check the input arguments ======================================== */
    assert(stats);

    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
        stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    assert(A);
    assert(p);

    if (n < 0)          /* n must be >= 0 */
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
      stats [COLAMD_INFO1] = n ;
      return 0;
    }

    nnz = p [n] ;
    if (nnz < 0)        /* nnz must be >= 0 */
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
      stats [COLAMD_INFO1] = nnz ;
      0;
    }

    if (p [0] != 0)
    {
      stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero ;
      stats [COLAMD_INFO1] = p [0] ;
      return 0;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
      COLAMD_set_defaults (default_knobs) ;
      knobs = default_knobs ;
    }

    /* === Allocate count and mark ========================================== */

    count = (Int *) ((*allocate) (n+1, sizeof (Int))) ;
    if (!count)
    {
        stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
        return 0 ;
    }

    mark = (Int *) ((*allocate) (n+1, sizeof (Int))) ;
    if (!mark)
    {
        stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
        (*release) ((void *) count) ;
        return 0 ;
    }

    /* === Compute column counts of M, check if A is valid ================== */

    stats [COLAMD_INFO3] = 0 ;  /* number of duplicate or unsorted row indices*/

    for (i = 0 ; i < n ; i++)
    {
        mark [i] = -1 ;
    }

    for (j = 0 ; j < n ; j++)
    {
        last_row = -1 ;

        length = p [j+1] - p [j] ;
        if (length < 0)
        {
            /* column pointers must be non-decreasing */
            stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
            stats [COLAMD_INFO1] = j ;
            stats [COLAMD_INFO2] = length ;
            (*release) ((void *) count) ;
            (*release) ((void *) mark) ;
            return 0 ;
        }

        for (pp = p [j] ; pp < p [j+1] ; pp++)
        {
            i = A [pp] ;
            if (i < 0 || i >= n)
            {
                /* row index i, in column j, is out of bounds */
                stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
                stats [COLAMD_INFO1] = j ;
                stats [COLAMD_INFO2] = i ;
                stats [COLAMD_INFO3] = n ;
                (*release) ((void *) count) ;
                (*release) ((void *) mark) ;
                return 0 ;
            }

            if (i <= last_row || mark [i] == j)
            {
                /* row index is unsorted or repeated (or both), thus col */
                /* is jumbled.  This is a notice, not an error condition. */
                stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
                stats [COLAMD_INFO1] = j ;
                stats [COLAMD_INFO2] = i ;
                (stats [COLAMD_INFO3]) ++ ;
            }

            if (i > j && mark [i] != j)
            {
                /* row k of M will contain column indices i and j */
                count [i]++ ;
                count [j]++ ;
            }

            /* mark the row as having been seen in this column */
            mark [i] = j ;

            last_row = i ;
        }
    }

    /* v2.4: removed free(mark) */

    /* === Compute column pointers of M ===================================== */

    /* use output permutation, perm, for column pointers of M */
    perm [0] = 0 ;
    for (j = 1 ; j <= n ; j++)
    {
        perm [j] = perm [j-1] + count [j-1] ;
    }
    for (j = 0 ; j < n ; j++)
    {
        count [j] = perm [j] ;
    }

    /* === Construct M ====================================================== */

    mnz = perm [n] ;
    n_row = mnz / 2 ;
    Mlen = COLAMD_recommended (mnz, n_row, n) ;
    M = (Int *) ((*allocate) (Mlen, sizeof (Int))) ;

    if (!M)
    {
        stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
        (*release) ((void *) count) ;
        (*release) ((void *) mark) ;
        return 0 ;
    }

    k = 0 ;

    if (stats [COLAMD_STATUS] == COLAMD_OK)
    {
        /* Matrix is OK */
        for (j = 0 ; j < n ; j++)
        {
            assert (p [j+1] - p [j] >= 0) ;
            for (pp = p [j] ; pp < p [j+1] ; pp++)
            {
                i = A [pp] ;
                assert (i >= 0 && i < n) ;
                if (i > j)
                {
                    /* row k of M contains column indices i and j */
                    M [count [i]++] = k ;
                    M [count [j]++] = k ;
                    k++ ;
                }
            }
        }
    }
    else
    {
        /* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
        for (i = 0 ; i < n ; i++)
        {
            mark [i] = -1 ;
        }
        for (j = 0 ; j < n ; j++)
        {
            assert (p [j+1] - p [j] >= 0) ;
            for (pp = p [j] ; pp < p [j+1] ; pp++)
            {
                i = A [pp] ;
                assert (i >= 0 && i < n) ;
                if (i > j && mark [i] != j)
                {
                    /* row k of M contains column indices i and j */
                    M [count [i]++] = k ;
                    M [count [j]++] = k ;
                    k++ ;

                    mark [i] = j ;
                }
            }
        }
        /* v2.4: free(mark) moved below */
    }

    /* count and mark no longer needed */
    (*release) ((void *) count) ;
    (*release) ((void *) mark) ;        /* v2.4: free (mark) moved here */
    assert (k == n_row) ;

    /* === Adjust the knobs for M =========================================== */

    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
        cknobs [i] = knobs [i] ;
    }

    /* there are no dense rows in M */
    cknobs [COLAMD_DENSE_ROW] = -1 ;
    cknobs [COLAMD_DENSE_COL] = knobs [COLAMD_DENSE_ROW] ;

    /* === Order the columns of M =========================================== */

    /* v2.4: colamd cannot fail here, so the error check is removed */
    (void) colamd (n_row, n, (Int) Mlen, M, perm, cknobs, stats) ;

    /* Note that the output permutation is now in perm */

    /* === get the statistics for symamd from colamd ======================== */

    /* a dense column in colamd means a dense row and col in symamd */
    stats [COLAMD_DENSE_ROW] = stats [COLAMD_DENSE_COL] ;

    /* === Free M =========================================================== */

    (*release) ((void *) M) ;
    return 1 ;

}

/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

/*
    Takes the column form of the matrix in A and creates the row form of the
    matrix.  Also, row and column attributes are stored in the Col and Row
    structs.  If the columns are un-sorted or contain duplicate row indices,
    this routine will also sort and remove duplicate row indices from the
    column form of the matrix.  Returns FALSE if the matrix is invalid,
    TRUE otherwise.  Not user-callable.
*/

static Int init_rows_cols       /* returns TRUE if OK, or FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n_row,                  /* number of rows of A */
    Int n_col,                  /* number of columns of A */
    Colamd_Row Row [],          /* of size n_row+1 */
    Colamd_Col Col [],          /* of size n_col+1 */
    Int A [],                   /* row indices of A, of size Alen */
    Int p [],                   /* pointers to columns in A, of size n_col+1 */
    Int stats [COLAMD_STATS]    /* colamd statistics */ 
)
{
    /* === Local variables ================================================== */

    Int col ;                   /* a column index */
    Int row ;                   /* a row index */
    Int *cp ;                   /* a column pointer */
    Int *cp_end ;               /* a pointer to the end of a column */
    Int *rp ;                   /* a row pointer */
    Int *rp_end ;               /* a pointer to the end of a row */
    Int last_row ;              /* previous row */

    /* === Initialize columns, and check column pointers ==================== */

    for (col = 0 ; col < n_col ; col++)
    {
        Col [col].start = p [col] ;
        Col [col].length = p [col+1] - p [col] ;

        if (Col [col].length < 0)
        {
            /* column pointers must be non-decreasing */
            stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
            stats [COLAMD_INFO1] = col ;
            stats [COLAMD_INFO2] = Col [col].length ;
            return 0 ;
        }

        Col [col].shared1.thickness = 1 ;
        Col [col].shared2.score = 0 ;
        Col [col].shared3.prev = COLAMD_EMPTY ;
        Col [col].shared4.degree_next = COLAMD_EMPTY ;
    }

    /* p [0..n_col] no longer needed, used as "head" in subsequent routines */

    /* === Scan columns, compute row degrees, and check row indices ========= */

    stats [COLAMD_INFO3] = 0 ;  /* number of duplicate or unsorted row indices*/

    for (row = 0 ; row < n_row ; row++)
    {
        Row [row].length = 0 ;
        Row [row].shared2.mark = -1 ;
    }

    for (col = 0 ; col < n_col ; col++)
    {
        last_row = -1 ;

        cp = &A [p [col]] ;
        cp_end = &A [p [col+1]] ;

        while (cp < cp_end)
        {
            row = *cp++ ;

            /* make sure row indices within range */
            if (row < 0 || row >= n_row)
            {
                stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
                stats [COLAMD_INFO1] = col ;
                stats [COLAMD_INFO2] = row ;
                stats [COLAMD_INFO3] = n_row ;
                return 0 ;
            }

            if (row <= last_row || Row [row].shared2.mark == col)
            {
                /* row index are unsorted or repeated (or both), thus col */
                /* is jumbled.  This is a notice, not an error condition. */
                stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
                stats [COLAMD_INFO1] = col ;
                stats [COLAMD_INFO2] = row ;
                (stats [COLAMD_INFO3]) ++ ;
            }

            if (Row [row].shared2.mark != col)
            {
                Row [row].length++ ;
            }
            else
            {
                /* this is a repeated entry in the column, */
                /* it will be removed */
                Col [col].length-- ;
            }

            /* mark the row as having been seen in this column */
            Row [row].shared2.mark = col ;

            last_row = row ;
        }
    }

    /* === Compute row pointers ============================================= */

    /* row form of the matrix starts directly after the column */
    /* form of matrix in A */
    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    Row [0].shared2.mark = -1 ;
    for (row = 1 ; row < n_row ; row++)
    {
        Row [row].start = Row [row-1].start + Row [row-1].length ;
        Row [row].shared1.p = Row [row].start ;
        Row [row].shared2.mark = -1 ;
    }

    /* === Create row form ================================================== */

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
        /* if cols jumbled, watch for repeated row indices */
        for (col = 0 ; col < n_col ; col++)
        {
            cp = &A [p [col]] ;
            cp_end = &A [p [col+1]] ;
            while (cp < cp_end)
            {
                row = *cp++ ;
                if (Row [row].shared2.mark != col)
                {
                    A [(Row [row].shared1.p)++] = col ;
                    Row [row].shared2.mark = col ;
                }
            }
        }
    }
    else
    {
        /* if cols not jumbled, we don't need the mark (this is faster) */
        for (col = 0 ; col < n_col ; col++)
        {
            cp = &A [p [col]] ;
            cp_end = &A [p [col+1]] ;
            while (cp < cp_end)
            {
                A [(Row [*cp++].shared1.p)++] = col ;
            }
        }
    }

    /* === Clear the row marks and set row degrees ========================== */

    for (row = 0 ; row < n_row ; row++)
    {
        Row [row].shared2.mark = 0 ;
        Row [row].shared1.degree = Row [row].length ;
    }

    /* === See if we need to re-create columns ============================== */

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
        /* === Compute col pointers ========================================= */

        /* col form of the matrix starts at A [0]. */
        /* Note, we may have a gap between the col form and the row */
        /* form if there were duplicate entries, if so, it will be */
        /* removed upon the first garbage collection */
        Col [0].start = 0 ;
        p [0] = Col [0].start ;
        for (col = 1 ; col < n_col ; col++)
        {
            /* note that the lengths here are for pruned columns, i.e. */
            /* no duplicate row indices will exist for these columns */
            Col [col].start = Col [col-1].start + Col [col-1].length ;
            p [col] = Col [col].start ;
        }

        /* === Re-create col form =========================================== */

        for (row = 0 ; row < n_row ; row++)
        {
            rp = &A [Row [row].start] ;
            rp_end = rp + Row [row].length ;
            while (rp < rp_end)
            {
                A [(p [*rp++])++] = row ;
            }
        }
    }

    /* === Done.  Matrix is not (or no longer) jumbled ====================== */

    return 1;
}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

/*
    Kills dense or empty columns and rows, calculates an initial score for
    each column, and places all columns in the degree lists.  Not user-callable.
*/

static void init_scoring
(
    /* === Parameters ======================================================= */

    Int n_row,                  /* number of rows of A */
    Int n_col,                  /* number of columns of A */
    Colamd_Row Row [],          /* of size n_row+1 */
    Colamd_Col Col [],          /* of size n_col+1 */
    Int A [],                   /* column form and row form of A */
    Int head [],                /* of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameters */
    Int *p_n_row2,              /* number of non-dense, non-empty rows */
    Int *p_n_col2,              /* number of non-dense, non-empty columns */
    Int *p_max_deg              /* maximum row degree */
)
{
    /* === Local variables ================================================== */

    Int c ;                     /* a column index */
    Int r, row ;                /* a row index */
    Int *cp ;                   /* a column pointer */
    Int deg ;                   /* degree of a row or column */
    Int *cp_end ;               /* a pointer to the end of a column */
    Int *new_cp ;               /* new column pointer */
    Int col_length ;            /* length of pruned column */
    Int score ;                 /* current column score */
    Int n_col2 ;                /* number of non-dense, non-empty columns */
    Int n_row2 ;                /* number of non-dense, non-empty rows */
    Int dense_row_count ;       /* remove rows with more entries than this */
    Int dense_col_count ;       /* remove cols with more entries than this */
    Int min_score ;             /* smallest column score */
    Int max_deg ;               /* maximum row degree */
    Int next_col ;              /* Used to add to degree list.*/

    /* === Extract knobs ==================================================== */

    /* Note: if knobs contains a NaN, this is undefined: */
    if (knobs [COLAMD_DENSE_ROW] < 0)
    {
        /* only remove completely dense rows */
        dense_row_count = n_col-1 ;
    }
    else
    {
        dense_row_count = DENSE_DEGREE (knobs [COLAMD_DENSE_ROW], n_col) ;
    }
    if (knobs [COLAMD_DENSE_COL] < 0)
    {
        /* only remove completely dense columns */
        dense_col_count = n_row-1 ;
    }
    else
    {
        dense_col_count =
          DENSE_DEGREE (knobs [COLAMD_DENSE_COL], Min (n_row, n_col)) ;
    }

    max_deg = 0 ;
    n_col2 = n_col ;
    n_row2 = n_row ;

    /* === Kill empty columns =============================================== */

    /* Put the empty columns at the end in their natural order, so that LU */
    /* factorization can proceed as far as possible. */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
        deg = Col [c].length ;
        if (deg == 0)
        {
            /* this is a empty column, kill and order it last */
            Col [c].shared2.order = --n_col2 ;
            KILL_PRINCIPAL_COL (c) ;
        }
    }

    /* === Kill dense columns =============================================== */

    /* Put the dense columns at the end, in their natural order */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
        /* skip any dead columns */
        if (COL_IS_DEAD (c))
        {
            continue ;
        }
        deg = Col [c].length ;
        if (deg > dense_col_count)
        {
            /* this is a dense column, kill and order it last */
            Col [c].shared2.order = --n_col2 ;
            /* decrement the row degrees */
            cp = &A [Col [c].start] ;
            cp_end = cp + Col [c].length ;
            while (cp < cp_end)
            {
                Row [*cp++].shared1.degree-- ;
            }
            KILL_PRINCIPAL_COL (c) ;
        }
    }

    /* === Kill dense and empty rows ======================================== */

    for (r = 0 ; r < n_row ; r++)
    {
        deg = Row [r].shared1.degree ;
        assert (deg >= 0 && deg <= n_col) ;
        if (deg > dense_row_count || deg == 0)
        {
            /* kill a dense or empty row */
            KILL_ROW (r) ;
            --n_row2 ;
        }
        else
        {
            /* keep track of max degree of remaining rows */
          max_deg = Max (max_deg, deg) ;
        }
    }

    /* === Compute initial column scores ==================================== */

    /* At this point the row degrees are accurate.  They reflect the number */
    /* of "live" (non-dense) columns in each row.  No empty rows exist. */
    /* Some "live" columns may contain only dead rows, however.  These are */
    /* pruned in the code below. */

    /* now find the initial matlab score for each column */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
        /* skip dead column */
        if (COL_IS_DEAD (c))
        {
            continue ;
        }
        score = 0 ;
        cp = &A [Col [c].start] ;
        new_cp = cp ;
        cp_end = cp + Col [c].length ;
        while (cp < cp_end)
        {
            /* get a row */
            row = *cp++ ;
            /* skip if dead */
            if (ROW_IS_DEAD (row))
            {
                continue ;
            }
            /* compact the column */
            *new_cp++ = row ;
            /* add row's external degree */
            score += Row [row].shared1.degree - 1 ;
            /* guard against integer overflow */
            score = Min (score, n_col) ;
        }
        /* determine pruned column length */
        col_length = (Int) (new_cp - &A [Col [c].start]) ;
        if (col_length == 0)
        {
            /* a newly-made null column (all rows in this col are "dense" */
            /* and have already been killed) */
            Col [c].shared2.order = --n_col2 ;
            KILL_PRINCIPAL_COL (c) ;
        }
        else
        {
            /* set column length and set score */
            assert (score >= 0) ;
            assert (score <= n_col) ;
            Col [c].length = col_length ;
            Col [c].shared2.score = score ;
        }
    }

    /* At this point, all empty rows and columns are dead.  All live columns */
    /* are "clean" (containing no dead rows) and simplicial (no supercolumns */
    /* yet).  Rows may contain dead columns, but all live rows contain at */
    /* least one live column. */

    /* === Initialize degree lists ========================================== */

    /* clear the hash buckets */
    for (c = 0 ; c <= n_col ; c++)
    {
        head [c] = COLAMD_EMPTY ;
    }
    min_score = n_col ;
    /* place in reverse order, so low column indices are at the front */
    /* of the lists.  This is to encourage natural tie-breaking */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
        /* only add principal columns to degree lists */
        if (COL_IS_ALIVE (c))
        {

            /* === Add columns score to DList =============================== */

            score = Col [c].shared2.score ;

            assert (min_score >= 0) ;
            assert (min_score <= n_col) ;
            assert (score >= 0) ;
            assert (score <= n_col) ;
            assert (head [score] >= COLAMD_EMPTY) ;

            /* now add this column to dList at proper score location */
            next_col = head [score] ;
            Col [c].shared3.prev = COLAMD_EMPTY ;
            Col [c].shared4.degree_next = next_col ;

            /* if there already was a column with the same score, set its */
            /* previous pointer to this new column */
            if (next_col != COLAMD_EMPTY)
            {
                Col [next_col].shared3.prev = c ;
            }
            head [score] = c ;

            /* see if this score is less than current min */
            min_score = Min (min_score, score) ;
        }
    }

    /* === Return number of remaining columns, and max row degree =========== */

    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;
}

/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

/*
    Clears the Row [].shared2.mark array, and returns the new tag_mark.
    Return value is the new tag_mark.  Not user-callable.
*/

static Int clear_mark   /* return the new value for tag_mark */
(
    /* === Parameters ======================================================= */

    Int tag_mark,       /* new value of tag_mark */
    Int max_mark,       /* max allowed value of tag_mark */

    Int n_row,          /* number of rows in A */
    Colamd_Row Row []   /* Row [0 ... n_row-1].shared2.mark is set to zero */
)
{
    /* === Local variables ================================================== */

    Int r ;

    if (tag_mark <= 0 || tag_mark >= max_mark)
    {
        for (r = 0 ; r < n_row ; r++)
        {
            if (ROW_IS_ALIVE (r))
            {
                Row [r].shared2.mark = 0 ;
            }
        }
        tag_mark = 1 ;
    }

    return (tag_mark) ;
}

/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

/*
    Defragments and compacts columns and rows in the workspace A.  Used when
    all avaliable memory has been used while performing row merging.  Returns
    the index of the first free position in A, after garbage collection.  The
    time taken by this routine is linear is the size of the array A, which is
    itself linear in the number of nonzeros in the input matrix.
    Not user-callable.
*/

static Int garbage_collection  /* returns the new value of pfree */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows */
    Int n_col,			/* number of columns */
    Colamd_Row Row [],		/* row info */
    Colamd_Col Col [],		/* column info */
    Int A [],			/* A [0 ... Alen-1] holds the matrix */
    Int *pfree			/* &A [0] ... pfree is in use */
)
{
    /* === Local variables ================================================== */

    Int *psrc ;			/* source pointer */
    Int *pdest ;		/* destination pointer */
    Int j ;			/* counter */
    Int r ;			/* a row index */
    Int c ;			/* a column index */
    Int length ;		/* length of a row or column */

    /* === Defragment the columns =========================================== */

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    /* move and compact the column */
	    assert (pdest <= psrc) ;
	    Col [c].start = (Int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (Int) (pdest - &A [Col [c].start]) ;
	}
    }

    /* === Prepare to defragment the rows =================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_DEAD (r) || (Row [r].length == 0))
	{
	    /* This row is already dead, or is of zero length.  Cannot compact
	     * a row of zero length, so kill it.  NOTE: in the current version,
	     * there are no zero-length live rows.  Kill the row (for the first
	     * time, or again) just to be safe. */
	    KILL_ROW (r) ;
	}
	else
	{
	    /* save first column index in Row [r].shared2.first_column */
	    psrc = &A [Row [r].start] ;
	    Row [r].shared2.first_column = *psrc ;
	    assert (ROW_IS_ALIVE (r)) ;
	    /* flag the start of the row with the one's complement of row */
	    *psrc = ONES_COMPLEMENT (r) ;
	}
    }

    /* === Defragment the rows ============================================== */

    psrc = pdest ;
    while (psrc < pfree)
    {
	/* find a negative number ... the start of a row */
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    /* get the row index */
	    r = ONES_COMPLEMENT (*psrc) ;
	    assert (r >= 0 && r < n_row) ;
	    /* restore first column index */
	    *psrc = Row [r].shared2.first_column ;
	    assert (ROW_IS_ALIVE (r)) ;
	    assert (Row [r].length > 0) ;
	    /* move and compact the row */
	    assert (pdest <= psrc) ;
	    Row [r].start = (Int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (Int) (pdest - &A [Row [r].start]) ;
	    assert (Row [r].length > 0) ;
	}
    }
    /* ensure we found all the rows */
    //assert (debug_rows == 0) ;

    /* === Return the new value of pfree ==================================== */

    return ((Int) (pdest - &A [0])) ;
}

/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

/*
    Detects supercolumns by finding matches between columns in the hash buckets.
    Check amongst columns in the set A [row_start ... row_start + row_length-1].
    The columns under consideration are currently *not* in the degree lists,
    and have already been placed in the hash buckets.

    The hash bucket for columns whose hash function is equal to h is stored
    as follows:

	if head [h] is >= 0, then head [h] contains a degree list, so:

		head [h] is the first column in degree bucket h.
		Col [head [h]].headhash gives the first column in hash bucket h.

	otherwise, the degree list is empty, and:

		-(head [h] + 2) is the first column in hash bucket h.

    For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
    column" pointer.  Col [c].shared3.hash is used instead as the hash number
    for that column.  The value of Col [c].shared4.hash_next is the next column
    in the same hash bucket.

    Assuming no, or "few" hash collisions, the time taken by this routine is
    linear in the sum of the sizes (lengths) of each column whose score has
    just been computed in the approximate degree computation.
    Not user-callable.
*/

static void detect_super_cols
(
    /* === Parameters ======================================================= */

    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A */
    Int head [],		/* head of degree lists and hash buckets */
    Int row_start,		/* pointer to set of columns to check */
    Int row_length		/* number of columns to check */
)
{
    /* === Local variables ================================================== */

    Int hash ;			/* hash value for a column */
    Int *rp ;			/* pointer to a row */
    Int c ;			/* a column index */
    Int super_c ;		/* column index of the column to absorb into */
    Int *cp1 ;			/* column pointer for column super_c */
    Int *cp2 ;			/* column pointer for column c */
    Int length ;		/* length of column super_c */
    Int prev_c ;		/* column preceding c in hash bucket */
    Int i ;			/* loop counter */
    Int *rp_end ;		/* pointer to the end of the row */
    Int col ;			/* a column index in the row to check */
    Int head_column ;		/* first column in hash bucket or degree list */
    Int first_col ;		/* first column in hash bucket */

    /* === Consider each column in the row ================================== */

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	/* get hash number for this column */
	hash = Col [col].shared3.hash ;
	//assert (hash <= n_col) ;

	/* === Get the first column in this hash bucket ===================== */

	head_column = head [hash] ;
	if (head_column > COLAMD_EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}

	/* === Consider each column in the hash bucket ====================== */

	for (super_c = first_col ; super_c != COLAMD_EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    assert (COL_IS_ALIVE (super_c)) ;
	    assert (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    /* prev_c is the column preceding column c in the hash bucket */
	    prev_c = super_c ;

	    /* === Compare super_c with all columns after it ================ */

	    for (c = Col [super_c].shared4.hash_next ;
		 c != COLAMD_EMPTY ; c = Col [c].shared4.hash_next)
	    {
		assert (c != super_c) ;
		assert (COL_IS_ALIVE (c)) ;
		assert (Col [c].shared3.hash == hash) ;

		/* not identical if lengths or scores are different */
		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score)
		{
		    prev_c = c ;
		    continue ;
		}

		/* compare the two columns */
		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    /* the columns are "clean" (no dead rows) */
		    //assert (ROW_IS_ALIVE (*cp1))  ;
		    //assert (ROW_IS_ALIVE (*cp2))  ;
		    /* row indices will same order for both supercols, */
		    /* no gather scatter nessasary */
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		/* the two columns are different if the for-loop "broke" */
		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		/* === Got it!  two columns are identical =================== */

		assert (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;
		/* order c later, in order_children() */
		Col [c].shared2.order = COLAMD_EMPTY ;
		/* remove c from hash bucket */
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;
	    }
	}

	/* === Empty this hash bucket ======================================= */

	if (head_column > COLAMD_EMPTY)
	{
	    /* corresponding degree list "hash" is not empty */
	    Col [head_column].shared3.headhash = COLAMD_EMPTY ;
	}
	else
	{
	    /* corresponding degree list "hash" is empty */
	    head [hash] = COLAMD_EMPTY ;
	}
    }
}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

/*
    Order the principal columns of the supercolumn form of the matrix
    (no supercolumns on input).  Uses a minimum approximate column minimum
    degree ordering method.  Not user-callable.
*/

static Int find_ordering       /* return the number of garbage collections */
(
    /* === Parameters ======================================================= */

    Int n_row,                  /* number of rows of A */
    Int n_col,                  /* number of columns of A */
    Int Alen,                   /* size of A, 2*nnz + n_col or larger */
    Colamd_Row Row [],          /* of size n_row+1 */
    Colamd_Col Col [],          /* of size n_col+1 */
    Int A [],                   /* column form and row form of A */
    Int head [],                /* of size n_col+1 */
    Int n_col2,                 /* Remaining columns to order */
    Int max_deg,                /* Maximum row degree */
    Int pfree,                  /* index of first free slot (2*nnz on entry) */
    Int aggressive
)
{
    /* === Local variables ================================================== */

    Int k ;                     /* current pivot ordering step */
    Int pivot_col ;             /* current pivot column */
    Int *cp ;                   /* a column pointer */
    Int *rp ;                   /* a row pointer */
    Int pivot_row ;             /* current pivot row */
    Int *new_cp ;               /* modified column pointer */
    Int *new_rp ;               /* modified row pointer */
    Int pivot_row_start ;       /* pointer to start of pivot row */
    Int pivot_row_degree ;      /* number of columns in pivot row */
    Int pivot_row_length ;      /* number of supercolumns in pivot row */
    Int pivot_col_score ;       /* score of pivot column */
    Int needed_memory ;         /* free space needed for pivot row */
    Int *cp_end ;               /* pointer to the end of a column */
    Int *rp_end ;               /* pointer to the end of a row */
    Int row ;                   /* a row index */
    Int col ;                   /* a column index */
    Int max_score ;             /* maximum possible score */
    Int cur_score ;             /* score of current column */
    unsigned Int hash ;         /* hash value for supernode detection */
    Int head_column ;           /* head of hash bucket */
    Int first_col ;             /* first column in hash bucket */
    Int tag_mark ;              /* marker value for mark array */
    Int row_mark ;              /* Row [row].shared2.mark */
    Int set_difference ;        /* set difference size of row with pivot row */
    Int min_score ;             /* smallest column score */
    Int col_thickness ;         /* "thickness" (no. of columns in a supercol) */
    Int max_mark ;              /* maximum value of tag_mark */
    Int pivot_col_thickness ;   /* number of columns represented by pivot col */
    Int prev_col ;              /* Used by Dlist operations. */
    Int next_col ;              /* Used by Dlist operations. */
    Int ngarbage ;              /* number of garbage collections performed */

    /* === Initialization and clear mark ==================================== */

    max_mark = INT_MAX - n_col ;        /* INT_MAX defined in <limits.h> */
    tag_mark = clear_mark (0, max_mark, n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;

    /* === Order the columns ================================================ */

    for (k = 0 ; k < n_col2 ; /* 'k' is incremented below */)
    {
        /* === Select pivot column, and order it ============================ */

        /* make sure degree list isn't empty */
        assert (min_score >= 0) ;
        assert (min_score <= n_col) ;
        assert (head [min_score] >= COLAMD_EMPTY) ;

        /* get pivot column from head of minimum degree list */
        while (head [min_score] == COLAMD_EMPTY && min_score < n_col)
        {
            min_score++ ;
        }
        pivot_col = head [min_score] ;
        assert (pivot_col >= 0 && pivot_col <= n_col) ;
        next_col = Col [pivot_col].shared4.degree_next ;
        head [min_score] = next_col ;
        if (next_col != COLAMD_EMPTY)
        {
            Col [next_col].shared3.prev = COLAMD_EMPTY ;
        }

        assert (COL_IS_ALIVE (pivot_col)) ;

        /* remember score for defrag check */
        pivot_col_score = Col [pivot_col].shared2.score ;

        /* the pivot column is the kth column in the pivot order */
        Col [pivot_col].shared2.order = k ;

        /* increment order count by column thickness */
        pivot_col_thickness = Col [pivot_col].shared1.thickness ;
        k += pivot_col_thickness ;
        assert (pivot_col_thickness > 0) ;

        /* === Garbage_collection, if necessary ============================= */

        needed_memory = Min (pivot_col_score, n_col - k) ;
        if (pfree + needed_memory >= Alen)
        {
            pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
            ngarbage++ ;
            /* after garbage collection we will have enough */
            assert (pfree + needed_memory < Alen) ;
            /* garbage collection has wiped out the Row[].shared2.mark array */
            tag_mark = clear_mark (0, max_mark, n_row, Row) ;
        }

        /* === Compute pivot row pattern ==================================== */

        /* get starting location for this new merged row */
        pivot_row_start = pfree ;

        /* initialize new row counts to zero */
        pivot_row_degree = 0 ;

        /* tag pivot column as having been visited so it isn't included */
        /* in merged pivot row */
        Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

        /* pivot row is the union of all rows in the pivot column pattern */
        cp = &A [Col [pivot_col].start] ;
        cp_end = cp + Col [pivot_col].length ;
        while (cp < cp_end)
        {
            /* get a row */
            row = *cp++ ;
            /* skip if row is dead */
            if (ROW_IS_ALIVE (row))
            {
                rp = &A [Row [row].start] ;
                rp_end = rp + Row [row].length ;
                while (rp < rp_end)
                {
                    /* get a column */
                    col = *rp++ ;
                    /* add the column, if alive and untagged */
                    col_thickness = Col [col].shared1.thickness ;
                    if (col_thickness > 0 && COL_IS_ALIVE (col))
                    {
                        /* tag column in pivot row */
                        Col [col].shared1.thickness = -col_thickness ;
                        assert (pfree < Alen) ;
                        /* place column in pivot row */
                        A [pfree++] = col ;
                        pivot_row_degree += col_thickness ;
                    }
                }
            }
        }

        /* clear tag on pivot column */
        Col [pivot_col].shared1.thickness = pivot_col_thickness ;
        max_deg = Max (max_deg, pivot_row_degree) ;

        /* === Kill all rows used to construct pivot row ==================== */

        /* also kill pivot row, temporarily */
        cp = &A [Col [pivot_col].start] ;
        cp_end = cp + Col [pivot_col].length ;
        while (cp < cp_end)
        {
            /* may be killing an already dead row */
            row = *cp++ ;
            KILL_ROW (row) ;
        }

        /* === Select a row index to use as the new pivot row =============== */

        pivot_row_length = pfree - pivot_row_start ;
        if (pivot_row_length > 0)
        {
            /* pick the "pivot" row arbitrarily (first row in col) */
            pivot_row = A [Col [pivot_col].start] ;
        }
        else
        {
            /* there is no pivot row, since it is of zero length */
            pivot_row = COLAMD_EMPTY ;
            assert (pivot_row_length == 0) ;
        }
        assert (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

        /* === Approximate degree computation =============================== */

        /* Here begins the computation of the approximate degree.  The column */
        /* score is the sum of the pivot row "length", plus the size of the */
        /* set differences of each row in the column minus the pattern of the */
        /* pivot row itself.  The column ("thickness") itself is also */
        /* excluded from the column score (we thus use an approximate */
        /* external degree). */

        /* The time taken by the following code (compute set differences, and */
        /* add them up) is proportional to the size of the data structure */
        /* being scanned - that is, the sum of the sizes of each column in */
        /* the pivot row.  Thus, the amortized time to compute a column score */
        /* is proportional to the size of that column (where size, in this */
        /* context, is the column "length", or the number of row indices */
        /* in that column).  The number of row indices in a column is */
        /* monotonically non-decreasing, from the length of the original */
        /* column on input to colamd. */

        /* === Compute set differences ====================================== */

        /* pivot row is currently dead - it will be revived later. */

        /* for each column in pivot row */
        rp = &A [pivot_row_start] ;
        rp_end = rp + pivot_row_length ;
        while (rp < rp_end)
        {
            col = *rp++ ;
            assert (COL_IS_ALIVE (col) && col != pivot_col) ;

            /* clear tags used to construct pivot row pattern */
            col_thickness = -Col [col].shared1.thickness ;
            assert (col_thickness > 0) ;
            Col [col].shared1.thickness = col_thickness ;

            /* === Remove column from degree list =========================== */

            cur_score = Col [col].shared2.score ;
            prev_col = Col [col].shared3.prev ;
            next_col = Col [col].shared4.degree_next ;
            assert (cur_score >= 0) ;
            assert (cur_score <= n_col) ;
            assert (cur_score >= COLAMD_EMPTY) ;
            if (prev_col == COLAMD_EMPTY)
            {
                head [cur_score] = next_col ;
            }
            else
            {
                Col [prev_col].shared4.degree_next = next_col ;
            }
            if (next_col != COLAMD_EMPTY)
            {
                Col [next_col].shared3.prev = prev_col ;
            }

            /* === Scan the column ========================================== */

            cp = &A [Col [col].start] ;
            cp_end = cp + Col [col].length ;
            while (cp < cp_end)
            {
                /* get a row */
                row = *cp++ ;
                row_mark = Row [row].shared2.mark ;
                /* skip if dead */
                if (ROW_IS_MARKED_DEAD (row_mark))
                {
                    continue ;
                }
                assert (row != pivot_row) ;
                set_difference = row_mark - tag_mark ;
                /* check if the row has been seen yet */
                if (set_difference < 0)
                {
                    assert (Row [row].shared1.degree <= max_deg) ;
                    set_difference = Row [row].shared1.degree ;
                }
                /* subtract column thickness from this row's set difference */
                set_difference -= col_thickness ;
                assert (set_difference >= 0) ;
                /* absorb this row if the set difference becomes zero */
                if (set_difference == 0 && aggressive)
                {
                    KILL_ROW (row) ;
                }
                else
                {
                    /* save the new mark */
                    Row [row].shared2.mark = set_difference + tag_mark ;
                }
            }
        }

        /* === Add up set differences for each column ======================= */

        /* for each column in pivot row */
        rp = &A [pivot_row_start] ;
        rp_end = rp + pivot_row_length ;
        while (rp < rp_end)
        {
            /* get a column */
            col = *rp++ ;
            assert (COL_IS_ALIVE (col) && col != pivot_col) ;
            hash = 0 ;
            cur_score = 0 ;
            cp = &A [Col [col].start] ;
            /* compact the column */
            new_cp = cp ;
            cp_end = cp + Col [col].length ;

            while (cp < cp_end)
            {
                /* get a row */
                row = *cp++ ;
                assert(row >= 0 && row < n_row) ;
                row_mark = Row [row].shared2.mark ;
                /* skip if dead */
                if (ROW_IS_MARKED_DEAD (row_mark))
                {
                    continue ;
                }
                assert (row_mark >= tag_mark) ;
                /* compact the column */
                *new_cp++ = row ;
                /* compute hash function */
                hash += row ;
                /* add set difference */
                cur_score += row_mark - tag_mark ;
                /* integer overflow... */
                cur_score = Min (cur_score, n_col) ;
            }

            /* recompute the column's length */
            Col [col].length = (Int) (new_cp - &A [Col [col].start]) ;

            /* === Further mass elimination ================================= */

            if (Col [col].length == 0)
            {
                /* nothing left but the pivot row in this column */
                KILL_PRINCIPAL_COL (col) ;
                pivot_row_degree -= Col [col].shared1.thickness ;
                assert (pivot_row_degree >= 0) ;
                /* order it */
                Col [col].shared2.order = k ;
                /* increment order count by column thickness */
                k += Col [col].shared1.thickness ;
            }
            else
            {
                /* === Prepare for supercolumn detection ==================== */

                /* save score so far */
                Col [col].shared2.score = cur_score ;

                /* add column to hash table, for supercolumn detection */
                hash %= n_col + 1 ;

                assert (((Int) hash) <= n_col) ;

                head_column = head [hash] ;
                if (head_column > COLAMD_EMPTY)
                {
                    /* degree list "hash" is non-empty, use prev (shared3) of */
                    /* first column in degree list as head of hash bucket */
                    first_col = Col [head_column].shared3.headhash ;
                    Col [head_column].shared3.headhash = col ;
                }
                else
                {
                    /* degree list "hash" is empty, use head as hash bucket */
                    first_col = - (head_column + 2) ;
                    head [hash] = - (col + 2) ;
                }
                Col [col].shared4.hash_next = first_col ;

                /* save hash function in Col [col].shared3.hash */
                Col [col].shared3.hash = (Int) hash ;
                assert (COL_IS_ALIVE (col)) ;
            }
        }

        /* The approximate external column degree is now computed.  */

        /* === Supercolumn detection ======================================== */

        detect_super_cols (

                Col, A, head, pivot_row_start, pivot_row_length) ;

        /* === Kill the pivotal column ====================================== */

        KILL_PRINCIPAL_COL (pivot_col) ;

        /* === Clear mark =================================================== */

        tag_mark = clear_mark (tag_mark+max_deg+1, max_mark, n_row, Row) ;

        /* === Finalize the new pivot row, and column scores ================ */

        /* for each column in pivot row */
        rp = &A [pivot_row_start] ;
        /* compact the pivot row */
        new_rp = rp ;
        rp_end = rp + pivot_row_length ;
        while (rp < rp_end)
        {
            col = *rp++ ;
            /* skip dead columns */
            if (COL_IS_DEAD (col))
            {
                continue ;
            }
            *new_rp++ = col ;
            /* add new pivot row to column */
            A [Col [col].start + (Col [col].length++)] = pivot_row ;

            /* retrieve score so far and add on pivot row's degree. */
            /* (we wait until here for this in case the pivot */
            /* row's degree was reduced due to mass elimination). */
            cur_score = Col [col].shared2.score + pivot_row_degree ;

            /* calculate the max possible score as the number of */
            /* external columns minus the 'k' value minus the */
            /* columns thickness */
            max_score = n_col - k - Col [col].shared1.thickness ;

            /* make the score the external degree of the union-of-rows */
            cur_score -= Col [col].shared1.thickness ;

            /* make sure score is less or equal than the max score */
            cur_score = Min (cur_score, max_score) ;
            assert (cur_score >= 0) ;

            /* store updated score */
            Col [col].shared2.score = cur_score ;

            /* === Place column back in degree list ========================= */

            assert (min_score >= 0) ;
            assert (min_score <= n_col) ;
            assert (cur_score >= 0) ;
            assert (cur_score <= n_col) ;
            assert (head [cur_score] >= COLAMD_EMPTY) ;
            next_col = head [cur_score] ;
            Col [col].shared4.degree_next = next_col ;
            Col [col].shared3.prev = COLAMD_EMPTY ;
            if (next_col != COLAMD_EMPTY)
            {
                Col [next_col].shared3.prev = col ;
            }
            head [cur_score] = col ;

            /* see if this score is less than current min */
            min_score = Min (min_score, cur_score) ;

        }

        /* === Resurrect the new pivot row ================================== */

        if (pivot_row_degree > 0)
        {
            /* update pivot row length to reflect any cols that were killed */
            /* during super-col detection and mass elimination */
            Row [pivot_row].start  = pivot_row_start ;
            Row [pivot_row].length = (Int) (new_rp - &A[pivot_row_start]) ;
            assert (Row [pivot_row].length > 0) ;
            Row [pivot_row].shared1.degree = pivot_row_degree ;
            Row [pivot_row].shared2.mark = 0 ;
            /* pivot row is no longer dead */
        }
    }

    /* === All principal columns have now been ordered ====================== */

    return (ngarbage) ;
}


/* ========================================================================== */
/* === order_children ======================================================= */
/* ========================================================================== */

/*
    The find_ordering routine has ordered all of the principal columns (the
    representatives of the supercolumns).  The non-principal columns have not
    yet been ordered.  This routine orders those columns by walking up the
    parent tree (a column is a child of the column which absorbed it).  The
    final permutation vector is then placed in p [0 ... n_col-1], with p [0]
    being the first column, and p [n_col-1] being the last.  It doesn't look
    like it at first glance, but be assured that this routine takes time linear
    in the number of columns.  Although not immediately obvious, the time
    taken by this routine is O (n_col), that is, linear in the number of
    columns.  Not user-callable.
*/

static void order_children
(
    /* === Parameters ======================================================= */

    Int n_col,                  /* number of columns of A */
    Colamd_Col Col [],          /* of size n_col+1 */
    Int p []                    /* p [0 ... n_col-1] is the column permutation*/
)
{
    /* === Local variables ================================================== */

    Int i ;                     /* loop counter for all columns */
    Int c ;                     /* column index */
    Int parent ;                /* index of column's parent */
    Int order ;                 /* column's order */

    /* === Order each non-principal column ================================== */

    for (i = 0 ; i < n_col ; i++)
    {
        /* find an un-ordered non-principal column */
        assert (COL_IS_DEAD (i)) ;
        if (!COL_IS_DEAD_PRINCIPAL (i) && Col [i].shared2.order == COLAMD_EMPTY)
        {
            parent = i ;
            /* once found, find its principal parent */
            do
            {
                parent = Col [parent].shared1.parent ;
            } while (!COL_IS_DEAD_PRINCIPAL (parent)) ;

            /* now, order all un-ordered non-principal columns along path */
            /* to this parent.  collapse tree at the same time */
            c = i ;
            /* get order of parent */
            order = Col [parent].shared2.order ;

            do
            {
                assert (Col [c].shared2.order == COLAMD_EMPTY) ;

                /* order this column */
                Col [c].shared2.order = order++ ;
                /* collaps tree */
                Col [c].shared1.parent = parent ;

                /* get immediate parent of this column */
                c = Col [c].shared1.parent ;

                /* continue until we hit an ordered column.  There are */
                /* guarranteed not to be anymore unordered columns */
                /* above an ordered column */
            } while (Col [c].shared2.order == COLAMD_EMPTY) ;

            /* re-order the super_col parent to largest order for this group */
            Col [parent].shared2.order = order ;
        }
    }

    /* === Generate the permutation ========================================= */

    for (c = 0 ; c < n_col ; c++)
    {
        p [Col [c].shared2.order] = c ;
    }
}
