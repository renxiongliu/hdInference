#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    
#define     max(x,y)                (((x)>(y))?(x):(y))
#define     min(x,y)                (((x)<(y))?(x):(y))
#define     abs(x  )                (((x) > 0)?(x):(-x))
#define     sign(x )                (((x) > 0)?(1):(-1))
#define     zero(x )                (((x) == 0)?(1):(0))
    enum penalty_types
    {
        L1,
        Truncated_L1,
        MCP
    };
    
    
    typedef struct
    {
        double  *Delta, *Delta_old, *tmp1, *tmp2, *Lambda, *lam_mat;
        bool *constrained_idx_matrix;
        const int *zero_para_index, *num_of_zero_para;
    } tmpvars;
    
    
    
    
    
#ifdef __cplusplus
}
#endif
