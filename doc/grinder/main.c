
#include "geom_physical_const.h"
#include "common_structures.h"

main() {
struct Scope Grid ;
void plot_template(char*,struct Scope,int,int);

/*  Parameters for construction of family of templates:
    obviosly these parameters should be read in from 
    a file at run time, but what the hell
*/

Grid.m_mn   = 0.8   ;      /* Max (total) mass allowed in search */
Grid.m_mx   = 50.0    ;    /* Min (total) mass allowed in search */
Grid.theta  = 0.978   ;    /* Sathya's theta */
Grid.dp     = 0.00321 ;    /* */
Grid.dq     = 0.0498  ;    /* */
Grid.f_start= 120.   ;    /* */

template_grid(&Grid);
plot_template("temp_list.ps",Grid,15,1);

return 0;
}
