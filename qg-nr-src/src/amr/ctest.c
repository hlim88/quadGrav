#include <signal.h>

extern void fint2_();

/* void def_int2__()  for linux                */
/* void def_int2_()   for sgi                  */
/* void def_int2_()   for Intel C++ Compiler   */
/* void def_int2_()  for Mac?                 */
/*                    doesn't recognize sigset */
void def_int2_()
{
/*   sigset(2,fint2_); */
}
