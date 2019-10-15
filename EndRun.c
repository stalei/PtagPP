// © Shahram @ 2019
// Initial setup and configuration
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "GlobalVars.h"

void EndRun(int ERRcode, char file[100])
{
  //switch ERRcode{
	//case 0:
        //#ifdef DoParallel
        //if(ThisTask==0)
        //#endif
	printf("!!! fatal error in tagPP, file:%s - line:%d . I better stop!!! \n ",file,ERRcode);
	//}
exit(0);

}
