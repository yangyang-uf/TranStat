/*================================
  Basic statistic summary functions
  ================================*/
#define functionx(name,type) name##_##type
#define function(name,type) functionx(name,type)

#define DATATYPE unsigned int
#define BASE uint
#include "statistics_source.h"
#undef BASE
#undef DATATYPE

#define DATATYPE unsigned long 
#define BASE ulong  
#include "statistics_source.h"
#undef BASE
#undef DATATYPE 

#define DATATYPE int 
#define BASE int  
#include "statistics_source.h"
#undef BASE
#undef DATATYPE

#define DATATYPE long 
#define BASE long
#include "statistics_source.h"
#undef BASE
#undef DATATYPE

#define DATATYPE float 
#define BASE float  
#include "statistics_source.h"
#undef BASE
#undef DATATYPE

#define DATATYPE double
#define BASE double  
#include "statistics_source.h"
#undef BASE
#undef DATATYPE 

#undef function
#undef functionx
