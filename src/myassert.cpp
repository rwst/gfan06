/*
 * myassert.cpp
 *
 *  Created on: Jul 10, 2011
 *      Author: anders
 */

#include "myassert.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#ifdef STACKDUMP_ENABLED
#include <execinfo.h>
//#include <cxxabi.h>

void dumpStackTrace()
{
  //link with -rdynamic??
  fprintf(stderr,"\nCALL STACK:\n");
  void *stackframes[100];
  int n=backtrace(stackframes,100);
  backtrace_symbols_fd(stackframes,n,2/*stderr*/);

 // fprintf(stderr,abi::__cxa_demangle(stack
}

/**
 * When the gcc compiler is called with
 * -D__assert_fail=__assert_fail2
 * then every failing assertion is tricked into calling the following function instead of the usual assertion fail.
 * This of course is implementation dependent.
 * Therefore, the above define should only be used together with -DSTACKDUMP_ENABLED .
 */


//extern void __assert_fail2 (__const char *__assertion, __const char *__file,
//                           unsigned int __line, __const char *__function)
extern void __assert_fail2 (__const char *__assertion, __const char *__file,
                           unsigned int __line, __const char *__function)throw()
{
  fprintf(stderr,"gfan:\n%s:%i:\n%s:\nAssertion `%s' failed.\n",__file,__line,__function,__assertion);
  dumpStackTrace();
//  *((char *)0)=0;//enable this to make gdb stop at assertions.
  exit(1);
}


#else

void dumpStackTrace()
{
  fprintf(stderr,"\nSTACK DUMPING NOT ENABLED. RECOMPILE myassert.cpp WITH stackdump=true\n");
}
#endif
