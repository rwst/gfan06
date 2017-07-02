/*
 * myassert.h
 *
 *  Created on: Jul 10, 2011
 *      Author: anders
 */

#ifndef MYASSERT_H_
#define MYASSERT_H_

/*
 * This header does not have to be included, unless an explicit dumpStackTrace call is needed.
 * Compiling (using gcc) with -DSTACKDUMP_ENABLED -D__assert_fail=__assert_fail2 (and linking the
 * .o file) will make the assertions fail with a stack dump automatically.
 */
void dumpStackTrace();


#endif /* MYASSERT_H_ */
