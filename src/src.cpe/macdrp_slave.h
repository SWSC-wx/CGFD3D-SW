/*
 * macdrp_slave.h
 *
 *  Created on: Jan 29, 2018
 *      Author: bingo
 */

#ifndef SRC_MACDRP_SLAVE_H_
#define SRC_MACDRP_SLAVE_H_


#define BWD -1
#define FWD 1

/*
#define vec_L(var,i,FLAG) ((FLAG==1)?vec_LF(var,i):vec_LB(var,i))
#define vec_LF(var,i) ((-0.30874 * rDH)*var[i-1]+(-0.6326*rDH)*var[i]+(1.2330*rDH)*var[i+1]+(-0.3334*rDH)*var[i+2]+(0.04168*rDH)*var[i+3])
#define vec_LB(var,i) ((-0.04168 * rDH)*var[i-3]+(0.3334*rDH)*var[i-2]+(-1.2330*rDH)*var[i-1]+(0.6326*rDH)*var[i]+(0.30874*rDH)*var[i+1])
*/
#define vec_LL(var,i,FLAG) ((FLAG==1)?vec_LLF(var,i):vec_LLB(var,i))
#define vec_LLF(var,i) ((-0.30874*rDH)*var[i]+(-0.6326*rDH)*var[i+1]+(1.2330*rDH)*var[i+2]+(-0.3334*rDH)*var[i+3]+(0.04168*rDH)*var[i+4])
#define vec_LLB(var,i) ((-0.04168 * rDH)*var[i]+(0.3334*rDH)*var[i+1]+(-1.2330*rDH)*var[i+2]+(0.6326*rDH)*var[i+3]+(0.30874*rDH)*var[i+4])

#define L(var,idx,stride,FLAG) ((FLAG==1)?LF(var,idx,stride):LB(var,idx, stride))
#define LF(var,idx,stride) ((-0.30874*rDH)*var[idx - stride]+(-0.6326*rDH)*var[idx]+(1.2330*rDH)*var[idx + stride]+(-0.3334*rDH)*var[idx + 2 * stride]+(0.04168*rDH)*var[idx + 3 * stride])
#define LB(var,idx,stride) ((-0.04168*rDH)*var[idx - 3 * stride]+(0.3334*rDH)*var[idx - 2 * stride]+(-1.2330*rDH)*var[idx - stride]+(0.6326*rDH)*var[idx]+(0.30874*rDH)*var[idx + stride])


#define L22(var,idx,stride,FLAG) ((FLAG==1)?L22F(var,idx,stride):L22B(var,idx,stride))
#define L22F(var,idx,stride) (var[idx + stride]-var[idx])*rDH
#define L22B(var,idx,stride) (var[idx]-var[idx - stride])*rDH

#define L24(var,idx,stride,FLAG) ((FLAG==1)?L24F(var,idx,stride):L24B(var,idx,stride))
#define L24F(var,idx,stride) ((-7.0/6.0)*rDH*var[idx]+(8.0/6.0)*rDH*var[idx + stride]+(-1.0/6.0)*rDH*var[idx + 2 * stride])
#define L24B(var,idx,stride) ((1.0/6.0)*rDH*var[idx - 2 * stride]+(-8.0/6.0)*rDH*var[idx - stride]+(7.0/6.0)*rDH*var[idx])


#endif /* SRC_MACDRP_SLAVE_H_ */
