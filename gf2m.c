// GF(2^1279) = GF(2)[x]/(x^1279 + x^319 + x^127 + x^63 + 1)
#include <string.h>
#include <stdio.h>
#include "gf2m.h"

void ff_rnd(ff_t a) {
  FILE* urandom = fopen("/dev/urandom", "r");
  size_t bytes_read = fread(a, NUMWORDS, sizeof(uint64_t), urandom);
  fclose(urandom);
  a[NUMWORDS-1] &= 0x7FFFFFFFFFFFFFFFULL;
} 

void ff_print(ff_t a) {
  for(int i=NUMWORDS-1; i>=0; i--) printf("%016lX ", a[i]);
  printf("\n");
}

void ff_print_poly(ff_t a) {
  for(int i=0; i<NUMWORDS; i++)
    for(int j=0; j<64; j++)
      if(a[i] >> j & 1)
        printf("+x^%d", i*64+j);
  printf("\n");
}

// c = a*b. Don't clobber a or b.
void ff_mul(ff_t c, const ff_t a, const ff_t b) 
{ // initializing, assigning and copying data is done here
  ff_t a_copy,b_copy,reduce_poly;
  ff_copy(a_copy,a); 
  ff_copy(b_copy,b);
  ff_zero(reduce_poly); 
  reduce_poly[0] = 0x8000000000000001ULL; // assigning the bits for the reduction polynimial
  reduce_poly[1] = 0x8000000000000000ULL; // ie., x^1279 = x^319 + x^127 + x^63 + 1
  reduce_poly[4] = 0x8000000000000000ULL; // so assigning the respective bits to 1 is done here
  if (a_copy[0] & 0x0000000000000001ULL) { // checking weather a0 = 0, if true then copy b to c, else copy 0 to c
    ff_copy(c,b_copy);                       
  } 
  else {
    ff_zero(c);
  }    // loop from the bit_position 0 to m-1. Here m=1279
  for (int bit_position = 1; bit_position < 1279; bit_position++) { 

    uint64_t bits[20] = {0}; // initializing of temporary variable to handle bits on shifting
    for(int word_position=0; word_position< 20; word_position++) {
      if(b_copy[word_position]&0x8000000000000000ULL) // Checking the MSB
         bits[word_position]= 0x0000000000000001ULL; 
      b_copy[word_position] = b_copy[word_position] << 1; // shifting is done here
    }
    b_copy[1] |= bits[0]; // Replacement of Bits is done here
    b_copy[2] |= bits[1];
    b_copy[3] |= bits[2];
    b_copy[4] |= bits[3];
    b_copy[5] |= bits[4];
    b_copy[6] |= bits[5];
    b_copy[7] |= bits[6];
    b_copy[8] |= bits[7];
    b_copy[9] |= bits[8];
    b_copy[10] |= bits[9];
    b_copy[11] |= bits[10];
    b_copy[12] |= bits[11];
    b_copy[13] |= bits[12];
    b_copy[14] |= bits[13];
    b_copy[15] |= bits[14];
    b_copy[16] |= bits[15];
    b_copy[17] |= bits[16];
    b_copy[18] |= bits[17];
    b_copy[19] |= bits[18];
    b_copy[20] |= bits[19];

    if(b_copy[19]&0x8000000000000000ULL) { //checking weather MSB is set
      b_copy[19] &= 0x7FFFFFFFFFFFFFFFULL;
      // Reduction of the polynomial to a polynomial within the field is done here
      for(int word_position=0; word_position<20; word_position++) { 
        b_copy[word_position] ^= reduce_poly[word_position]; // reduction is done in order to bring the polynomial within the field
      }
     }

    if(a_copy[bit_position/64] & (0x0000000000000001ULL << (bit_position%64))) { //Checking weather the ith bit of a is 1, if true copy the reduced values 
     for (int word_position=0; word_position<20; word_position++) {
       c[word_position] ^= b_copy[word_position]; // copying the reduced values which are within the field now to c
     } 
    }
   }
 }
   // c = b^e. Don't clobber b or e.
void ff_exp(ff_t c, const ff_t b, const ff_t ee) {  
  ff_t a,g,e; // initialization and copying of the generator and exponent elements is done here
  ff_copy(g,b);
  ff_copy(e,ee);
  ff_zero(c); c[0] = 0x1; // initializing c to 1

  int i;
  for(i=1279; i>= 0; i--) {
          ff_mul(a,c,c); // from t to 0 multiplying a = c*c
          if (e[i/64] & (0x0000000000000001ULL << (i%64))) { //checking the ith bit of e for 1, if true c = a*g (generator element) else copy the output to c
                  ff_mul(c,a,g); 
          } else {                        
                  ff_copy(c,a);
          }
  }
}

