#ifndef _NR_RANDOM_H
#define _NR_RANDOM_H

/**
* Long period (> 2 ¡Á 1018) random number generator of L¡¯Ecuyer with Bays-Durham shu?e and added safeguards. 
* @return a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). 
* @note Call with idum a negative integer to initialize; thereafter, do not alter
* idum between successive deviates in a sequence. RNMX should approximate the largest ?oating
* value that is less than 1.	
*/
float ran2(long *idum);

#endif

