/*
   ************************************************
	Author: Ali Raza
	Date: 01/15/2019.
    Crypto Scheme:
	An Efficient Predicate Encryption with
	Constant Pairing Computations and Minimum Costs
	https://ieeexplore.ieee.org/document/7399733
   ************************************************

   Implemented on Type-3 pairing

   Compile with modules as specified below

   For MR_PAIRING_CP curve
   cp_pair.cpp zzn2.cpp big.cpp zzn.cpp ecn.cpp miracl.lib

   For MR_PAIRING_MNT curve
   mnt_pair.cpp zzn6a.cpp ecn3.cpp zzn3.cpp zzn2.cpp big.cpp zzn.cpp ecn.cpp miracl.lib
	
   For MR_PAIRING_BN curve
   bn_pair.cpp zzn12a.cpp ecn2.cpp zzn4.cpp zzn2.cpp big.cpp zzn.cpp ecn.cpp miracl.lib

   For MR_PAIRING_KSS curve
   kss_pair.cpp zzn18.cpp zzn6.cpp ecn3.cpp zzn3.cpp big.cpp zzn.cpp ecn.cpp miracl.lib

   For MR_PAIRING_BLS curve
   bls_pair.cpp zzn24.cpp zzn8.cpp zzn4.cpp zzn2.cpp ecn4.cpp big.cpp zzn.cpp ecn.cpp miracl.lib

	Test Program.
*/

#include <iostream>
#include <ctime>

//********* choose just one of these pairs **********
//#define MR_PAIRING_CP      // AES-80 security   
//#define AES_SECURITY 80

//#define MR_PAIRING_MNT	// AES-80 security
//#define AES_SECURITY 80

#define MR_PAIRING_BN    // AES-128 or AES-192 security
#define AES_SECURITY 128
//#define AES_SECURITY 192

//#define MR_PAIRING_KSS    // AES-192 security
//#define AES_SECURITY 192

//#define MR_PAIRING_BLS    // AES-256 security
//#define AES_SECURITY 256
//*********************************************

//pairing type 3 library
#include "pairing_3.h"
//library for hash length
#define HASH_LEN 20
/*
This hash function takes an input character string and order of Big data type and outputs
a hash of Big data type. 
*/
Big H(char *string,Big n)
{ // Hash a zero-terminated string to a number h < modulus n
  // Then increment this number until  (h/n)=+1
    Big h;
    char s[HASH_LEN];
    int i,j; 
    sha sh;

    shs_init(&sh);

    for (i=0;;i++)
    {
        if (string[i]==0) break;
        shs_process(&sh,string[i]);
    }
    shs_hash(&sh,s);

    h=1; j=0; i=1;
    forever
    {
        h*=256; 
        if (j==HASH_LEN)  {h+=i++; j=0;}
        else         h+=s[j++];
        if (h>=n) break;
    }
    while (jacobi(h,n)!=1) h+=1;
    h%=n;
    return h;
}
int main()
{
	// initialise pairing-friendly curve
	PFC pfc(AES_SECURITY); 
	// get handle on mip (Miracl Instance Pointer)
	miracl *mip=get_mip();  
	/*Initialize  time seed. The time_t datatype is a data type in the ISO C library defined for storing
	system time values. Such values are returned from the standard time() library function.
	*/
	time_t seed;
	cout << "********************* Predicate Encryption ******************* " << endl;
	//////////////
	/// Setup
	/////////////
	cout << "Starting Setup" << endl;

	//get pairing-friendly group order
	Big det,order=pfc.order();  // 
	/*mip is the Miracl Instance Pointer. mip->IOBASE=256 simply changes the base to 256.
	We take input in base 256 to componsate all the real world letters and special characters.
	Which are easy to represented in base 256 system of numbers.
	*/
	mip->IOBASE=256;

	// Get the id and hash it using function "H".
	Big id= H("ALICE",order);
	
		
	/*mip is the Miracl Instance Pointer. mip->IOBASE=16 simply changes the base to 16(Hexadecimal).
	We use the hexadecimal numbers to make coding for microprocessor. But it converts that to binary
	for computation. After the computation the result will be in hexadecimal format by inverse conversion.
	*/
	mip->IOBASE=16;
	//initialize elements of type big (Setup parameters).
	Big alpha,beta,gama,a1,a2,z,r,R,s,mm;
	//Vector V. The reason for selecting the vectors v and x has been explained in above Note.
	Big v[2]={1,id};
	//Vector x.
	Big x[2]={id,-1};
	//initialize Group G2 elements 
	G2 h,h1,h21,h22,k1,k2,k3,k0,k01,k02,k03;
	//initialize Group G1 elements 
	G1 g,g1,g21,g22, g0,c1,c21,c22,c3;
	//initialize Group GT elements 
	GT m1,Y,M,c0,m2,m3;
	// initialise random numbers
	time(&seed);
	 //creat random number of long type from the seed value.
    irand((long)seed);
	//randomly generate gama
	pfc.random(gama);
	//randomly generate alpha
	pfc.random(alpha);
	//randomly generate beta
	pfc.random(beta);
	//randomly generate a1
	pfc.random(a1);
	//randomly generate a2
	pfc.random(a2);
	//randomly generate z
	pfc.random(z);
	//randomly generate g
	pfc.random(g);
	//precompute g
	pfc.precomp_for_mult(g);
	//randomly generate h
	pfc.random(h);
	//precompute g
	pfc.precomp_for_mult(h);
	/*
	gama*g instead of g^gama, We changed multiplicative group operations
	into additive group operations, becasue there is no function to calculate G1^Big.
	Here G1 and Big are elements of data type G1 and Big respectively. Similarly we 
	changed multiplicative group operations into additive group operations for g21, g22,
	g0, h1, h21 and h22*/
	g1=pfc.mult(g,gama);
	//g*a1 instead of g^a1.
	g21=pfc.mult(g,a1);
	//g*a2 instead of g^a2.
	g22=pfc.mult(g,a2);
	//g*z instead of g^z.
	g0=pfc.mult(g,z);
	//h*gama instead of h^gama.
	h1=pfc.mult(h,gama);
	//h*a1 instead of h^a1.
	h21=pfc.mult(h,a1);
	//h*a2 instead of h^a2.
	h22=pfc.mult(h,a2);
	//precompute g0
	pfc.precomp_for_mult(g0);
	//precompute g1
	pfc.precomp_for_mult(g1);
	//precompute g21
	pfc.precomp_for_mult(g21);
	//precompute g22
	pfc.precomp_for_mult(g22);	
	//precompute h1
	pfc.precomp_for_mult(h1);
	//precompute h21
	pfc.precomp_for_mult(h21);
	//precompute h22
	pfc.precomp_for_mult(h22);
	/*Y=e(g,h*(modmult(alpha,beta,order)). We used modmult() to keep the product of alpha and beta 
	with in the order. Otherwise if we direclty calculate alpha*beta, it will be a very large number.
	Then the compiler will give the error "number too large". In order to avoid this error, we use modmult().
	We calculated h*(modmult(alpha,beta,order)) instead of h^modmult(alpha,beta,order). Here  We changed
	multiplicative group operations into additive group operations, becasue there is no function to calculate
	G1^Big. Here G1 and Big are elements of data type G1 and Big respectively. Then we calculated
	e(g,h*(modmult(alpha,beta,order)).

	*/
	Y=pfc.pairing(pfc.mult(h,modmult(alpha,beta,order)),g); 
	//precomute Y.
	pfc.precomp_for_power(Y);
	
	cout << "Setup Completed" << endl;
	cout << "**************************************************************" << endl;

	cout << "Starting KeyGen" << endl;
	//////////////////////////////
	// kEY Gen
	/////////////////////////////
	////randomly generate r
	pfc.random(r);
	//randomly generate R
	pfc.random(R);
	/*.
	k0=	h*(modmult(alpha,beta,order))
	k01=h21*(modmult(v[0],r,order))
	k02=h22*(modmult(v[1],r,order))
	k03=pfc.mult(h1,R)
	k1=k0+k01+k02+k03
	We can also write k1=(h*(odmult(alpha,beta,order)))+(h2i*(odmult(Vj,r,order)))+(h1*R).
	Here i=1,2 and j=0,1. We used modmult() to keep the product Big*Big (Big is an element of Big data type) 
	with in the order. Otherwise if we direclty calculate Big*Big, it will be a very large number.
	Then the compiler will give the error "number too large". In order to avoid this error, we use modmult().
	We calculate G1*Big instead of G1^Big and G2+G2 insted of G2*G2. Multiplicative group operations are 
	changed into additive group operations, becasue there is no function to calculate G1^Big and G2*G2.
	Here G1, G2 and Big are elements of data type G1, G2 and Big respectively. 
	*/
	k0=	pfc.mult(h,(modmult(alpha,beta,order)));
	k01=pfc.mult(h21,(modmult(v[0],r,order)));
	k02=pfc.mult(h22,(modmult(v[1],r,order)));
	k03=(pfc.mult(h1,R));
	k1=k0+k01+k02+k03;

	/*We calculate h*r instead of h^r. Multiplicative group operations are changed into additive group operations,
	becasue there is no function to calculate G2^Big. Here G2 and Big are elements of data type G2 and Big respectively.
	*/
	k2=pfc.mult(h,r);

	/*We calculate h*R instead of h^R. Multiplicative group operations are changed into additive group operations,
	becasue there is no function to calculate G2^Big. Here G2 and Big are elements of data type G2 and Big respectively.
	*/
	k3=pfc.mult(h,R);
	cout << "KeyGen Completed" << endl;
	cout << "**************************************************************" << endl;

	//////////////////
	//Encryption
	/////////////////

	cout << "Starting Encryption" << endl;
	//declare ran of Big data type
	Big ran;
	//randomly generate ran
	pfc.random(ran);
	/*M=Y^ran is session key. We use M as session key. Because M is of GT data type and there is no
	function to convert GT type data in to plain text.
	*/
	M=pfc.power(Y,ran); 
	/*mip is the Miracl Instance Pointer. mip->IOBASE=256 simply changes the base to 256.
	We take input in base 256 to componsate all the real world letters and special characters.
	Which are easy to represented in base 256 system of numbers.
	*/
	mip->IOBASE=256;
	/*Plaintext is actually the message to encrypted. Take char data type input and convrt it to Big data type.
	We took Plaintext of Big data type because it is convenient to print Big data type variables.
	*/
	Big Plaintext=(char *)"Test message"; 
	// Print Plaintext.
	cout << "Message to be Encrypted= " << Plaintext << endl;  
	/*mip is the Miracl Instance Pointer. mip->IOBASE=16 simply changes the base to 16(Hexadecimal).
	We use the hexadecimal numbers to make coding for microprocessor. But it converts that to binary
	for computation. After the computation the result will be in hexadecimal format by inverse conversion.
	*/
	mip->IOBASE=16;
	//randomly generate s
	pfc.random(s);
	// c0=M*(Y^s).
	c0=M*pfc.power(Y,s);
	/*
	We calculate g*s instead of g^s. Multiplicative group operations are changed into additive group operations,
	becasue there is no function to calculate G1^Big. Here G1 and Big are elements of data type G1 and Big respectively.
	*/
	c1=pfc.mult(g,s);
	/*
	c2i=(g0*modmult(xj,s,order))+(g2i^s).Here i=1,2 and j=0,1. We used modmult() to keep the product Big*Big (Big is and
	element of Big data type) with in the order. Otherwise if we direclty calculate Big*Big, it will be a very large number.
	Then the compiler will give the error "number too large". In order to avoid this error, we use modmult(). We calculate
	G1*Big instead of G1^Big and G2+G2 insted of G2*G2. Multiplicative group operations are changed into additive group
	operations, becasue there is no function to calculate G1^Big and G2*G2. Here G1, G2 and Big are elements of data type
	G1, G2 and Big respectively.
	*/
	c21=pfc.mult(g0,(modmult(x[0],s,order)))+(pfc.mult(g21,s));
	c22=pfc.mult(g0,(modmult(x[1],s,order)))+(pfc.mult(g22,s));
	/*
	c3=g1*s. We calculate g1*s instead of g1^s. Multiplicative group operations are changed into additive group operations,
	becasue there is no function to calculate G1^Big. Here G1 and Big are elements of data type G1 and Big respectively.
	*/
	c3=pfc.mult(g1,s);
	// CT =Plaintext xor H(M). We hashed M to convert it in to Big data type. Because xor takes input of Big data type. 
	Big CT= lxor(Plaintext,pfc.hash_to_aes_key(M));

	cout << "Encryption Completed" << endl;

	cout << "**************************************************************" << endl;
	//////////////////
	//Decryption
	/////////////////
	cout << "Starting Decryption" << endl;

	/*
	m1=e((c21*v0)+(c22*v1),k2)
	m2=m1*e(k3,c3)/e(k1,c1)
	m3=c0*m2
	We can also write m3=((e((c21*v0)+(c22*v1),k2)*e(k3,c3))/e(k1,c1))*c0. We calculate C2i*Vj instead of  C2i^Vj.Here i=1,2
	and j=0,1. Multiplicative group operations are changed into additive group operations,becasue there is no function to
	calculate G1^Big. Here G1 and Big are elements of data type G1 and Big respectively. Here m1 is decrypted session key
	*/
	m1=pfc.pairing(k2,(pfc.mult(c21,v[0])+pfc.mult(c22,v[1])));
	m2=(m1*pfc.pairing(k3,c3))/pfc.pairing(k1,c1);
	m3=c0*m2;
	
	
	/*mip is the Miracl Instance Pointer mip->IOBASE=256. It simply changes the base to work with, we will converted
	the result back to base 256 from base 16, because if we output the result in base 16 it will not be same as 
	the ipnut. Becuase the input was in base 256. So we need to change back the base of number system to 256. So that
	we can get the same output and input display.
	*/
	mip->IOBASE=256;
	/* Here we output the decrypted value of plaintext. For this purpose we xor CT with hash of m1. Because in encryption
	we xor Plaintext with hash of M. 
	*/
	cout << "Decrypted Message= " << lxor(CT,pfc.hash_to_aes_key(m3)) << endl;
	cout << "Decryption Completed" << endl;
	cout << "**************************************************************" << endl;
	system("pause");
	return 0;
}
