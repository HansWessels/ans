/*
**
** Test om de beste priemgetallen te vinden voor het vullen van een ans-tabel
**
*/


/*
resultaten:
Table(4)
Abs err, priem = 3, err = 0
Rel err, priem = 3, err = 0.000000
Rel errsq, priem = 3, err = 0.000000
Table(8)
Abs err, priem = 3, err = 8
Rel err, priem = 3, err = 1.000000
Rel errsq, priem = 3, err = 1.000000
Table(16)
Abs err, priem = 3, err = 33
Rel err, priem = 3, err = 3.650000
Rel errsq, priem = 3, err = 2.461250
Table(32)
Abs err, priem = 23, err = 78
Rel err, priem = 23, err = 5.741666
Rel errsq, priem = 23, err = 3.046979
Table(64)
Abs err, priem = 19, err = 304
Rel err, priem = 5, err = 28.490276
Rel errsq, priem = 23, err = 16.294859
Table(128)
Abs err, priem = 47, err = 1063
Rel err, priem = 331, err = 105.414497
Rel errsq, priem = 47, err = 50.583920
Table(256)
Abs err, priem = 181, err = 4010
Rel err, priem = 1657, err = 397.973755
Rel errsq, priem = 181, err = 208.831421
Table(512)
Abs err, priem = 881, err = 13852
Rel err, priem = 881, err = 1456.786011
Rel errsq, priem = 881, err = 729.684570
Table(1024)
Abs err, priem = 139, err = 57108
Rel err, priem = 139, err = 5964.334961
Rel errsq, priem = 8053, err = 3046.951416
Table(2048)
Abs err, priem = 449, err = 227729
Rel err, priem = 5261, err = 23755.970703
Rel errsq, priem = 883, err = 12254.189453
Table(4096)
Abs err, priem = 1607, err = 859253
Rel err, priem = 1607, err = 93559.109375
Rel errsq, priem = 1607, err = 46985.398438
Table(8192)
Abs err, priem = 6427, err = 3485166
Rel err, priem = 6427, err = 375175.281250
Rel errsq, priem = 18149, err = 189977.968750
Table(16384)
Abs err, priem = 11807, err = 13934098
Rel err, priem = 3467, err = 1531829.250000
Rel errsq, priem = 11807, err = 762292.062500
Table(32768)
Abs err, priem = 1399, err = 54638264
Rel err, priem = 5323, err = 5984361.500000
Rel errsq, priem = 1399, err = 2796467.750000
Table(65536)
Abs err, priem = 63601, err = 219498354
Rel err, priem = 28463, err = 16251780.000000
Rel errsq, priem = 29333, err = 8780226.000000

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define MAX_TABLE_SIZE 16

int make_primes(int primes[])
{
	int count=0;
	int i;
	for(i=3; i<(1<<MAX_TABLE_SIZE); i++)
	{
		int j;
		int is_priem=1;
		for(j=2; j<i; j++)
		{
			if((i%j)==0)
			{
				is_priem=0;
				break;
			}
		}
		if(is_priem!=0)
		{
			primes[count]=i;
			count++;
		}
	}
	return count;
}

int afstand(int i, int j, int mod)
{
	int a=(i-j)&mod;
	int b=(j-i)&mod;
	if(a<b)
	{
		return a;
	}
	return b;
}

int test(int table_size, int primes[], int max_prime)
{
	int table[1<<MAX_TABLE_SIZE];
	int start=0;
	uint64_t min_abs_err=(uint64_t)-1;
	uint64_t min_abs_err2=(uint64_t)-1;
	float min_rel_err=1e18;
	float min_rel_err2=1e18;
	float min_rel_errsq=1e18;
	float min_rel_errsq2=1e18;
	int best_abs_priem=0;
	int best_abs2_priem=0;
	int best_rel_priem=0;
	int best_rel2_priem=0;
	int best_relsq_priem=0;
	int best_relsq2_priem=0;
	for(start=0;start<max_prime; start++)
	{
		int i;
		int pos=0;
		int mod=table_size-1;
		int mod2=mod>>1;
		int priem=primes[start];
		int min_a;
		float rel_err=0;
		float rel_errsq=0;
		int abs_err=0;
		int max_i;
		max_i=table_size/4;
		for(i=0; i<table_size; i++)
		{
			table[i]=-1;
		}
		for(i=0; i<table_size; i++)
		{
			if(table[pos]>=0)
			{
				printf("Error! Collission table %i met prime %i\n", table_size, priem);
				exit(-1);
			}
			table[pos]=i;
			pos+=priem;
			pos&=mod;
		}
		for(i=1; i<max_i; i++)
		{
			int pos;
			int j;
			int opt_a=table_size/i;
			int start=table_size-1;
			for(j=0; j<table_size; j++)
			{
				table[j]=0;
			}
			pos=start;
			for(j=0; j<i; j++)
			{
				table[pos]=1;
				pos+=priem;
				pos&=mod;
			}
			pos=0;
			while(pos<table_size)
			{
				if(table[pos]!=0)
				{
					int a;
					a=afstand(start, pos, mod2);
					if((a<opt_a) || (a>(opt_a+1)))
					{
						float rel=fabs((float)(opt_a-a)/(float)opt_a);
						abs_err+=abs(opt_a-a);
						rel_err+=rel;
						rel_errsq+=rel*rel;
					}
					start=pos;
				}
				pos++;
			}
		}
		if(abs_err<min_abs_err)
		{
			min_abs_err=abs_err;
			best_abs_priem=priem;
		}
		if(rel_err<min_rel_err)
		{
			min_rel_err=rel_err;
			best_rel_priem=priem;
		}
		if(rel_errsq<min_rel_errsq)
		{
			min_rel_errsq=rel_errsq;
			best_relsq_priem=priem;
		}
	}
	printf("Table(%i)\n", table_size);
	printf("Abs err, priem = %i, err = %li\n", best_abs_priem, min_abs_err);
	printf("Rel err, priem = %i, err = %f\n", best_rel_priem, min_rel_err);
	printf("Rel errsq, priem = %i, err = %f\n", best_relsq_priem, min_rel_errsq);
}

int main(void)
{
	int primes[1<<MAX_TABLE_SIZE];
	int count;
	int i;
	count=make_primes(primes);
	for(i=2; i<=MAX_TABLE_SIZE; i++)
	{
		test(1<<i, primes, count);
	}
	return 0;
}
	