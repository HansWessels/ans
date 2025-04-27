/*
**
** Test om de beste priemgetallen te vinden voor het vullen van een ans-tabel
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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
	int b=(j-i)& mod;
	if(a<b)
	{
		return a;
	}
	return b;
}

int test(int table_size, int primes[], int max_prime)
{
	int table[1<<MAX_TABLE_SIZE];
	int min_afstand[1<<MAX_TABLE_SIZE];
	int start=0;
	uint64_t min_fail=(uint64_t)-1;
	int best_priem=0;
	for(start=0;start<max_prime; start++)
	{
		int i;
		int pos=0;
		int mod=table_size-1;
		int priem=primes[start];
		int min_a;
		int fail_count=0;
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
		for(i=0; i<table_size/4; i++)
		{
			min_afstand[i]=table_size;
		}
		for(i=1; i<table_size/4; i++)
		{
			int j;
			for(j=0; j<table_size; j++)
			{
				int k;
				int pos;
				int start;
				min_a=table_size;
				for(k=0; k<table_size; k++)
				{
					table[k]=0;
				}
				pos=j;
				for(k=0; k<i; k++)
				{
					table[pos]=1;
					pos+=priem;
					pos&=mod;
				}
				k=0;
				while(table[k]==0)
				{
					k++;
				}
				start=k;
				k++;
				while(k<table_size)
				{
					if(table[k]!=0)
					{
						int a;
						a=afstand(start, k, mod);
						if(a<(table_size/(i)))
						{
							fail_count++;	
						}
					}
					k++;
				}
			}
		}
		if(fail_count<min_fail)
		{
			min_fail=fail_count;
			best_priem=priem;
		}
	}
	printf("Table(%i), priem = %i, fail = %li\n", table_size, best_priem, min_fail);
}

int main(void)
{
	int primes[1<<MAX_TABLE_SIZE];
	int count;
	int i;
	count=make_primes(primes);
	for(i=8; i<=MAX_TABLE_SIZE; i++)
	{
		test(1<<i, primes, count);
	}
	return 0;
}
	