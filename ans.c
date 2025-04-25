/*
**
** ANS coder test
**
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int64_t load_file(char* infile, uint8_t** data_in)
{
	FILE* f;
	uint8_t* data;
	int64_t size;
	*data_in=NULL;
	f = fopen(infile, "rb");
	if (f == NULL)
	{
		fprintf(stderr, "File open error %s!\n", infile);
		return -1;
	}
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	data = (uint8_t*)malloc(size);
	if (data == NULL)
	{
		fprintf(stderr, "Malloc error voor file data %s !\n", infile);
		fclose(f);
		return -1;
	}
	fseek(f, 0, SEEK_SET);
	if(fread(data, 1, size, f)!=size)
	{
		fprintf(stderr, "Read error %s\n", infile);
		fclose(f);
		return -1;
	}
	fclose(f);
	*data_in=data;
	return size;
}

void freq_count(uint8_t *data, int64_t size, int64_t freq[], int symbol_size)
{
	bzero(freq, symbol_size*sizeof(int64_t));
	while(size>0)
	{
		size--;
		freq[data[size]]++;
	}
}

void verdeel_symbols(int64_t freq[], int symbols_count[], int symbol_size, int table_size)
{
	int i;
	bzero(symbols_count, symbol_size*sizeof(int));
	for(i=0; i<symbol_size; i++)
	{
		if(freq[i]>0)
		{
			symbols_count[i]=1;
			table_size--;
		}
	}
	while(table_size!=0)
	{
		int i;
		float worst=0;
		int worst_i=0;
		for(i=0; i<symbol_size; i++)
		if(freq[i]>0)
		{
			float cost;
			cost=(float)freq[i]/(float)symbols_count[i];
			if(cost>worst)
			{
				worst=cost;
				worst_i=i;
			}
		}
		symbols_count[worst_i]++;
		table_size--;
	}
}

void make_rng_table(int table[], int symbols_count[], int symbol_size, int table_size)
{
	uint32_t seed=0;
	int symbol;
	int i;
	for(i=0; i<table_size; i++)
	{
		table[i]=-1;
	}
	for(symbol=0; symbol<symbol_size; symbol++)
	{
		int n;
		n=symbols_count[symbol];
		while(n>0)
		{
			int pos;
			int k=0;
			seed=1103515245UL*seed+12345UL;
			pos=(table_size*((seed>>15)&0xffff))>>16; /* calculate insertion pos */
			table_size--;
			for(;;)
			{
				if(table[k]>=0) /* find insertion pos */
				{
					k++;
				}
				else if(pos>0)
				{
					pos--;
					k++;
				}
				else
				{
					break;
				}
			}
			table[k]=symbol; /* insert symbol in ans_table */
			n--;
		}
	}
}

void make_simple_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
	int i;
	int table_pos=0;
	for(i=0; i<symbol_size; i++)
	{
		int n;
		n=symbols_count[i];
		while(n>0)
		{
			ans_table[table_pos]=i;
			table_pos++;
			n--;
		}
	}
}

void make_sorted_simple_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
	int i;
	int table_pos=0;
	int max=0;
	for(i=0; i<symbol_size; i++)
	{
		int n;
		n=symbols_count[i];
		if(n>max)
		{
			max=n;
		}
	}
	while(max>0)
	{
		int sub_max=0;	
		for(i=0; i<symbol_size; i++)
		{
			int n;
			n=symbols_count[i];
			if(n>sub_max)
			{
				if(n==max)
				{
					while(n>0)
					{
						ans_table[table_pos]=i;
						table_pos++;
						n--;
					}
				}
				else if(n<max)
				{
					sub_max=n;
				}
			}
		}
		max=sub_max;
	}
}

void make_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
//	make_rng_table(ans_table, symbols_count, symbol_size, table_size);
//	make_simple_table(ans_table, symbols_count, symbol_size, table_size);
	make_sorted_simple_table(ans_table, symbols_count, symbol_size, table_size);
}

#define TABLE_SIZE 128
#define SYMBOL_SIZE 256

int main(int argc, char* argv[])
{
	int i;
	i=1;
	while(i<argc)
	{
		int64_t freq[SYMBOL_SIZE];
		int symbols_count[SYMBOL_SIZE];
		int table[TABLE_SIZE];
		int64_t size;
		uint8_t* data;
		printf("Loding file %s\n",argv[i]);
		size=load_file(argv[i], &data);
		if(data==NULL)
		{
			printf("Error loading %s\n", argv[i]);
			return -1;
		}
		printf("File size = %li\n", size);
		freq_count(data, size, freq, SYMBOL_SIZE);
		verdeel_symbols(freq, symbols_count, SYMBOL_SIZE, TABLE_SIZE);
		make_table(table, symbols_count, SYMBOL_SIZE, TABLE_SIZE);
		{
			int i;
			float total_bits=0;
			for(i=0; i<256; i++)
			{
				if(freq[i]!=0)
				{
					float bits;
					bits=log2f((float)TABLE_SIZE/(float)symbols_count[i]);
					total_bits+=bits*freq[i];
					if((i>32) && (i<127))
					{
						printf("freq[%4c] = %3li = %4i = %f bits\n", i, freq[i], symbols_count[i], bits);
					}
					else
					{
						printf("freq[0x%02X] = %3li = %4i = %f bits\n", i, freq[i], symbols_count[i], bits);
					}
				}
			}
			printf("Total_bits %li = %f\n", size*8, total_bits);
			for(i=0; i<TABLE_SIZE; i++)
			{
				if((table[i]>32) && (table[i]<127))
				{
					printf("%c ", table[i]);
				}
				else
				{
					printf("0x%02X ", table[i]);
				}
			}
			printf("\n");
		}
		free(data);
		i++;
	}
	return 0;
}
