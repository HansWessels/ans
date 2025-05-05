/*
**
** ANS coder test
**
*/

/* 
** gcc ../q/c_code/ans/ans.c -lm -lgmp -o ans
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

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
{ /* init table using a Linear congruential generator, x = (ax+c) mod m, m = 2^32, a = 1103515245, c = 12345 */
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
{ /* init table by just inserting symbols, smallest symbol first */
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
{ /* init table by just inserting symbols, most frequent symbol first */
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

void make_precise_lff_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* according to "Jarek Duda" paper, least frequent first  */
	float* precise_pos;
	int i;
	int pos;
	precise_pos = malloc(symbol_size*sizeof(float));
	if(precise_pos == NULL)
	{
		printf("Malloc error for %s!\n", __func__);
		exit(-1);
	}
	for(i=0; i<symbol_size; i++)
	{
		float pp;
		if(symbols_count[i]!=0)
		{
			pp=(float)table_size/(float)(2*symbols_count[i]);
		}
		else
		{
			pp=(float)2*table_size;
		}
		precise_pos[i]=pp;
	}
	for(pos=0; pos<table_size; pos++)
	{
		float smallest_pos=(float)table_size;
		int smallest_count=symbol_size;
		int smallest_symbol=0;
		int i;
		for(i=0; i<symbol_size; i++)
		{
			float pp;
			pp=precise_pos[i];
			if(pp<smallest_pos)
			{
				smallest_pos=pp;
				smallest_count=symbols_count[i];
				smallest_symbol=i;
			}
			else if((pp==smallest_pos) && (symbols_count[i]<smallest_count))
			{
				smallest_count=symbols_count[i];
				smallest_symbol=i;
			}
		}
		ans_table[pos]=smallest_symbol;
		precise_pos[smallest_symbol]+=(float)table_size/(float)symbols_count[smallest_symbol];
	}
	free(precise_pos);
}

/*
**
** Priem table:
Table(4): Rel errsq, priem = 3, err = 0.000000
Table(8): Rel errsq, priem = 3, err = 1.000000
Table(16): Rel errsq, priem = 3, err = 2.461250
Table(32): Rel errsq, priem = 23, err = 3.046979
Table(64): Rel errsq, priem = 23, err = 16.294859
Table(128): Rel errsq, priem = 47, err = 50.583920
Table(256): Rel errsq, priem = 181, err = 208.831421
Table(512): Rel errsq, priem = 881, err = 729.684570
Table(1024): Rel errsq, priem = 8053, err = 3046.951416
Table(2048): Rel errsq, priem = 883, err = 12254.189453
Table(4096): Rel errsq, priem = 1607, err = 46985.398438
Table(8192): Rel errsq, priem = 18149, err = 189977.968750
Table(16384): Rel errsq, priem = 11807, err = 762292.062500
Table(32768): Rel errsq, priem = 1399, err = 2796467.750000
Table(65536): Rel errsq, priem = 29333, err = 8780226.000000
**
*/
void make_priem_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* init table by using a prime to spave out the data */
	int i;
	int table_pos=0;
	int const priem=1607;
	for(i=0; i<symbol_size; i++)
	{
		int n;
		n=symbols_count[i];
		while(n>0)
		{
			ans_table[table_pos]=i;
			table_pos+=priem;
			table_pos%=table_size;
			n--;
		}
	}
}

void make_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
//	make_rng_table(ans_table, symbols_count, symbol_size, table_size);
//	make_simple_table(ans_table, symbols_count, symbol_size, table_size);
//	make_sorted_simple_table(ans_table, symbols_count, symbol_size, table_size);
	make_precise_lff_table(ans_table, symbols_count, symbol_size, table_size);
//	make_priem_table(ans_table, symbols_count, symbol_size, table_size);
	{ /* Sanety check */
		int count[4096]={0};
		int i;
		for(i=0; i<table_size; i++)
		{
			count[ans_table[i]]++;
		}
		for(i=0; i<symbol_size; i++)
		{
			if(symbols_count[i]!=count[i])
			{
				printf("Error! %02X : count = %i, tab_count=%i\n", i, symbols_count[i], count[i]);
				exit(-1);
			}
		}
	}
}

void encode_symbol(int symbol, mpz_t x, int ans_table[], int symbols_count[], int table_size)
{
	int r;
	int pos;
	mpz_add_ui(x, x, 1);
	r=(int)mpz_tdiv_q_ui(x, x, symbols_count[symbol]);
	mpz_mul_ui(x, x, table_size);
	pos=0;
	for(;;)
	{ /* zoek r'th occurence of symbol in table */
		if(ans_table[pos]==symbol)
		{
			if(r==0)
			{
				break;
			}
			else
			{
				r--;
			}
		}
		pos++;
	}
	mpz_add_ui(x, x, pos);
}

void encode_ans_uint8_t(mpz_t x, uint8_t *data, unsigned long data_size, int ans_table[], int symbols_count[], int table_size)
{
	while(data_size!=0)
	{
		data_size--;
		encode_symbol(data[data_size], x, ans_table, symbols_count, table_size);
	}
}

int decode_symbol(mpz_t x, int ans_table[], int symbols_count[], int table_size)
{
	int symbol;
	int index;
	int pos;
	index=(int)mpz_tdiv_q_ui(x, x, table_size);
	symbol=ans_table[index];
	pos=0;
	while(--index>=0)
	{
		if(ans_table[index]==symbol)
		{
			pos++;
		}
	}
	mpz_mul_ui(x, x, symbols_count[symbol]);
	pos--;
	if(pos>0)
	{
		mpz_add_ui(x, x, pos);
	}
	else if(pos<0)
	{
		mpz_sub_ui(x, x, 1);
	}
	return symbol;
}

void decode_ans_ascii_t(mpz_t x, int ans_table[], int symbols_count[], int table_size)
{
	while(mpz_sgn(x)!=0)
	{
		int symbol;
		symbol=decode_symbol(x, ans_table, symbols_count, table_size);
//		printf("%c", symbol);
	}
}

#define TABLE_SIZE 4096
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
		printf("Make table\n");
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
					#if 0
					if((i>32) && (i<127))
					{
						printf("freq[%4c] = %3li = %4i = %f bits\n", i, freq[i], symbols_count[i], bits);
					}
					else
					{
						printf("freq[0x%02X] = %3li = %4i = %f bits\n", i, freq[i], symbols_count[i], bits);
					}
					#endif
				}
			}
			printf("Total_bits %li = %f\n", size*8, total_bits);
			#if 0
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
			#endif
			{ /* encodeding and decoding */
				mpz_t x;
				mpz_init(x);
				unsigned long bits;
//				printf("Result:\n");
				encode_ans_uint8_t(x, data, size, table, symbols_count, TABLE_SIZE);
//				nibbles=mpz_out_str (stdout, 16, x);
				bits=mpz_sizeinbase (x, 2);
				printf("Compressed = %li bits\n", bits);
//				decode_ans_ascii_t(x, table, symbols_count, TABLE_SIZE);
				mpz_clear(x);
			}
		}
		free(data);
		i++;
	}
	return 0;
}
