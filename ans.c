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

#define MAX_SYMBOL_COUNT 512
#define MAX_ANS_TABLE 65536

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
	double precise_pos[512];
	int i;
	int pos;
	for(i=0; i<symbol_size; i++)
	{
		double pp;
		if(symbols_count[i]!=0)
		{
			pp=(double)table_size/(double)(2*symbols_count[i]);
		}
		else
		{
			pp=(double)(2*table_size);
		}
		precise_pos[i]=pp;
	}
	for(pos=0; pos<table_size; pos++)
	{
		double smallest_pos=(double)(2*table_size);
		int smallest_count=symbol_size;
		int smallest_symbol=0;
		int i;
		for(i=0; i<symbol_size; i++)
		{
			double pp;
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
		precise_pos[smallest_symbol]+=(double)table_size/(double)symbols_count[smallest_symbol];
	}
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
	int const priem=47;
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

int add_symbol(int table[], int size, int symbol, int count)
{ /* Bresenham's line algorithm */
	int d;
	int i;
	int j=size-1;
	size+=count;
	d=2*count-size;
	for(i=size; i>0;)
	{
		i--;
		if(d>0)
		{
			table[i]=symbol;
			d-=2*size;
		}
		else
		{
			table[i]=table[j];
			j--;
		}
		d+=2*count;
	}
	return size;
}

int add_symbols(int table[], int size, int symbols[], int symbols_count, int count)
{ /* Bresenham's line algorithm, using symbol sequence */
	int d;
	int i;
	int pos=0;
	int j=size-1;
	size+=count;
	d=2*count-size;
	for(i=size; i>0;)
	{
		i--;
		if(d>0)
		{
			table[i]=symbols[pos++];
			if(pos>=symbols_count)
			{
				pos=0;
			}
			d-=2*size;
		}
		else
		{
			table[i]=table[j];
			j--;
		}
		d+=2*count;
	}
	return size;
}

void make_us_bh_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* init table by using Bresenham's line algorithm */
	int size=0;
	int i;
	for(i=0; i<symbol_size; i++)
	{
		if(symbols_count[i]>0)
		{
			size=add_symbol(ans_table, size, i, symbols_count[i]);
		}
	}
}

void make_sb_bh_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* init table by using Bresenham's line algorithm, biggest first */
	int size=0;
	int i;
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
			if(symbols_count[i]>sub_max)
			{
				if(symbols_count[i]==max)
				{
					size=add_symbol(ans_table, size, i, symbols_count[i]);
				}
				else if(symbols_count[i]<max)
				{
					sub_max=symbols_count[i];
				}
			}
		}
		max=sub_max;
	}
}

void make_sb_bl_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* init table by using Bresenham's line algorithm, biggest first, alle symbolen met zelfde frequentie tegelijk toevoegen */
	int size=0;
	int i;
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
		int line[MAX_SYMBOL_COUNT];
		int pos=0;
		for(i=0; i<symbol_size; i++)
		{
			if(symbols_count[i]>sub_max)
			{
				if(symbols_count[i]==max)
				{
					line[pos++]=i;
				}
				else if(symbols_count[i]<max)
				{
					sub_max=symbols_count[i];
				}
			}
		}
		size=add_symbols(ans_table, size, line, pos, pos*max);
		max=sub_max;
	}
}

void make_ss_bh_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* init table by using Bresenham's line algorithm, smallest first */
	int size=0;
	int i;
	int min=1;
	while(min<table_size)
	{
		int sub_min=table_size;
		for(i=0; i<symbol_size; i++)
		{
			if(symbols_count[i]<sub_min)
			{
				if(symbols_count[i]==min)
				{
					size=add_symbol(ans_table, size, i, symbols_count[i]);
				}
				else if(symbols_count[i]>min)
				{
					sub_min=symbols_count[i];
				}
			}
		}
		min=sub_min;
	}
}

void make_ss_bl_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{ /* init table by using Bresenham's line algorithm, smallest first, alle symbolen met zelfde frequentie tegelijk toevoegen */
	int size=0;
	int i;
	int min=1;
	while(min<table_size)
	{
		int sub_min=table_size;
		int pos=0;
		int line[MAX_SYMBOL_COUNT];
		for(i=0; i<symbol_size; i++)
		{
			if(symbols_count[i]<sub_min)
			{
				if(symbols_count[i]==min)
				{
					line[pos++]=i;
				}
				else if(symbols_count[i]>min)
				{
					sub_min=symbols_count[i];
				}
			}
		}
		if(pos>0)
		{
			size=add_symbols(ans_table, size, line, pos, pos*min);
		}
		min=sub_min;
	}
}

void make_hdl_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
	int used[MAX_SYMBOL_COUNT]={0};
	int start[MAX_SYMBOL_COUNT];
	int pos;
	{ /* zoek meest frequente symbool */
		int i;
		int max=0;
		int max_symbol=0;
		for(i=0; i<symbol_size; i++)
		{
			if(symbols_count[i]>max)
			{
				max=symbols_count[i];
				max_symbol=i;
			}
		}
		ans_table[0]=max_symbol;
		used[max_symbol]++;
		start[max_symbol]=0;
	}
	for(pos=1; pos<table_size; pos++)
	{
		int i;
		float max=0.0;
		int max_symbol=0;
		for(i=0; i<symbol_size; i++)
		{
			if(used[i]!=0)
			{
				float delta=(float)symbols_count[i]/(float)table_size-(float)used[i]/(float)(pos-start[i]);
				if(delta>=max)
				{
					if(used[i]<symbols_count[i])
					{
						if(delta==max)
						{ /* bij gelijk spel heeft minst frequente symbool voorrang */
							if(symbols_count[i]<symbols_count[max_symbol])
							{
								max_symbol=i;
							}
						}
						else
						{
							max_symbol=i;
							max=delta;
						}
					}
				}
			}
		}
		if(max==0.0)
		{
			int max=0;
			for(i=0; i<symbol_size; i++)
			{
				if(used[i]==0)
				{
					if(symbols_count[i]>max)
					{
						max_symbol=i;
						max=symbols_count[i];
					}
				}
			}
			if(max==0)
			{
				float max=0.0;
				for(i=0; i<symbol_size; i++)
				{
					if(used[i]!=0)
					{
						float delta=(float)used[i]/(float)(pos-start[i])-(float)symbols_count[i]/(float)table_size;
						if(delta>=max)
						{
							if(used[i]<symbols_count[i])
							{
								if(delta==max)
								{ /* bij gelijk spel heeft minst frequente symbool voorrang */
									if(symbols_count[i]<symbols_count[max_symbol])
									{
										max_symbol=i;
									}
								}
								else
								{
									max_symbol=i;
									max=delta;
								}
							}
						}
					}
				}
			}
			else
			{
				start[max_symbol]=pos;
			}
		}
		ans_table[pos]=max_symbol;
		used[max_symbol]++;
	}	
}

void make_table(int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
//	make_rng_table(ans_table, symbols_count, symbol_size, table_size);
//	make_simple_table(ans_table, symbols_count, symbol_size, table_size);
//	make_sorted_simple_table(ans_table, symbols_count, symbol_size, table_size);
//	make_precise_lff_table(ans_table, symbols_count, symbol_size, table_size);
//	make_priem_table(ans_table, symbols_count, symbol_size, table_size);
//	make_us_bh_table(ans_table, symbols_count, symbol_size, table_size);
//	make_sb_bh_table(ans_table, symbols_count, symbol_size, table_size);
//	make_ss_bh_table(ans_table, symbols_count, symbol_size, table_size);
//	make_sb_bl_table(ans_table, symbols_count, symbol_size, table_size);
//	make_ss_bl_table(ans_table, symbols_count, symbol_size, table_size);
	make_hdl_table(ans_table, symbols_count, symbol_size, table_size);
	{ /* Sanety check */
		int count[MAX_ANS_TABLE]={0};
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
	if(1)
	{ /* print table */
		int i;
		for(i=0; i<symbol_size; i++)
		{
			if(symbols_count[i]>0)
			{
				int j;
				printf("%02X(%4i):", i, symbols_count[i]);
				for(j=0; j<table_size; j++)
				{
					if(ans_table[j]==i)
					{
						printf("@");
					}
					else
					{
						printf("_");
					}
				}
				printf("\n");
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
