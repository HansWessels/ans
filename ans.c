/*
**
** ANS coder test
**
*/

/* 
** gcc ../q/c_code/ans/ans.c -lm -lgmp -o ans
*/

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

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#define TEST_MODE
#define MAX_SYMBOL_COUNT 512
#define MAX_ANS_TABLE 65536
#define DELTA_STATE_OFFSET 0 //1607


#include <stdint.h>

#ifdef TEST_MODE
uint64_t global_bit_count=0;
#endif

/*
** ARJ CRC32 routines
**
** 2022 Hans Wessels
*/

/*
** this function makes a CRC32 table for ARJ CRC32
** call function with pointer to an unsigned long[256] array
*/

void make_crc32_table(uint32_t crc_table[])
{
  int max_bits=8;
  int bits=1;
  crc_table[0]=0;
  crc_table[1<<(max_bits-1)]=0xEDB88320UL;
  do
  {
    int i=(1<<bits)-1;
    do
    {
      int j=(1<<bits)-1;
      int done=i<<(max_bits-2*bits);
      uint32_t crc=crc_table[done<<bits];
      int offset=(int)((crc&j)<<(max_bits-bits));
      crc>>=bits;
      do
      {
        int tmp=j<<(max_bits-bits);
        crc_table[tmp+done]=crc^crc_table[tmp^offset];
      }
      while(j--);
    }
    while(--i);
    bits+=bits;
  }
  while(bits<max_bits);
}

/*
** this function calculates the ARJ CRC32 checksum over count bytes in data
** call function with count: nuber of bytes to CRC, data: pointer to data bytes
*/

uint32_t crc32(unsigned long count, uint8_t* data, uint32_t crc_table[])
{
	uint32_t crc=(uint32_t)-1;
	while(count--!=0)
	{
		uint8_t c = *data++;
		c ^= (uint8_t) crc;
		crc >>= 8;
		crc ^= crc_table[c];
	}
	crc=~crc;
	return crc;
}

uint32_t crc32init(void)
{
	uint32_t crc=(uint32_t)-1;
	return crc;
}

uint32_t crc32byte(uint8_t c, uint32_t crc_table[], uint32_t crc)
{
	c ^= (uint8_t) crc;
	crc >>= 8;
	crc ^= crc_table[c];
	return crc;
}

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
	data = (uint8_t*)malloc(size+1024);
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
	int priem;
	if(table_size<=16)
	{
		priem=3;
	}
	else if(table_size<=64)
	{
		priem=23;
	}
	else if(table_size<=128)
	{
		priem=47;
	}
	else if(table_size<=256)
	{
		priem=181;
	}
	else if(table_size<=512)
	{
		priem=881;
	}
	else if(table_size<=1024)
	{
		priem=8053;
	}
	else if(table_size<=2048)
	{
		priem=883;
	}
	else if(table_size<=4096)
	{
		priem=1607;
	}
	else if(table_size<=8192)
	{
		priem=18149;
	}
	else if(table_size<=16384)
	{
		priem=11807;
	}
	else if(table_size<=32768)
	{
		priem=1399;
	}
	else
	{
		priem=29333;
	}
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

#define MAX_TABLE (13)
char make_table_names[MAX_TABLE][128]=
{
	"entropy ",
	"huffman ",
	"rng     ",
	"simple  ",
	"s.simple",
	"precise ",
	"priem   ",
	"us. bh. ",
	"sb. bh. ",
	"ss. bh. ",
	"sb. bl. ",
	"ss. bl. ",
	"HDL     "
};

void make_table(int nr, int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
	switch(nr)
	{
	case 2:
	default:
		make_rng_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 3:
		make_simple_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 4:
		make_sorted_simple_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 5:
		make_precise_lff_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 6:
		make_priem_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 7:
		make_us_bh_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 8:
		make_sb_bh_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 9:
		make_ss_bh_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 10:
		make_sb_bl_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 11:
		make_ss_bl_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	case 12:
		make_hdl_table(ans_table, symbols_count, symbol_size, table_size);
		break;
	}
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
}

void encode_symbol(int symbol, mpz_t x, int ans_table[], int symbols_count[], int table_size, int s_offset[], int count_2_pos[])
{
	int r;
	int pos;
	mpz_add_ui(x, x, 1);
	r=(int)mpz_tdiv_q_ui(x, x, symbols_count[symbol]);
	mpz_mul_ui(x, x, table_size);
	mpz_add_ui(x, x, count_2_pos[s_offset[symbol]+r]);
}

void encode_ans_uint8_t(mpz_t x, uint8_t *data, unsigned long data_size, int ans_table[], int symbols_count[], int table_size)
{
	int s_offset[MAX_SYMBOL_COUNT];
	int count_2_pos[MAX_ANS_TABLE];
	{
		int i;
		int s_count[MAX_SYMBOL_COUNT]={0};
		for(i=0; i<table_size; i++)
		{
			int s=ans_table[i];
			s_count[s]++;
		}
		s_offset[0]=0;
		for(i=0; i<(MAX_SYMBOL_COUNT-1); i++)
		{
			s_offset[i+1]=s_offset[i]+s_count[i];
			s_count[i]=0;
		}
		s_count[MAX_SYMBOL_COUNT]=0;
		for(i=0; i<table_size; i++)
		{
			int s=ans_table[i];
			count_2_pos[s_offset[s]+s_count[s]]=i;
			s_count[s]++;
		}
	}
	while(data_size!=0)
	{
		data_size--;
		encode_symbol(data[data_size], x, ans_table, symbols_count, table_size, s_offset, count_2_pos);
	}
}

int decode_symbol(mpz_t x, int ans_table[], int symbols_count[], int table_size, int pos_2_count[])
{
	int symbol;
	int index;
	int pos;
	index=(int)mpz_tdiv_q_ui(x, x, table_size);
	symbol=ans_table[index];
	pos=pos_2_count[index];
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

uint32_t decode_ans_ascii_t(mpz_t x, int ans_table[], int symbols_count[], int table_size, uint32_t crc_table[])
{
	uint32_t crc=crc32init();
	int s_count[MAX_SYMBOL_COUNT]={0};
	int pos_2_count[MAX_ANS_TABLE];
	int i;
	for(i=0; i<table_size; i++)
	{
		int s=ans_table[i];
		pos_2_count[i]=s_count[s];
		s_count[s]++;
	}
	while(mpz_sgn(x)!=0)
	{
		int symbol;
		symbol=decode_symbol(x, ans_table, symbols_count, table_size, pos_2_count);
		crc=crc32byte(symbol, crc_table, crc);
//		printf("%c", symbol);
	}
	crc=~crc;
	return crc;
}

void make_clen(int node, int c_len[], int c_left[], int c_right[], int bits)
{
	if(node<MAX_SYMBOL_COUNT)
	{
		c_len[node]=bits;
	}
	else
	{
		make_clen(c_left[node], c_len, c_left, c_right, bits+1);
		make_clen(c_right[node], c_len, c_left, c_right, bits+1);
	}
}

uint64_t huffman_entropy(uint8_t *data, unsigned long data_size)
{
	int64_t symbol_count[2*MAX_SYMBOL_COUNT]={0};
	int c_left[2*MAX_SYMBOL_COUNT];
	int c_right[2*MAX_SYMBOL_COUNT];
	int c_len[2*MAX_SYMBOL_COUNT]={0};
	int node=MAX_SYMBOL_COUNT;
	if(data_size==0)
	{
		return 0;
	}
	{
		unsigned long i=data_size;
		while(i!=0)
		{
			i--;
			symbol_count[data[i]]++;
		}
	}
	for(;;)
	{
		uint64_t f0;
		int s0;
		uint64_t f1;
		int s1;
		int i;
		f0=data_size;
		for(i=0; i<2*MAX_SYMBOL_COUNT; i++)
		{
			if((symbol_count[i]>0) && (symbol_count[i]<=f0))
			{
				f0=symbol_count[i];
				s0=i;
			}
		}
		symbol_count[s0]=-f0;
		if(f0==data_size)
		{
			node=s0;
			break;
		}
		f1=data_size;
		for(i=0; i<2*MAX_SYMBOL_COUNT; i++)
		{
			if((symbol_count[i]>0) && (symbol_count[i]<f1))
			{
				f1=symbol_count[i];
				s1=i;
			}
		}
		symbol_count[s1]=-f1;
		c_left[node]=s0;
		c_right[node]=s1;
		symbol_count[node]=f0+f1;
		node++;
	}
	make_clen(node, c_len, c_left, c_right, 0);
	{
		uint64_t total_bits=0;
		int i;
		for(i=0; i<MAX_SYMBOL_COUNT; i++)
		{
			total_bits+=c_len[i]*(-symbol_count[i]);
		}
		return total_bits;
	}
	
}

uint64_t tans_entropy(uint8_t *data, unsigned long data_size, int symbols_count[], int symbol_size, int table_size)
{
	double symbol_bits[MAX_SYMBOL_COUNT]={0.0};
	int symbol;
	int table_bits=0;
	int tmp=table_size;
	tmp>>=1;
	while(tmp>0)
	{
		table_bits++;
		tmp>>=1;
	}
	for(symbol=0; symbol<symbol_size; symbol++)
	{
		int symbol_count=symbols_count[symbol];
		if(symbol_count>0)
		{
			int tmp=symbol_count;
			int bits=0;
			int total_bits=0;
			int x_bits;
			tmp--;
			while(tmp>0)
			{
				bits++;
				tmp>>=1;
			}
			x_bits=(1<<bits)-symbol_count;
			symbol_count-=x_bits;
			total_bits=symbol_count*(table_bits-bits);
			bits--;
			total_bits+=2*x_bits*(table_bits-bits);
			symbol_bits[symbol]=(double)total_bits/(double)(symbol_count+2*x_bits);
		}
	}
	{
		double bit_count=0.0;
		while(data_size!=0)
		{
			int symbol;
			data_size--;
			symbol=(int)data[data_size];
			bit_count+=symbol_bits[symbol];
		}
		return (uint64_t)ceil(bit_count);
	}
}		

void make_tans_table(int tans_table[], int tans_bits[], int ans_table[], int symbols_count[], int symbol_size, int table_size)
{
	int symbol;
	int table_bits=0;
	int tmp=table_size;
	tmp>>=1;
	while(tmp>0)
	{
		table_bits++;
		tmp>>=1;
	}
	for(symbol=0; symbol<symbol_size; symbol++)
	{
		int symbol_count=symbols_count[symbol];
		if(symbol_count>0)
		{
			int tmp=symbol_count;
			int pos=0;
			int bit_pattern;
			int bits=0;
			int x_bits;
			tmp--;
			while(tmp>0)
			{
				bits++;
				tmp>>=1;
			}
			x_bits=(1<<bits)-symbol_count;
			if(01)
			{
				if(x_bits==0)
				{
					bit_pattern=table_size-(1<<(table_bits-bits));
				}
				else
				{
					bit_pattern=table_size-(1<<(table_bits-bits+1));
				}
			}
			else
			{
				bit_pattern=0;
			}
			symbol_count-=x_bits;
			while(symbol_count>0)
			{
				if(ans_table[pos]==symbol)
				{
					bit_pattern&=table_size-1;
					tans_table[pos]=bit_pattern;
					tans_bits[pos]=table_bits-bits;
					bit_pattern+=1<<(table_bits-bits);
					symbol_count--;
				}
				pos++;
			}
			bits--;
			while(x_bits>0)
			{
				if(ans_table[pos]==symbol)
				{
					bit_pattern&=table_size-1;
					tans_table[pos]=bit_pattern;
					tans_bits[pos]=table_bits-bits;
					bit_pattern+=1<<(table_bits-bits);
					x_bits--;
				}
				pos++;
			}
			//printf("symbol = %02x, bit_pattern=%i\n", symbol, bit_pattern);
		}
	}	
}		

unsigned int get_bits(mpz_t x, int count)
{
	unsigned int bits=(unsigned int)mpz_tdiv_q_ui(x, x, 1<<count);
	return bits;
}

void put_bits(mpz_t x, unsigned int bits, int count)
{
	#ifdef TEST_MODE
	global_bit_count+=count;
	#else
	mpz_mul_ui(x, x, 1<<count);
	mpz_add_ui(x, x, bits);
	#endif
}

uint32_t decode_tans(mpz_t x, unsigned long data_size, int ans_table[], int symbols_count[], int table_size, int symbol_size, uint32_t crc_table[])
{
	int tans_table[MAX_ANS_TABLE];
	int tans_bits[MAX_ANS_TABLE];
	unsigned int state;
	unsigned int state_offset=0;
	uint32_t crc=crc32init();
	make_tans_table(tans_table, tans_bits, ans_table, symbols_count, symbol_size, table_size);
	table_size--;
	{
		int table_bits=0;
		int tmp=table_size;
		while(tmp>0)
		{
			table_bits++;
			tmp>>=1;
		}
		state=get_bits(x, table_bits);
	}
	while(--data_size!=0)
	{
		int symbol;
		symbol=ans_table[state];
		state=tans_table[state]+get_bits(x, tans_bits[state]);
		state_offset+=DELTA_STATE_OFFSET;
		state+=state_offset;
		state&=table_size;
		crc=crc32byte(symbol, crc_table, crc);
//		printf("%c", symbol);
	}
	{
		int symbol=ans_table[state];
		crc=crc32byte(symbol, crc_table, crc);
//		printf("%c", symbol);
	}
	crc=~crc;
	return crc;
}

void encode_tans(mpz_t x, uint8_t *data, unsigned long data_size, int ans_table[], int symbols_count[], int table_size, int symbol_size, uint32_t crc_table[])
{
	int tans_table[MAX_ANS_TABLE];
	int tans_bits[MAX_ANS_TABLE];
	unsigned int state=0;
	unsigned int state_offset=data_size*DELTA_STATE_OFFSET;
	make_tans_table(tans_table, tans_bits, ans_table, symbols_count, symbol_size, table_size);
	if(data_size!=0)
	{
		int symbol;
		data_size--;
		symbol=(int)data[data_size];
		while(ans_table[state]!=symbol)
		{
			state++;
		}
	}
	while(data_size!=0)
	{
		int symbol;
		int new_state=0;
		data_size--;
		symbol=(int)data[data_size];
		state_offset-=DELTA_STATE_OFFSET;
		state-=state_offset;
		state&=table_size-1;
		for(;;)
		{
			if(ans_table[new_state]==symbol)
			{
				if((tans_table[new_state]<=state) && ((tans_table[new_state]+(1<<tans_bits[new_state]))>state))
				{
					break;
				}
			}
			new_state++;
		}
		{
			int bits;
			bits=1<<tans_bits[new_state];
			bits--;
			bits&=state;
			put_bits(x, bits, tans_bits[new_state]);
			state=new_state;
		}
	}
	{
		int table_bits=0;
		table_size>>=1;
		while(table_size>0)
		{
			table_bits++;
			table_size>>=1;
		}
		put_bits(x, state, table_bits);
	}
}

void statistiek(int table_size, int64_t freq[], int symbols_count[], int table[], int symbol_size, int64_t size)
{
	int i;
	float total_bits=0;
	for(i=0; i<symbol_size; i++)
	{
		if(freq[i]!=0)
		{
			float bits;
			bits=log2f((float)table_size/(float)symbols_count[i]);
			total_bits+=bits*freq[i];
			#if 0
			if((i>32) && (i<127))
			{
				printf("freq[%4c] = %3li = %4i = %.1f bits\n", i, freq[i], symbols_count[i], bits);
			}
			else
			{
				printf("freq[0x%02X] = %3li = %4i = %.1f bits\n", i, freq[i], symbols_count[i], bits);
			}
			#endif
		}
	}
	printf("Total_bits %.1f = %.1f, ", (float)size*log2f(symbol_size), total_bits);
	#if 0
	for(i=0; i<table_size; i++)
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
	if(0)
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
					if(table[j]==i)
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

#define TABLE_SIZE 4096
#define SYMBOL_SIZE 256
#define BIG_DATA_SIZE (256*1024)

int main(int argc, char* argv[])
{
	uint32_t crc_table[256];
	int i;
	i=1;
	make_crc32_table(crc_table);
	double ratio[MAX_TABLE]={0.0};
	while(i<argc)
	{
		uint32_t crc;
		int64_t freq[SYMBOL_SIZE];
		int symbols_count[SYMBOL_SIZE];
		int table[MAX_ANS_TABLE];
		int64_t size;
		uint8_t* data;
		uint8_t* big_data;
		uint8_t* data_to_be_compressed;
		int symbol_size=SYMBOL_SIZE;
		big_data=(uint8_t*)malloc(BIG_DATA_SIZE);
		if(big_data==NULL)
		{
			printf("Malloc error big_data!\n");
			exit(-1);
		}
		printf("Loading file %s\n",argv[i]);
		size=load_file(argv[i], &data);
		if(data==NULL)
		{
			printf("Error loading %s\n", argv[i]);
			return -1;
		}
		freq_count(data, size, freq, symbol_size);
		if(0)
		{ /* vul big data met random bytes met symbols_count frequentie, we willen altijd de zelfde random data hebben */
			int i;
			verdeel_symbols(freq, symbols_count, symbol_size, 65536);
			make_table(1, table, symbols_count, symbol_size, 65536);
			verdeel_symbols(freq, symbols_count, symbol_size, TABLE_SIZE);
			srand48(31415);
			for(i=0; i<BIG_DATA_SIZE; i++)
			{
				big_data[i]=table[lrand48()&0x7FFF];
			}
			size=BIG_DATA_SIZE;
			data_to_be_compressed=big_data;
			freq_count(data_to_be_compressed, size, freq, symbol_size);
		}
		else
		{
			verdeel_symbols(freq, symbols_count, symbol_size, TABLE_SIZE);
			data_to_be_compressed=data;
		}
		crc=crc32(size, data_to_be_compressed, crc_table);
		printf("File size = %li\n", size);
		{
			int table_no;
			int best_table=0;
			int64_t best_size=size*8;
			{
				uint64_t bits;
				uint64_t huff_bits;
				bits=tans_entropy(data, size, symbols_count, symbol_size, TABLE_SIZE);
				huff_bits=huffman_entropy(data, size);
				ratio[0]+=(double)bits/(double)size;
				ratio[1]+=(double)huff_bits/(double)size;
				printf("Total_bits %li, entropy %li, huffman %li\n", size*8, bits, huff_bits);
			}
			for(table_no=2; table_no<MAX_TABLE; table_no++)
			{
				make_table(table_no, table, symbols_count, symbol_size, TABLE_SIZE);
				printf("Table: %s, table_size=%i ", make_table_names[table_no], TABLE_SIZE);
				statistiek(TABLE_SIZE, freq, symbols_count, table, symbol_size, size);
				{ /* encoding and decoding */
					mpz_t x;
					mpz_init(x);
					#ifdef TEST_MODE
					global_bit_count=0;
					#endif
					uint32_t new_crc;
					unsigned long bits;
//					printf("Result:\n");
					encode_tans(x, data_to_be_compressed, size, table, symbols_count, TABLE_SIZE, symbol_size, crc_table);
//					encode_ans_uint8_t(x, data_to_be_compressed, size, table, symbols_count, TABLE_SIZE);
					#ifdef TEST_MODE
					bits=global_bit_count;
					#else
					bits=mpz_sizeinbase (x, 2);
					#endif
					printf("compressed = %li bits, ", bits);
					if(bits<best_size)
					{
						best_size=bits;
						best_table=table_no;
					}
					ratio[table_no]+=(double)bits/(double)size;
					#ifndef TEST_MODE
               new_crc=decode_tans(x, size, table, symbols_count, TABLE_SIZE, symbol_size, crc_table);
//               new_crc=decode_ans_ascii_t(x, table, symbols_count, TABLE_SIZE, crc_table);
					if(crc!=new_crc)
					{
						printf("CRC error!\n");
						exit(-1);
					}
					else
					{
						printf("CRC OK\n");
					}
					#else
					printf("No CRC\n");
               #endif
					mpz_clear(x);
				}
			}
			printf("Best table: %s, %li bits\n", make_table_names[best_table], best_size);
		}
		free(big_data);
		free(data);
		i++;
	}
	{
		int table_no;
		int best=0;
		double best_ratio=ratio[best];
		for(table_no=0; table_no<MAX_TABLE; table_no++)
		{
//			printf("table: %s = %lf\n", make_table_names[table_no], ratio[table_no]);
			if(ratio[table_no]<best_ratio)
			{
				best_ratio=ratio[table_no];
				best=table_no;
			}
		}
		table_no=0;
		do
		{
			double second_best=1e9;
			int sb=0;
			int i;
			int current=table_no+1;
			for(i=0; i<MAX_TABLE; i++)
			{
				if(ratio[i]==best_ratio)
				{
					printf("%02i: %s = %lf\n", current, make_table_names[i], best_ratio);
					table_no++;
				}
			}
			for(i=0; i<MAX_TABLE; i++)
			{
				if((ratio[i]<second_best) && (ratio[i]>best_ratio))
				{
					second_best=ratio[i];
					sb=i;
				}
			}
			best_ratio=second_best;
			best=sb;
		}	while(table_no<MAX_TABLE);
	}
	return 0;
}
