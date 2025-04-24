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

void verdeel_bits(int64_t freq[], int bits_count[], int symbol_size, int total_bits)
{
	int i;
	bzero(bits_count, symbol_size*sizeof(int));
	for(i=0; i<symbol_size; i++)
	{
		if(freq[i]>0)
		{
			bits_count[i]=1;
			total_bits--;
		}
	}
	while(total_bits!=0)
	{
		int i;
		float worst=0;
		int worst_i=0;
		for(i=0; i<symbol_size; i++)
		if(freq[i]>0)
		{
			float cost;
			cost=(float)freq[i]/(float)bits_count[i];
			if(cost>worst)
			{
				worst=cost;
				worst_i=i;
			}
		}
		bits_count[worst_i]++;
		total_bits--;
	}
}

#define BITS_COUNT 4096
#define SYMBOL_SIZE 256

int main(int argc, char* argv[])
{
	int i;
	i=1;
	while(i<argc)
	{
		int64_t freq[SYMBOL_SIZE];
		int bits_count[BITS_COUNT];
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
		verdeel_bits(freq, bits_count, SYMBOL_SIZE, BITS_COUNT);
		{
			int i;
			float total_bits=0;
			for(i=0; i<256; i++)
			{
				if(freq[i]!=0)
				{
					float bits;
					bits=log2f((float)BITS_COUNT/(float)bits_count[i]);
					total_bits+=bits*freq[i];
					if((i>32) && (i<127))
					{
						printf("freq[%4c] = %3li = %4i = %f bits\n", i, freq[i], bits_count[i], bits);
					}
					else
					{
						printf("freq[0x%02X] = %3li = %4i = %f bits\n", i, freq[i], bits_count[i], bits);
					}
				}
			}
			printf("Total_bits %li = %f\n", size*8, total_bits);
		}
		free(data);
		i++;
	}
	return 0;
}
