#include <stdio.h>
#include <stdlib.h>

void writeInputMatrix(char *fileName, float value, int width)
{
	FILE *fp = fopen(fileName, "w");
	int i, j;

	for(i=0; i<width; i++)
	{
		for(j=0; j<width; j++)
		{
			fprintf(fp, "%f ", value);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

int main(int argc, char **argv)
{
	int width;
	float *mat;
	float value;
	int i, j;

	if(argc != 4)
	{
		printf("Incorrect Usage!\n");
		printf("./writeInputMatrix <fileName> <value> <width>\n");
		exit(1);
	}

	width = atoi(argv[2]);
	sscanf(argv[3], "%f", &value);

	writeInputMatrix(argv[1], width, value);

	return 0;
}
