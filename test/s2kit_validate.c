#include <stdio.h>


#include <s2kit/FST_semi_memo.h>
#include <sphradiometer/sh_series.h>


int main(int argc, char *argv[])
{
	int l_max = 10;
	int l, m;

	for(l = 0; l <= l_max; l++)
		for(m = -l; m <= l; m++)
			if(seanindex(m, l, l_max + 1) != sh_series_params_lmoffset(l_max, l, m)) {
				fprintf(stderr, "seanindex(%d, %d, %d) = %d but sh_series_params_lmoffset(%d, %d, %d) = %d\n", m, l, l_max + 1, seanindex(m, l, l_max + 1), l_max, l, m, sh_series_params_lmoffset(l_max, l, m));
				exit(1);
			}

	return 0;
}
