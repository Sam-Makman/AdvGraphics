#include <time.h>
#include <stdlib.h>
#include <stdio.h>

int main(){
	int n = 1000;
	int i;
	int count = 0;
	srand(time(NULL));
	for(i=0 i< n; i++){
		int knight = rand()*8;
		printf("%d\n",knight );
	}

}