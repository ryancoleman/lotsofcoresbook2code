#define _GNU_SOURCE
#include<stdio.h>
#include<pthread.h>
#include<sched.h>

void *dosomething(void *data) {

int i=100;
	while(1);

	return NULL;

}

void main()
{

pthread_t thread1, thread2;

cpu_set_t cs1, cs2;
int rc;
int status;
rc = pthread_create(&thread1, NULL, dosomething, NULL);
	
	if (rc) {
				printf("ERROR; return code from pthread_create() is %d\n", rc);
		
	} 
	else {
				//cout << "Thread Created Successfully\n";
					}
CPU_ZERO (&cs1);
    CPU_SET (6, &cs1);
    status = sched_setaffinity(0, sizeof(cs1), &cs1);
    if (status < 0) {
	printf ("Error: unable to bind thread to core");
    }

/*rc = pthread_create(&thread2, NULL, dosomething, NULL);

CPU_ZERO(&cs2);
CPU_SET(5, &cs2);
status = sched_setaffinity(0, sizeof(cs2), &cs2);



*/
pthread_join(thread1, NULL);
//pthread_join(thread2, NULL);

    printf("Bind successful to core 1");


}
