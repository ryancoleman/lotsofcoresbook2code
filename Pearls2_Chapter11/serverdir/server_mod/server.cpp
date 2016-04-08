#define _MULTI_THREADED
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include "defs.h"
#include <sys/sendfile.h>
#include <iostream>
#include <pthread.h>
#include <flann/flann.h>
#include <time.h>
#include <flann/flann.hpp>
#include <assert.h>
#include <sys/time.h>
#include <signal.h>
#include <sched.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include <sys/mman.h>

#define TOP_FIVE 5
#define TOP_FOUR 4
#define MAX_CORES 256
#define NUM_PROCESSES     16
#define MAX_FILE_NAME_LENGTH 60
#define SIZE_LINE_ARRAY 5
#define MAX_BUFFER 100	
#define MAX_IMAGES 400
#define MAX_FILE_LEN 1000
#define MAX_SIFT_CMD_LEN 5000
#define PORT_NO 27015
#define DIM 128
#define MAX_KEYPOINTS		1000
#define LOWE_MATCH_RATIO	0.6
#define NNIDX_0				0
#define NNIDX_1				1
#define NN					3
#define MATCH_PERCENT		0.3
#define EPSILON 3.0
#define RANDOM_WALK	0
#define MAX_BUFFER_SIZE 1000
#define CPU_0 4
#define ARGCOUNT 4
#define LISTENPORT 5
#define BUFFER 255
#define MAX_MSG 18
#define NUM_TREES 8
#define NUM_CHECKS 32
#define PRECISION 0.6
#define VERY_LARGE_NUMBER 99999
#define BUFFERSIZE 256

#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define BLAH alloc_if(0) free_if(0)

using namespace std;

//struct to hold query image info
struct packet {
	int num_keys;
	int total_query_images;
	int img_id;
}__attribute__((packed));

// __attribute__((target(mic))) struct packet client_packet;
struct packet client_packet;

char *query_images[MAX_IMAGES];
char *query_key_names[MAX_IMAGES];
char *query_key_names_structured[MAX_IMAGES];
// __attribute__((target(mic))) float *hdf5_query_files;
// __attribute__((target(mic))) float *hdf5_data_files;
float *hdf5_query_files;
float *hdf5_data_files[MAX_IMAGES];
float *proxy_data_files;
int key_count_query[MAX_IMAGES];
float *query_file;

int t_sum = 0;
int t_create = 0;
int t_search = 0;
int th_create = 0;
int file_size;

char *data_images[MAX_IMAGES];
char *data_key_names[MAX_IMAGES];
char *data_key_names_structured[MAX_IMAGES];

int num_query_images = 0;
int num_data_images = 0;
int g = 0;
int total_cpus_present;
int KDT,FT;
Keypoint keydarr[MAX_IMAGES];
Keypoint keyqarr[MAX_IMAGES];

const char *sift1 = "./sift <";
const char *sift2 = " >";
const char *sift3_2 = "/mnt/pmfs/data/";
const char *sift4 = ".key";

//cmd line params
int thread_count;
char db_path_query[MAX_FILE_LEN];
char db_path_data[MAX_FILE_LEN];

char buffer[MAX_BUFFER_SIZE];
int keys_data[MAX_IMAGES] = { 0 };

//struct to hold matched image data
typedef struct img_data {
	int num;
	int match_found;
} im_data;

//struct to hold data to be passed to a thread
typedef struct str_thdata {
	int threadcount;
	int image_count;
	int from;
	int to;
	int top_four;
	int from_value;
	int to_value;

} thdata;

int BindToCpu(int cpu_num);

int flann_threads = 1;

thdata th_data[MAX_CORES];
// __attribute__((target(mic))) im_data imageInfo[MAX_IMAGES];
// __attribute__((target(mic))) int keys_data[MAX_IMAGES] = { 0 };
im_data imageInfo[MAX_IMAGES];

/*
 * 	Function: error
	Parameters: const char *msg - Error message to be printed
	Purpose: Prints error message
	Returns: none
 */
void error(const char *msg) {
	perror(msg);
	exit(1);
}

/* qsort struct comparision function */
/*
 * 	Function: int_cmp
	Parameters: const void *a - pointer to one matched image structure
				const void *b - pointer to one matched image structure
	Purpose: Qsort structure comparison function. Function helps sort the macthed image structure img_data to obtain top 4 matches found.
			 It returns negative if b>a and positive if a>b 
	Returns: int
 */
int int_cmp(const void *a, const void *b) {
	im_data *ia = (im_data*) a;
	im_data *ib = (im_data *) b;
	return (ib->match_found - ia->match_found);

}

//Allocation of memory
/*
 * 	Function: init_server
	Parameters: none
	Purpose: Mallocs memory for certain variables used throughout the program
	Returns: void
 */
void init_server() {
	int i;

	for (i = 0; i < MAX_IMAGES; i++) {
		hdf5_data_files[i] = (float *) malloc(MAX_IMAGES);
		data_images[i] = (char *) malloc(MAX_FILE_LEN);
		data_key_names[i] = (char *) malloc(MAX_FILE_LEN);
		data_key_names_structured[i] = (char *) malloc(MAX_FILE_LEN);
		query_images[i] = (char *) malloc(MAX_FILE_LEN);
		query_key_names[i] = (char *) malloc(MAX_FILE_LEN);
		query_key_names_structured[i] = (char *) malloc(MAX_FILE_LEN);
	}
	
		hdf5_query_files = (float *) malloc(MAX_IMAGES * MAX_IMAGES);
		proxy_data_files = (float *) malloc(MAX_IMAGES * MAX_IMAGES);
		query_file= (float*) malloc(MAX_IMAGES);
}

//Function which mmaps the key points
/*
 *  Function: read_points_q
	Parameters: const char* filename - location of the structured key file of the database image
	Purpose: This function mmaps the structured binary key files from storage areas specified in the command line param (PMFS, PMBD, SATA) into the memory.
			 It returns a pointer to the new mapping created in the virtual address space of the process.
	Returns: float* (pointer to float array)
 */
float* read_points_q(const char* filename) {
	float* data;
	float *p;
	int fin;
	int i, j;

	fin = open(filename, O_RDWR);

	if (!fin) {
		printf("Cannot open input file.\n");
		exit(1);
	}

	struct stat sb;
	size_t len;
	fstat(fin, &sb);

	//printf("Size of %s : %lu \n", filename, len);

	float *file_memory = (float*) mmap(NULL, sb.st_size,
			PROT_READ | PROT_WRITE, MAP_SHARED, fin, 0);
	file_size = sb.st_size;
	close(fin);

	return file_memory;
}

//Function which generates the structured key files from the key files and
//mmaps the database images key points into suitable storage options
/*  Function: generateKey_mmap
	Parameters: int count - no of database images
				int runsift - whether the SIFT needs to run or not (0/1) on the database images. (Supplied through the command line param)
	Purpose: This function generated the structured key files from the .key files produced from the SIFT algorithm. The SIFT algorithm produces .key files which has too many paramters. This function creates a file of only the 128 bit features. FLANN algorithm needs the input in this format.
	Returns: void
 */
inline void generateKey_mmap(int count, int runsift) {

	int i, j, p;
	float p1;
	Keypoint k1 = NULL;
	FILE *fp;
	Keypoint k;
	unsigned char *pk1;

	char key_file[MAX_FILE_LEN];

	timeval start;
	timeval end;
	gettimeofday(&start, NULL);

	float elapsed = 0;

	for (i = 0; i < count; i++) {

		int key_count = 0;
		sprintf(key_file, "%s%d", db_path_data, i);
		sprintf(data_key_names[i], "%s%d%s", db_path_data, i, sift4);
		//sprintf(data_key_names_structured[i], "%s%d", db_path_data, i);

		FILE *key;
		key = fopen(data_key_names[i], "r");
		fscanf(key, "%d", &keys_data[i]);
		fclose(key);

		if (runsift) {
			fp = fopen(key_file, "w");
			k1 = ReadKeyFile(data_key_names[i]);

			for (k = k1; k != NULL; k = k->next) {
				key_count++;

				pk1 = k->descrip;
				for (j = 0; j < DIM; j++) {
					p1 = *pk1++;
					fwrite(&p1, sizeof(float), 1, fp);
				}
			}
			pk1 = NULL;
			fclose(fp);
		}

		hdf5_data_files[i] = read_points_q(data_key_names_structured[i]);
		float * tmp = read_points_q(data_key_names_structured[i]);
		memcpy(&proxy_data_files[MAX_IMAGES * i], tmp, file_size);
	}
	gettimeofday(&end, NULL);

	double elapsedTime;
	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0; // sec to ms
	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0; // us to ms

	printf("Total time in init thread: %g ms\n", elapsedTime);
}

//The db init module which creates all the key files and then creates them in the
//format suitable for FLANN function. It also mmaps these points into memory defined
//in the command line parameter (SATA, PMFS, PMBD)
/*
 *  Function: db_init
	Parameters: char *dbimgpath - database file having the locations of all the database images
				int runsift - whether the SIFT needs to run or not (0/1) on the database images. (Supplied through the command line param)
	Purpose: This function gets the locations of all the database images and stores
	 	 	 them into an array. If runsift parameter is 1, then it runs the SIFT 
	 	 	 algorithm to kind the key points of all the database images. 
	 	 	 It stores them as sift key files. (This parameter needs to be 1 
	 	 	 everytime the database images have changes or being loaded for 
	 	 	 the first time)
	 	 	 It then calls generateKey_mmap to map the key files to process's memory from specified storage areas (PMFS, PMBD, SATA)
	Returns: void
 */
void db_init(char *dbimgpath, int runsift) {

	init_server();

	FILE *fp;
	fp = fopen(dbimgpath, "r");
	int i = 0;
	char *c;
	num_data_images = 0;
	if (fp != NULL) {
		while (fgets(data_images[i], MAX_BUFFER, fp) != NULL) {
			c = strchr(data_images[i], '\n');
			if (c) {
				*c = '\0';
			}
			num_data_images++;
			i++;
		}
	}

	else {
		cout << "Error opening file, check location\n";
	}
	fclose(fp);
	cout<<"Number of data images"<<num_data_images<<endl;

	//generate sift key files for the data_images
	char siftcommand[MAX_SIFT_CMD_LEN];

	for (int i = 0; i < num_data_images; i++) {

		sprintf(siftcommand, "%s%s%s%s%d%s", sift1, data_images[i], sift2,
				db_path_data, i, sift4);
		sprintf(data_key_names[i], "%s%d%s", db_path_data, i, sift4);
		if (runsift) {
			system(siftcommand);
		}

	}
	/*
	 //restructure the file into ANN format
	 FILE *shellfp;
	 shellfp = popen("./../findnum_images.sh", "r");
	 char result[10];
	 fscanf(shellfp, "%s", result);*/

	//num_data_images = atoi(result);
	//num_data_images = 7794;

	timeval start;
	timeval end;
	gettimeofday(&start, NULL);

	generateKey_mmap(num_data_images, runsift);

	gettimeofday(&end, NULL);

	double elapsedTime;

	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0; // sec to ms
	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0; // us to ms

}

//Function which processes the query image data obtained
/*
 *  Function: process_query_data
	Parameters: none
	Purpose: Post processes the query image data obtained from the client. It also calls the read_points_q to mmap the query files too.
	Returns: void
 */
void process_query_data() {

	Keypoint k;
	unsigned char *pk1;
	Keypoint k1 = NULL;
	float p1;

	num_query_images = client_packet.total_query_images;

	char key_file[MAX_FILE_LEN];
	float *temp;
	int total = DIM * client_packet.num_keys;

	temp = (float*) malloc(total * sizeof(float));
	sprintf(query_key_names_structured[client_packet.img_id], "%s%d",
			db_path_query, client_packet.img_id);

	FILE *fp;

	k1 = ReadKeyFile(query_key_names[client_packet.img_id]);

	fp = fopen(query_key_names_structured[client_packet.img_id], "w");
	if (fp == NULL) {
		cout << "Error opening file\n";
	}

	int iterate = 0;
	int keep_count = 0;

	for (k = k1; k != NULL; k = k->next) {
		pk1 = k->descrip;
		iterate++;
		for (int j = 0; j < DIM; j++) {
			p1 = *pk1++;
			fwrite(&p1, sizeof(float), 1, fp);
		}
	}
	client_packet.num_keys = iterate;

	pk1 = NULL;
	fclose(fp);

	float * tmp = read_points_q(query_key_names_structured[client_packet.img_id]);
	memcpy(&hdf5_query_files[MAX_IMAGES * client_packet.img_id], tmp, file_size);
	// assert(hdf5_query_files[MAX_IMAGES * client_packet.img_id] != NULL);
	memcpy(query_file, &hdf5_query_files[MAX_IMAGES * client_packet.img_id], MAX_IMAGES);
}

//Fucntion to write results to files. Unused fucntion. Used for debugging only
/*
 *  Function: write_results
	Parameters: const char* filename - Fielname to write the results to
				int *data - Return array from the FLANN search function
				int rows - No of keys in the client query image
				int cols - Dimension of the key data (128)
				float *dists - Distances obtained from the FLANN search function
				int image_no - Image id of the query image whose results are being written
	Purpose: This function is for debugging purposes only to write the image results locally to a file.
	Returns: void
 */
void write_results(const char* filename, int *data, int rows,
		float *dists, int image_no) {
	FILE* fout;
	float *p;
	int i, j;

	fout = fopen(filename, "a+");
	if (!fout) {
		printf("Cannot open output file.\n");
		exit(1);
	}

	int matched = 0;
	p = dists;
	for (i = 0; i < rows; ++i) {
		if (*p < (LOWE_MATCH_RATIO * *(p + 1))) {

			matched++;
		}
		p = p + 2;
	}
	cout<<"\nMatched"<<matched;
	//int xc;
	//cin>>xc;
	if (!(matched < 0.15 * keys_data[image_no])) {

		fprintf(fout, "Query Image match count with: %d \n", image_no);
		fprintf(fout, "Matches: %d \n", matched);
	}
	fclose(fout);

}

//Function which is the actual FLANN create and search function. Here each image
//index is created and then the query points are searched against this index. 
//The correct matches are also determined in this function. And the top 4 matches 
//are collected
/*
 * 	Function: call_flannC
	Parameters: int from - Since each thread is given subset of images. from indicates the start val of that image subset.
				int to - to indicates the to val of the image subset. 
				(from and to range from 0 to total_number_database_images)
				int from_value - Each thread maintains top 4 image matches it found. This from_value is an indicator of the array index it needs to store the top 4 matched image data.
				int to_value - end value of the array index it needs to store top 4 matched image data.
				(these integer are used as array indices of the img_data structure which stores all the information about the matches images)
	Purpose: This function does the actual image search.
			 IT DOES THE FOLLOWING THINGS FOR EACH THREAD AND ITS ASSIGNED SUBSET OF IMAGES:
				A. Build kd tree for each database image in its subset
				b. Search query key points in that kd tree
				c. Determines how many matches were obatined for each database image in its subset
				d. Obtains the top 4 matches once it finishes processing all the images in its subset.
	Returns: void
 */
void call_flannC(int id, int from, int to, int from_value, int to_value) {

	timeval start;
	timeval end;

	timeval createstart;
	timeval createend;
	double el_create;
	double sum_create = 0.0;

	timeval searchstart;
	timeval searchend;
	double el_search;
	double sum_search = 0.0;
	int images_searched = 0;

	gettimeofday(&start, NULL);

	float elapsed = 0;
	float *temp;

	int* result;
	flann_index_t index_id;
	struct FLANNParameters p;
	float* dists;

	p = DEFAULT_FLANN_PARAMETERS;
	p.algorithm = FLANN_INDEX_KDTREE;
	p.trees = KDT;
	p.checks = 32;
	p.target_precision = 0.6;
	p.cores = FT;
	cout<<"Flann_Thread_Cores  "<<p.cores;

	for (int i = from; i < to; i++) {
		result = (int*) malloc(client_packet.num_keys * NN * sizeof(int));
		dists = (float*) malloc(client_packet.num_keys * NN * sizeof(float));

		float speedup;
		gettimeofday(&createstart, NULL);
		//index_id = flann_build_index(hdf5_data_files[i], keys_data[i], DIM, &speedup, &p);
		index_id = flann_build_index(&proxy_data_files[MAX_IMAGES * i], keys_data[i], DIM, &speedup, &p);
		gettimeofday(&createend, NULL);
		// cout << "Cores= " << p.cores << " Speedup=" << speedup << endl;

		el_create = (createend.tv_sec - createstart.tv_sec) * 1000.0; // sec to ms
		el_create += (createend.tv_usec - createstart.tv_usec) / 1000.0; // us to ms
		sum_create += el_create;

		gettimeofday(&searchstart, NULL);
		//memcpy(query_file, &hdf5_query_files[MAX_IMAGES * client_packet.img_id], file_size);
		flann_find_nearest_neighbors_index(index_id,
		//		 &hdf5_query_files[MAX_IMAGES * client_packet.img_id], client_packet.num_keys,
			 	query_file, client_packet.num_keys,
				result, dists, NN, &p);
		gettimeofday(&searchend, NULL);

		el_search = (searchend.tv_sec - searchstart.tv_sec) * 1000.0; // sec to ms
		el_search += (searchend.tv_usec - searchstart.tv_usec) / 1000.0; // us to ms
		sum_search += el_search;
		images_searched++;

		temp = dists;
		int matched = 0;
		for (int s = 0; s < client_packet.num_keys; ++s) {
			if (*temp < (LOWE_MATCH_RATIO * *(temp + 1))) {
				matched++;
//				printf("\nMatching Matching Matching");
			}
			temp = temp + 2;
		}

        	//cout << "Thread: " << id << ". Image: " << i << ". Matched images: " << matched << endl;
		int min = VERY_LARGE_NUMBER;

		if (!(matched < 0.15 * client_packet.num_keys)) {
			if (matched < min) {
				min = matched;
			}

			if (from_value < to_value) {
				imageInfo[from_value].num = i;
				imageInfo[from_value].match_found = matched;
				from_value++;
			}

			else {
				size_t structs_len = sizeof(imageInfo) / sizeof(im_data);
				qsort(imageInfo, structs_len, sizeof(im_data), int_cmp);

				if (imageInfo[to_value - 1].match_found < matched) {
					imageInfo[to_value - 1].num = i;
					imageInfo[to_value - 1].match_found = matched;
				}
			}
			//printf("Exiting the match found loop %d\n",matched);
		}

	//	printf("tree searched %d\n", i);	
		//char result_file[MAX_FILE_LEN];
		//sprintf(result_file, "/home1/02857/prash/intel-img/visualsearch_orig/serverdir_phi/results/results%d.dat",pthread_self());
		//write_results(result_file, result, client_packet.num_keys, dists, i);

		flann_free_index(index_id, &p);
		free(result);
		free(dists);
	}

	gettimeofday(&end, NULL);

	double elapsedTime;

	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0; // sec to ms
	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0; // us to ms

	printf("Thread %d: Total time for create/thread: %g ms\n", id, sum_create);
	printf("Thread %d: Total time for search/thread: %g ms\n", id, sum_search);

	printf("Thread %d: Total time: %g ms\n", id, sum_create + sum_search);
	cout << "Thread " << id << ": Images searched = " << images_searched << endl << endl;
	//t_sum += elapsedTime ;
	t_create += sum_create;
	t_search += sum_search;

}

//Function which binds the threads to cores
/*
 * 	Function: BindToCpu
	Parameters: int cpu_num - core number (which core the thread should be bound to)
	Purpose: Function which binds threads to specified core. Affinitizes thread to a core.
	Returns: int
 */
int BindToCpu(int cpu_num) {
	cpu_set_t cs;
	CPU_ZERO(&cs);
	CPU_SET(cpu_num, &cs);
	pthread_t current_thread = pthread_self();
	int status;

	status = pthread_setaffinity_np(current_thread, sizeof(cs), &cs);
	if (status < 0) {
		printf("Error: unable to bind thread to core %d\n", cpu_num);
		perror(0);
		exit(1);
	}

	status = pthread_getaffinity_np(current_thread, sizeof(cs), &cs);

	for (int i = 0; i < CPU_SETSIZE; i++) {
		if (CPU_ISSET(i, &cs)) {
		}
	}
	return 1;
}

//Function which the threads calls to call the FLANN main functions for image search
//Hear the threads are also affinitized to different cores of the machine
/*
 * 	Function: call_image_matching_module
	Parameters: void *ptr - start_routine argument passed during pthread_create
	Purpose: This function first of all calls BindToCpu to affinitize each thread of a specific core. Then it calls call_flannC which does the FLANN building of kd tree and tree searching.
	Returns: NULL
 */
void *call_image_matching_module(void *ptr) {

	thdata *data;
	data = (thdata*) ptr;

	BindToCpu(data->threadcount);

	call_flannC(data->threadcount, data->from, data->to, data->from_value, data->to_value);

	return NULL;

}

//Function which sends all the image results back to the client.
/*
 *  Function: send_results_to_client
	Parameters: int sock - socket id
	Purpose: It does two things:
			 First it sends the number of result matches that the server has found for eacgh query image sent. (Top 4 matches for the query images are sent from the server).
			 If no matches are found, it suitably prints that no matches were found on the server database and also sends that number to the client.
			 Next, it send the actual key file information to the client on the TCP socket.
			 The key files of all the top 4 matches are sent to the client.
	Returns: void
 */
void send_results_to_client(int sock) {
	//get result number, defualt is top 4 image matches
	int image_matches = 0;
	for (int i = 0; i < TOP_FOUR; i++) {
		if (imageInfo[i].match_found > 0)
			image_matches++;
	}

	int converted_number = htonl(image_matches);
	if (image_matches == 0) {
		printf("\nFound no matches in server database\n");
		write(sock, &converted_number, sizeof(converted_number));
		return;
	}

	else {
		printf("\nSending top %d matches in server database\n", image_matches);
		write(sock, &converted_number, sizeof(converted_number));
	}

	int *sizes_images;
	sizes_images = (int*) malloc(TOP_FOUR * sizeof(int));

	/* Get file stats */
	for (int i = 0; i < image_matches; i++) {
		int fd;
		struct stat file_stat;
		fd = open(data_key_names[imageInfo[i].num], O_RDONLY);
		if (fd == -1) {
			printf("Error\n");
			exit(0);
		}
		if (fstat(fd, &file_stat) < 0) {
			printf("Error\n");
			exit(0);
		}

		//fprintf(stdout, "File Size: \n%d bytes\n", file_stat.st_size);

		sizes_images[i] = file_stat.st_size;

		char file_size[BUFFERSIZE];
		ssize_t len;
		sprintf(file_size, "%d", file_stat.st_size);
		/* Sending file size */

		len = send(sock, file_size, sizeof(file_size), 0);
		if (len < 0) {
			fprintf(stderr, "Error on sending greetings");
			exit(0);
		}

		//fprintf(stdout, "Server sent %d bytes for the size\n", len);
		close(fd);
	}

	for (int i = 0; i < image_matches; i++) {
		int fd;
		fd = open(data_key_names[imageInfo[i].num], O_RDONLY);

		off_t offset;
		int remain_data;
		int sent_bytes = 0;
		offset = 0;
		remain_data = sizes_images[i];
		/* Sending file data */
		while (((sent_bytes = sendfile(sock, fd, &offset, BUFFERSIZE)) > 0)
				&& (remain_data > 0)) {
			//fprintf(stdout, "1. Server sent %d bytes from file's data, offset is now : %d and remaining data = %d\n", sent_bytes, offset, remain_data);
			remain_data -= sent_bytes;
			//fprintf(stdout, "2. Server sent %d bytes from file's data, offset is now : %d and remaining data = %d\n", sent_bytes, offset, remain_data);
		}

		close(fd);
	}
}

//Function which dispatches work to each of the worker threads based on the number
//specified in the command line param. Each thread in turn calls the FLANN
//create and search function
/*
 * 	Function: dispatch_to_worker_threads
	Parameters: none
	Purpose: Function which dispatches work to each of the worker threads based on the number specified in the command line param. It divides the database image set equally among all the available threads for effective load balancing. 
			 Each thread in turn calls the FLANN create and search function.
			 After all the threads return with their matched image info, it obtains the consolidated top 4 image matches. 
	Returns: void
 */
void dispatch_to_worker_threads(int num_data_images, int thread_count, thdata  *th_data, im_data *imageInfo) {
	int num = num_data_images / thread_count;
	pthread_t *threads;

	if ((threads = (pthread_t*) malloc(sizeof(pthread_t) * thread_count))
			== NULL) {
		cout << "Pthread alloc error\n";
		return;
	}

	for (int i = 0; i < thread_count; i++) {
		th_data[i].image_count = num;
	}

	//assign tasks to threads
	int i = 0;
	if (num * thread_count != num_data_images) {
		int temp_count = num_data_images - (num * thread_count);
		while (temp_count != 0 && i < thread_count) {
			th_data[i].image_count++;
			i++;
			temp_count--;
		}
	}

	th_data[0].from = 0;
	th_data[0].to = th_data[0].image_count;
	int sum = th_data[0].image_count;

	th_data[0].from_value = 0;
	th_data[0].to_value = TOP_FOUR;
	int temp_sum = TOP_FOUR;

	for (int i = 1; i < thread_count; i++) {
		th_data[i].from = th_data[i - 1].to;
		sum += th_data[i].image_count;
		th_data[i].to = sum;

		th_data[i].from_value = th_data[i - 1].to;
		temp_sum += TOP_FOUR;
		th_data[i].to_value = temp_sum;
	}

	for (int t = 0; t < thread_count; t++) {
		th_data[t].threadcount = t;
		int rc = pthread_create(&threads[t], NULL, call_image_matching_module,
				(void *) &th_data[t]);

		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		} else {
			cout << "Thread " << t << " Created Successfully\n";
		}
      	}

	for (int t = 0; t < thread_count; t++) {
		pthread_join(threads[t], NULL);
	}

	size_t structs_len = sizeof(imageInfo) / sizeof(im_data);

	qsort(imageInfo, structs_len, sizeof(im_data), int_cmp);
        
        int image_matches = 0;
        for (int i = 0; i < TOP_FOUR; i++) {
                if (imageInfo[i].match_found > 0)
                        image_matches++;
        }

        //int converted_number = htonl(image_matches);
        if (image_matches == 0) {
                printf("\nFound no matches in server database from phi\n");
      //          write(sock, &converted_number, sizeof(converted_number));
              //  return;
        }

        else {
                printf("\nSending top %d phi matches in server database\n", image_matches);
        //        write(sock, &converted_number, sizeof(converted_number));
        }

   	printf("Thread Count: %d\n",thread_count); 
   	printf("Total Create: %d\n",t_create); 
   	printf("Total Search: %d\n",t_search); 
	printf("Create Average time: %d\n", t_create / thread_count);
	printf("Search Average time: %d\n", t_search / thread_count);

}

//Function which does all the query image post processing once it recieves
//the key data from the client
/*
 * 	Function: dostuff
	Parameters: int sock - socket id
	Purpose: This function receives all the query image data one by one from the server. Then calls dispatch_to_worker_threads for doing the actual image search. Finally calls send_results_to_client to send the results back to the client. 
	Returns: void
 */

void dostuff(int sock) 
{
	int n;

	bzero(buffer, BUFFER);

	int received_int = 0;

	int return_status = read(sock, &received_int, sizeof(received_int));
	int converted_int = ntohl(received_int);

	if (return_status > 0) {
		//fprintf(stdout, "Received int = %d\n", converted_int);
	} else {
		printf("Error receiving count\n");
	}

	int *sizes_images;
	sizes_images = (int*) malloc(converted_int * sizeof(int));

	for (int s = 0; s < converted_int; s++) {
		bzero(buffer, BUFFERSIZE);
		recv(sock, buffer, BUFFERSIZE, 0);
		sizes_images[s] = atoi(buffer);
		//fprintf(stdout, "\nFile size recieved : %d\n", sizes_images[s]);
	}

	for (int k = 0; k < converted_int; k++) {
		client_packet.total_query_images = converted_int;
		client_packet.img_id = k;
		char buffer[BUFFERSIZE];
		int file_size;
		FILE *received_file;
		ssize_t len;
		int remain_data = 0;
		char name[MAX_FILE_LEN];

		remain_data = sizes_images[k];
		sprintf(query_key_names[k], "%s%d%s", db_path_query, k, sift4);
		sprintf(name, "%s%d%s", db_path_query, k, sift4);
		received_file = fopen(query_key_names[k], "w");
		if (received_file == NULL)
			cout << "error\n";

		while (((len = recv(sock, buffer, BUFFERSIZE, 0)) > 0) && (remain_data
				> 0)) {
			fwrite(buffer, sizeof(char), len, received_file);

			remain_data -= len;
			if (remain_data < BUFFERSIZE) {
				recv(sock, buffer, remain_data, 0);
				fwrite(buffer, sizeof(char), remain_data, received_file);
				break;
			}
			//fprintf(stdout, "Receive %d bytes and we hope :- %d bytes\n", len, remain_data);

		}
		fclose(received_file);

		n = write(sock, "I got your message", MAX_MSG);
		if (n < 0)
			error("ERROR writing to socket");

		process_query_data();

		cout << "Num query images: " << num_query_images << endl;
		cout << "Num data images: " << num_data_images << endl;
		dispatch_to_worker_threads(num_data_images, thread_count, th_data, imageInfo);
		send_results_to_client(sock);
	}

}

//Entry point for the server program
/*
 * 	Function: main
	Parameters: int argc - command line argument count
				char *argv[] - command line argument vector
	Purpose: The main entry point for the server program. Waits for a client connection. Once found, it receives all the query image info. Processes and macthes it against the server database. And sends the results back to the client. Continues to wait for another client to connect.
	Returns: 0
 */
int main(int argc, char *argv[]) {
	int sockfd, newsockfd, portno;
	socklen_t clilen;

	total_cpus_present = sysconf(_SC_NPROCESSORS_CONF);
	flann_threads = atoi(argv[7]);
	//cout << "Flann Threads: " << flann_threads << endl;

	signal(SIGCHLD, SIG_IGN);
	struct sockaddr_in serv_addr, cli_addr;
	int n;

	if (argc < ARGCOUNT) {
		fprintf(stderr, "ERROR, cmd line should be port, original img path,"
			"run_sift,  thread count, db_path for query, db_path for database\n");
		exit(1);
	}

	thread_count = atoi(argv[4]);

	strcpy(db_path_query, argv[6]);
	strcpy(db_path_data, argv[5]);
	KDT = atoi(argv[7]);
	FT = atoi(argv[8]);
	sockfd = socket(AF_INET, SOCK_STREAM, 0);
	if (sockfd < 0)
		error("ERROR opening socket");

	bzero((char *) &serv_addr, sizeof(serv_addr));
	portno = atoi(argv[1]);
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = INADDR_ANY;
	serv_addr.sin_port = htons(portno);

	if (bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)
		error("ERROR on binding");

	int z;
	int so_reuseaddr = 1;
	int pid;

	z = setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &so_reuseaddr,
			sizeof so_reuseaddr);

	db_init(argv[2], atoi(argv[3]));
	cout << "DB init over\n";

	listen(sockfd, LISTENPORT);

	clilen = sizeof(cli_addr);
	cout<<"End of listen done";
	while (1) {
		newsockfd = accept(sockfd, (struct sockaddr *) &cli_addr, &clilen);
		if (newsockfd < 0)
			error("ERROR on accept");
		t_sum = 0;
		t_create = 0;
		t_search = 0;
		pid = fork();
		if (pid < 0)
			error("ERROR on fork");
		if (pid == 0) {
			close(sockfd);
			dostuff(newsockfd);
			exit(0);
		}

		else
			close(newsockfd);
	} /* end of while */
	close(sockfd);

	//return 0;
}
