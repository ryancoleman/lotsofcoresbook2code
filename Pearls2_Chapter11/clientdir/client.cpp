#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 
#include <iostream>
#include <sys/time.h>
#include <arpa/inet.h>
#include "defs.h"
#include <sys/sendfile.h>
#include <sys/stat.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include<fcntl.h>
using namespace std;

#define MAX_FILE_NAME_LENGTH 60
#define SMALL_BUFF_SIZE 18
#define TOP_FOUR 4
#define MAX_MSG 18
#define MAX_BUFFER 100
#define MAX_IMAGES 1000
#define MAX_FILE_LEN 100
#define MAX_SIFT_CMD_LEN 5000
#define PORT_NO 27015
#define DIM 128
#define MAX_KEYPOINTS		1000
#define LOWE_MATCH_RATIO	0.6
#define NNIDX_0				0
#define NNIDX_1				1
#define NN					2
#define MATCH_PERCENT		0.3
#define EPSILON 0.0
#define RANDOM_WALK_PERCENT	0.25
#define RANDOM_WALK	1
#define MAX_BUFFER_SIZE 1000
#define BUFFERSIZE 256

//global declarations
char db_path_query[MAX_FILE_LEN];
char *query_images[MAX_IMAGES];
char *query_key_names[MAX_IMAGES];
char input_file_name[MAX_FILE_NAME_LENGTH];

const char *sift1 = "./sift <";
const char *sift2 = " >";
const char *sift4 = ".key";

int num_query_images = 0;
//int offset = 0;
char temp_buffer[MAX_BUFFER_SIZE];

/*
 * 	Function: error
 Parameters: const char *msg - Error message to be printed
 Purpose: Prints error message
 Returns: void 
 */
void error(const char *msg) {
	perror(msg);
	exit(0);
}

/*
 * 	Function: init_client
    Parameters: none
    Purpose: Mallocs memory for certain variables used throughout the program
    Returns: void
 */
inline void init_client() {
	int i;

	for (i = 0; i < MAX_IMAGES; i++) {
		query_images[i] = (char *) malloc(MAX_FILE_LEN);
		query_key_names[i] = (char *) malloc(MAX_FILE_LEN);

	}
}

//Function to get all the query image locations from the input.txt file
/*
 *  Function: get_all_query_image_locations
    Parameters: char *input_file_name - input query file having the locations of all the query images
    Purpose: Based on the input file name specified in the command line param, reads all the query images locations and stores them into arrays.
    Returns: void
 */
inline void get_all_query_image_locations(char *input_file_name) {
	FILE *fp;

	fp = fopen(input_file_name, "r");
	int i = 0;
	char *c;
	num_query_images = 0;
	if (fp != NULL) {
		while (fgets(query_images[i], MAX_BUFFER, fp) != NULL) {
			c = strchr(query_images[i], '\n');
			if (c)
				*c = '\0';
			num_query_images++;
			i++;

		}

	}

	else {
		cout << "Error opening file, check location\n";
	}
	fclose(fp);
}

//Function to generate the query sift key files
/*
 * 	Function: generate_sift_key_files
    Parameters: none
    Purpose: Runs the SIFT algorithm on all the query images and generates the .key files. 
    Returns: void
 */
inline void generate_sift_key_files() {
	char siftcommand[MAX_SIFT_CMD_LEN];

	for (int i = 0; i < num_query_images; i++) {

		sprintf(siftcommand, "%s%s%s%s%d%s", sift1, query_images[i], sift2,
				db_path_query, i, sift4);
		sprintf(query_key_names[i], "%s%d%s", db_path_query, i, sift4);
		//offset += sprintf(temp_buffer + offset, "%s\t", query_key_names[i]);

		system(siftcommand);
	}

	//offset += sprintf(temp_buffer + offset, "\n");

}

/*
 * 	Function: construct_client_packet
 	Parameters: char *fpath[] - input query file having the locations of all the query images
    Purpose: Calls the many functions like init_client, get_all_query_image_locations and generate_sift_key_files to prepare to send the SIFT feature files to the server for image matching.
    Returns: void
 */
void construct_client_featurefile(char *fpath) {

	//initialize vairiables
	init_client();

	//open file pointer to actually get the query image locations in the database
	get_all_query_image_locations(fpath);

	//generate sift key files for the query_images
	generate_sift_key_files();

	//restructure the file into ANN format
	//generateKey(query_key_names, num_query_images);
}

//Receive all the key files one after the anotehr from the server
/*
 * 	Function: receive_picture_server
	Parameters: int sock - command line argument count
	 	 	 	int picid - command line argument vector
	Purpose: It does two things:
			 First it receives the number of result matches that the server has found for the query image sent. (Top 4 matches for the query images are received from the server).
	         If no matches are found, it suitably prints that no matches were found on the server database.
	         Next, it receives the actual key file information through the recv() function from the server.
	         All the key files of the matching images received are stored at the client. 
	Returns: void
 */
void receive_picture_server(int sock, int picid) {

	int received_int = 0;

	int return_status = read(sock, &received_int, sizeof(received_int));
	int converted_int = ntohl(received_int);
	if (return_status > 0) {
		//fprintf(stdout, "Received int = %d\n", converted_int);
	} else {
	//	printf("error receiving count\n");
	}

	if (converted_int == 0) {
		printf("Found no matches in server database printed from client\n");
		return;
	}

	int *sizes_images;
	sizes_images = (int*) malloc(converted_int * sizeof(int));

	/* Receiving file size */
	char buffer[BUFFERSIZE];
	int file_size;
	FILE *received_file;
	ssize_t len;
	int remain_data = 0;
	char name[MAX_FILE_LEN];

	for (int s = 0; s < converted_int; s++) {
		bzero(buffer, BUFFERSIZE);
		recv(sock, buffer, BUFFERSIZE, 0);
		sizes_images[s] = atoi(buffer);
		//fprintf(stdout, "\nFile size recieved : %d\n", sizes_images[s]);
	}

	struct stat st;
	char results_folder_name[MAX_FILE_LEN];

	sprintf(results_folder_name, "Result%d", picid);
	if (stat(results_folder_name, &st) == -1) {
		mkdir(results_folder_name, 777);
	}

	for (int s = 0; s < converted_int; s++) {
		bzero(buffer, BUFFERSIZE);
		sprintf(name, "ResultFile%d%d.key", picid, s);
		received_file = fopen(name, "w");
		if (received_file == NULL) {
			fprintf(stderr, "Failed to open file");
			exit(0);
		}

		remain_data = sizes_images[s];

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
		//fprintf(stdout, "\n");
	}

}

/*Entry point for the client program*/
/*
 * 	Function: main
	Parameters: int argc - command line argument count
	 	 	 	char *argv[] - command line argument vector
	Purpose: The main entry point for the client program. Creates a socket connection with the server and sends the query images key files one by one to the server.
	         Once all the results are obtained at the client, the socket connection is closed.
	Returns: 0
 */
int main(int argc, char *argv[]) {
	int sockfd, portno, n;
	struct sockaddr_in serv_addr;
	struct hostent *server;

	char buffer[BUFFERSIZE];
	if (argc < 4) {
		fprintf(stderr, "usage %s hostname port path_to_input_file\n", argv[0]);
		exit(0);
	}
	portno = atoi(argv[2]);

	sockfd = socket(AF_INET, SOCK_STREAM, 0);
	if (sockfd < 0)
		error("ERROR opening socket");

	bzero((char *) &serv_addr, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;

	if (inet_pton(AF_INET, argv[1], &serv_addr.sin_addr) <= 0) {
		printf("inet_pton error occured\n");
	}

	serv_addr.sin_port = htons(portno);
	if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)
		error("ERROR connecting");

	construct_client_featurefile(argv[3]);

	printf("Sending client packet... \n");

	int converted_query_count = htonl(num_query_images);
	write(sockfd, &converted_query_count, sizeof(converted_query_count));

	printf("There are %d query images\n", num_query_images);
	int *sizes_images;
	sizes_images = (int*) malloc(num_query_images * sizeof(int));

	for (int i = 0; i < num_query_images; i++) {
		int fd;
		char name_tmp[MAX_FILE_LEN];
		struct stat file_stat;
		sprintf(name_tmp, "%d", i);
		fd = open(query_key_names[i], O_RDONLY);
		if (fd == -1) {
			printf("Error\n");
			exit(0);
		}
		if (fstat(fd, &file_stat) < 0) {
			printf("Error\n");
			exit(0);
		}
		sizes_images[i] = file_stat.st_size;

		char file_size[BUFFERSIZE];
		ssize_t len;
		sprintf(file_size, "%d", file_stat.st_size);
		/* Sending file size */

		len = send(sockfd, file_size, sizeof(file_size), 0);
		if (len < 0) {
			fprintf(stderr, "Error on sending greetings");
			exit(0);
		}

		//fprintf(stdout, "Server sent %d bytes for the size\n", len);
		close(fd);
	}

	//Send all the query images one by one to the server
	
	for (int i = 0; i < num_query_images; i++) {

		timeval start;
		timeval end;
		gettimeofday(&start, NULL);

		int fd;
		fd = open(query_key_names[i], O_RDONLY);

		off_t offset;
		int remain_data;
		int sent_bytes = 0;
		offset = 0;
		remain_data = sizes_images[i];

		while (((sent_bytes = sendfile(sockfd, fd, &offset, BUFFERSIZE)) > 0)
				&& (remain_data > 0)) {
			//fprintf(stdout, "1. Server sent %d bytes from file's data, offset is now : %d and remaining data = %d\n", sent_bytes, offset, remain_data);
			remain_data -= sent_bytes;
			//fprintf(stdout, "2. Server sent %d bytes from file's data, offset is now : %d and remaining data = %d\n", sent_bytes, offset, remain_data);
		}

		close(fd);

		bzero(buffer, SMALL_BUFF_SIZE);
		n = read(sockfd, buffer, SMALL_BUFF_SIZE);
		if (n < 0)
			error("ERROR reading from socket");
		printf("%s\n", buffer);

		//Receive the key files from the server
		receive_picture_server(sockfd, i);

		gettimeofday(&end, NULL);

		double elapsedTime;

		elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0; // sec to ms
		elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0; // us to ms

		printf("\nResults for Image %d of Input file\n", i + 1);
		printf("Total time to get back to client: %g ms\n", elapsedTime);
		printf("Result key files are stored as Result/ResultFile* in current dir\n\n");
	}

	close(sockfd);

	return 0;
}
