Installation is a two step process, one is to download and build FLANN
library and the next is to build the actual visual search
application. If if you want to run the code on Xeon Phi, you have to
compile the FLANN library with -mmic flag. Instructions on building
FLANN library on both Xeon and Xeon Phi is given below,

Prerequisites:
Have the FLANN Library downloaded from 
http://www.cs.ubc.ca/research/flann/uploads/FLANN/flann-1.8.4-src.zip

These steps have been automated with the scripts:
#get and build flann
sh 1_GET_BUILD_FLANN.sh
#build the client and server
sh 2_BUILD_SERVER_CLIENT.sh
#Run on Xeon 
sh 3_RUN_XEON.sh
#Run on MIC
sh 4_RUN_MIC.sh

Build the library and export the LD_LIBRARY_PATH as follows,

STEP 1:

#for explicit commands
To change the default compilation from gcc to icc
mkdir build
cd build
cmake ../ -DCMAKE_C_COMPILER="/path/to/icc" -DCMAKE_C_FLAGS=" -O3 -mmic(incase of mic compilation)" -DCMAKE_CXX_COMPILER="/path/to/icpc" -DCMAKE_CXX_FLAGS=" -O3 -mmic(incase of mic compilation)"
make
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/path_to_flann-version/build/lib/

Like wise for the library to run on xeon phi, build the library with -mmic
flags and export MIC_LD_LIBRARY_PATH in the same fashion.

STEP 2:
There are three make files to build seperately, Pick on among the two serverdir or serverdir_phi accordingly, the former runs the code on Xeon and 
the latter runs a hybrid version of Xeon and Xeon phi

1.serverdir - Which is for xeon run (make sure you have completed STEP 1 by exporting the FLANN library path compiled for Xeon)
2.serverdir_phi - which is for mix of xeon and xeon phi runs (make sure you have completed STEP 1 by exporting the FLANN mic library path compiled for Xeon phi)
3.clientdir - builds the client which connects to the server

*******Change the FLANN and PROJECT PATHS in the MAKEFILE**********

Run the client or server
//Open Make File and change FLANN_PATH to the respective flann directory
To run the server:
1.      Make
2.      ./server <port-no> <original loc of pgm files> <run sift[0/1]> <no of DB threads[1-16]> <db_path of key files (/mnt/pmfs/data, /mnt/pmbd/data, /disk2/data)> <no of KD trees> <no of FLANN threads> &
 a. Ex:  ./server 5304 $(PROJ_DIR)/alldb.txt 0 8 $(PROJ_DIR)/data/ $(PROJ_DIR)/input.txt 1 1 &

There is an example run script in the server directory. 
WARNING! The server executable is very sensitive to command-line input and
must have absolute paths plus a '/' at the end of 'data/'!

To run the client:
1. make
2. ./client localhost <portno> <inputfile db_path of key files>
2a. Ex ../clientdir/client localhost 5304 $(PROJ_DIR)/input.txt Client.cpp

//Make sure that the paths in all_db.txt and input.txt are pointing to the right path

High Level Flow of the Client
1. Creates a new socket communication endpoint 
2. Connects to a server which is listening for connections
3. Client gets the query image location from the input file provided on command line
4. Generates the SIFT key files using the SIFT algorithm
5. Send the query key data the server
		a. First sends the number of query images server should expect
		b. Then sends the size of the image key data server should expect
6. For each query image in the input file
		a. Send the query key data to the server
		b. Receives all the matching results information from the server and post processes it
		The function invoked to do this does the following:
				a. First receives the number of result images server will send it
				b. Then receives the size of the result image key data client should expect, and finally receives the key data from the server
		c. All the matches results are stored as key files on the client
9. Client closes the TCP socket connection

Makefile - makefile of the client program. Contains all the compile information.

1.	Function: main
	Parameters: int argc - command line argument count
				char *argv[] - command line argument vector
	Purpose: The main entry point for the client program. Creates a socket connection with the server and sends the query images key files one by one to the server.
			 Once all the results are obtained at the client, the socket connection is closed.
			 
			 Brief description: 
			 1. socket()
			 2. connect()
			 3. generate key file for all query images()
			 4. send query image count to server()
			 5. send query key points size to server()
			 6. for each query image
					send bytes buffer (key data) to server
					receive ACK()
					receive matching images from server()
				end for
			7. close socket connection
			
	Returns: 0
	
2.	Function: construct_client_featurefile
	Parameters: char *fpath[] - input query file having the locations of all the query images
	Purpose: Calls the many functions like init_client, get_all_query_image_locations and generate_sift_key_files to prepare to send the SIFT feature files to the server for image matching.
	Returns: void
	
3.	Function: init_client
	Parameters: none
	Purpose: Mallocs memory for certain variables used throughout the program
	Returns: void
	
4.	Function: get_all_query_image_locations
	Parameters: char *input_file_name - input query file having the locations of all the query images
	Purpose: Based on the input file name specified in the command line param, reads all the query images locations and stores them into arrays.
	Returns: void 
	
5.	Function: generate_sift_key_files
	Parameters: none
	Purpose: Runs the SIFT algorithm on all the query images and generates the .key files. 
	Returns: void

6.	Function: receive_picture_server
	Parameters: int sock - 
				int picid - 
	Purpose: It does two things:
			 First it receives the number of result matches that the server has found for the query image sent. (Top 4 matches for the query images are received from the server).
			 If no matches are found, it suitably prints that no matches were found on the server database.
			 Next, it receives the actual key file information through the recv() function from the server.
			 All the key files of the matching images received are stored at the client. 
	Returns: void

7.	Function: error
	Parameters: const char *msg - Error message to be printed
	Purpose: Prints error message
	Returns: void 	

Server.cpp

High Level Flow of the Server
1. Creates a new socket communication endpoint
2. Server does DB init. In DB init, This function gets the locations of all the database images and stores them into an array.  database images and stores them into an array. If runsift parameter is 1, then it runs the SIFT algorithm to kind the key points of all the database images. It stores them as sift key files. (This parameter needs to be 1 everytime the database images have changes or being loaded for the first time) It then calls generateKey_mmap to map the key files to process's memory from specified storage areas (PMFS, PMBD, SATA)
2. Listen for incoming client connections
3. Once it receives a client request, it accepts the connection
4. Then it calls a function to do post processing of all the query images sent one by one by the client.
	a. Receives the number of query images to process
	b. Size of the key files of the query images
	c. Query image key data
5. For the query image key data received, it calls process_query_data(), generates the structured key file and mmaps them into process's memory.
6. dispatch_to_threads()
	a. Worker threads are spwaned to divide the image database set equally among all threads
7. The threads are bound to specific core in bind_to_cpu()
8. call_flannC() is called for the actul image matching
9. call_flannC()
	for each image in the subset assigned to each thread
		Build kd tree index of database image
		Search query points in this kd tree
		Find no of matches
		Determine if this image is counted towards a match
		Obtain the top 4 matches
10. Obatins results from all the threads
11. Get top 4 matches from them
12. Send the 4 matching images key files back to the client 
13. Wait for another client to connect

Makefile - makefile of the server program. Contains all the compile information.

1.	Function: main
	Parameters: int argc - command line argument count
				char *argv[] - command line argument vector
	Purpose: The main entry point for the server program. Waits for a client connection. Once found, it receives all the query image info. Processes and macthes it against the server database. And sends the results back to the client. Continues to wait for another client to connect.
	Returns: 0
	
2.	Function: int_cmp
	Parameters: const void *a - pointer to one matched image structure
				const void *b - pointer to one matched image structure
	Purpose: Qsort structure comparison function. Function helps sort the macthed image structure img_data to obtain top 4 matches found. It returns negative if b>a and positive if a>b 
	Returns: int 
	
3.	Function: init_server
	Parameters: none
	Purpose: Mallocs memory for certain variables used throughout the program
	Returns: void
	
4.	Function: read_points_q
	Parameters: const char* filename - location of the structured key file of the database image
	Purpose: This function mmaps the structured binary key files from storage areas specified in the command line param (PMFS, PMBD, SATA) into the memory. It returns a pointer to the new mapping created in the virtual address space of the process.
	Returns: float* (pointer to float array)
	
5.	Function: generateKey_mmap
	Parameters: int count - no of database images
				int runsift - whether the SIFT needs to run or not (0/1) on the database images. (Supplied through the command line param)
	Purpose: This function generated the structured key files from the .key files produced from the SIFT algorithm. The SIFT algorithm produces .key files which has too many paramters. This function creates a file of only the 128 bit features. FLANN algorithm needs the input in this format.
	Returns: void

6.	Function: db_init
	Parameters: char *dbimgpath - database file having the locations of all the database images
				int runsift - whether the SIFT needs to run or not (0/1) on the database images. (Supplied through the command line param)
	Purpose: This function gets the locations of all the database images and stores them into an array. If runsift parameter is 1, then it runs the SIFT algorithm to kind the key points of all the database images. It stores them as sift key files. (This parameter needs to be 1 everytime the database images have changes or being loaded for the first time) It then calls generateKey_mmap to map the key files to process's memory from specified storage areas (PMFS, PMBD, SATA)
	Returns: void
	
7.	Function: process_query_data
	Parameters: none
	Purpose: Post processes the query image data obtained from the client. It also calls the read_points_q to mmap the query files too.
	Returns: void
	
8.	Function: write_results
	Parameters: const char* filename - Fielname to write the results to
				int *data - Return array from the FLANN search function
				int rows - No of keys in the client query image
				int cols - Dimension of the key data (128)
				float *dists - Distances obtained from the FLANN search function
				int image_no - Image id of the query image whose results are being written
	Purpose: This function is for debugging purposes only to write the image results locally to a file.
	Returns: void
	
9.	Function: call_flannC
	Parameters: int from - Since each thread is given subset of images. from indicates the start val of that image subset.
				int to - to indicates the to val of the image subset. 
				(from and to range from 0 to total_number_database_images)
				int from_value - Each thread maintains top 4 image matches it found. This from_value is an indicator of the array index it needs to store the top 4 matched image data.
				int to_value - end value of the array index it needs to store top 4 matched image data.
				(these integer are used as array indices of the img_data structure which stores all the information about the matches images)
	Purpose: This function does the actual image search.
			 It does the following things for each thread and its assigned subset of images:
				a. Build kd tree for each database image in its subset
				b. Search query key points in that kd tree
				c. Determines how many matches were obatined for each database image in its subset
				d. Obtains the top 4 matches once it finishes processing all the images in its subset.
	Returns: void
	
10.	Function: BindToCpu
	Parameters: int cpu_num - core number (which core the thread should be bound to)
	Purpose: Function which binds threads to specified core. Affinitizes thread to a core.
	Returns: int
	
11.	Function: call_image_matching_module
	Parameters: void *ptr - start_routine argument passed during pthread_create
	Purpose: This function first of all calls BindToCpu to affinitize each thread of a specific core. Then it calls call_flannC which does the FLANN building of kd tree and tree searching.
	Brief description: 
	1. call Bind to CPU for affinitization()
	2. call call_flannC() for image matching()
	Returns: NULL
	
12.	Function: send_results_to_client
	Parameters: int sock - socket id
	Purpose: It does two things:
			 First it sends the number of result matches that the server has found for eacgh query image sent. (Top 4 matches for the query images are sent from the server).
			 If no matches are found, it suitably prints that no matches were found on the server database and also sends that number to the client.
			 Next, it send the actual key file information to the client on the TCP socket.
			 The key files of all the top 4 matches are sent to the client.
	Returns: void
	
13.	Function: dispatch_to_worker_threads
	Parameters: none
	Purpose: Function which dispatches work to each of the worker threads based on the number specified in the command line param. It divides the database image set equally among all the available threads for effective load balancing. Each thread in turn calls the FLANN create and search function. After all the threads return with their matched image info, it obtains the consolidated top 4 image matches. 
	Brief description:
	1. Divide image subset equally among all threads 
	2. call thread module call_image_matching_module()
	3. wait for threads to complete
	4. Get top 4 matches from all the matches returned by threads
	5. Send results back to client
	Returns: void
	
14.	Function: dostuff
	Parameters: int sock - socket id
	Purpose: This function receives all the query image data one by one from the server. Then calls dispatch_to_worker_threads for doing the actual image search. Finally calls send_results_to_client to send the results back to the client. 
	Bried description:
	1. recv() query image data
	2. call function process_query_data()
		Generate structured SIFT key files and mmap key file data from PMFS, PMBD etc to process's memory
	3. dispatch worker threads()
	4. send results to client()	
	Returns: void
	
15.	Function: error
	Parameters: const char *msg - Error message to be printed
	Purpose: Prints error message
	Returns: none 	

Learnings and improvements
1. Image processing is a very challenging and interesting area to work
2. ANN library by Mount and Arya cant be used in a multithreaded env, it isnt thread safe
3. Put images in an actual database to avoid numerous file IO I do right now. This could significantly reduce response time

What would I do differently
1. Automate the SIFT creation to a script. I did create a script to do that. to_sift.sh. Would want to make more changes to that
2. Use a better image database. I have a lot of repeated images. 
3. Make changes to FLANN library to use benefits of PMFS (create the kd trees itself in PM)
4. Current code can be optimized lot more to prevent IOs etc. Did not get much time to do the same.
5. 
4. Load the images to a NoSQL DBMS like (MongoDB) etc to create an unstructured database of images and fetch and retrive from that
5. Do some preprocessing of images like face detection
6. Make ANN library by Mount and Arya thread safe. Since ANN library is a much more stable, commonly used library for image matching with SIFT features
7. Collect more EMON results. 
