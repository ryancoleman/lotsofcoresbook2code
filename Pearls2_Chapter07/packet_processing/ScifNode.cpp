/*
 *
 * This file is provided under a dual BSD/GPLv2 license.  When using or
 * redistributing this file, you may do so under either license.
 *
 * GPL LICENSE SUMMARY
 *
 * Copyright(c) 2015 Intel Corporation.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of version 2 of the GNU General Public License as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * BSD LICENSE
 *
 * Copyright(c) 2015 Intel Corporation.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *  - Neither the name of Intel Corporation nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
//
// Abstraction for a SCIF Node
//

#include "ScifNode.h"
#include <cstring>
#include <iostream>

using namespace std;

unsigned int ScifNode::pop_epd_index = 0;

//-----------------------------------------------------------------------
//                                 SCIF Node 
//-----------------------------------------------------------------------
ScifNode::ScifNode(const char *node_name, ScifDeviceType type)
: deviceType(type), l_port(0)
{
	memset(name, 0, MAX_NAME_SIZE);
	strcpy(name, node_name);

	memset(partner, 0, MAX_NAME_SIZE);

	partnerType = (ScifDeviceType)(1 - (int)type);	
}

ScifNode::~ScifNode()
{
	for (size_t i = 0; i < epd.size(); ++i) {
		if (scif_close(epd[i]) != 0) {
			cerr << "ScifNode::~ScifNode close failed" << endl;
      		}
   	}
}

//
// A "sync" is simply send a message and receive one from the remote side, then compare the two 
// The send/recv happens over endpoint "epd". The message sent is "token".
//
bool ScifNode::sync(scif_epd_t epd, const char *token) const
{
	bool ret = true;
	int err;
	char message[SCIF_NODE_MAX_MSG_SIZE];

	memset(message, 0, SCIF_NODE_MAX_MSG_SIZE);

   	if ((err = scif_send(epd, (char*)token, SCIF_NODE_MAX_MSG_SIZE, 1)) <= 0) {
   		ret = false;
   	}

	if ((err = scif_recv(epd, message, SCIF_NODE_MAX_MSG_SIZE, 1)) <= 0) {
		ret = false;
	}

	if (strcmp(token, message) != 0) {
		cerr << "ScifNode::sync token: " << token << " message: " << message <<  endl;
		ret = false;
	}

	if (ret) {
		cerr << "ScifNode::sync sync with token " << token << "complete" << endl;
	} else {
		cerr << "ScifNode::sync sync with token " << token << "failed" << endl;
	}

	return ret;
}

//
// A "barrier" is simply a send/recv (like "sync" above), except it is with end point "index"
// Used internally as a true barrier across endpoints.
//
void ScifNode::barrier(unsigned int index) const
{
	int err;
	char message[SCIF_NODE_MAX_MSG_SIZE];
	memset(message, 0, SCIF_NODE_MAX_MSG_SIZE);

	if ((err = scif_send(epd[index], (char *)partner, SCIF_NODE_MAX_MSG_SIZE, 1)) <= 0) {
		cerr << "ScifNode::barrier send failed" << endl;
		// error handling omitted for simplicity 
	}

	if ((err = scif_recv(epd[index], message, SCIF_NODE_MAX_MSG_SIZE, 1)) <= 0) {
		cerr << "ScifNode::barrier recv failed" << endl;
		// error handling omitted for simplicity 
	}

	if (strcmp(name, message) != 0) {
		cerr << "ScifNode::barrier name: " << name << " message: " << message <<  endl;
		// error handling omitted for simplicity 
	}
}

//
// Return an the SCIF endpoint for "index" (stored in the vector)
//
scif_epd_t ScifNode::get_epd(size_t index) const
{
	if (epd.size() == 0) {
		cerr << "ScifNode::get_epd empty" << endl;
		// error handling omitted for simplicity 
	}

 	if (index > epd.size() - 1)
	{
		cerr << "ScifNode::get_epd index is greater list size" << endl;
		// error handling omitted for simplicity 
	}

	return epd[index];
}

scif_epd_t ScifNode::pop_epd() const 
{
	if (epd.size() == 0) {
		cerr << "ScifNode::pop_epd: empty" << endl;
		// error handling omitted for simplicity 
	}	

	if (pop_epd_index > epd.size() - 1) {
		cerr << "ScifNode::pop_epd nothing to pop" << endl;
		// error handling omitted for simplicity 
	}

	return epd[pop_epd_index++];
}

vector<scif_epd_t> ScifNode::get_epd() const
{
   if (epd.size() == 0) {
	cerr << "ScifNode::get_epd empty" << endl;
	// error handling omitted for simplicity 
   }

   return epd;
}

//---------------------------------------------------------------------
//                                 SCIF Client
//---------------------------------------------------------------------
ScifClient::ScifClient(const char *node_name, ScifDeviceType type)
: ScifNode(node_name, type), svr_port(0), num_try(200)
{
}

ScifClient::~ScifClient()
{
}

// 
// Create "n_conn" number of connections to a listening server.
//
void ScifClient::connect(const char *server_name, unsigned int server_port, unsigned int n_conn)
{
	strcpy(partner, server_name);
	svr_port = server_port;	

	int conn_port;
	scif_epd_t epd_t;
	
	if (l_port == 0) {
		l_port = server_port;
	}

	// Start by opening n_conn endpoints and binding them to a local port (l_port + i)
	for (unsigned int i = 0; i < n_conn; ++i) {
		// open end point
		if ((epd_t = scif_open()) < 0) {
			cerr << "ScifClient::connect open failed" << endl;
			// error handling omitted for simplicity 
		}
      
		// bind end point to available port, generated dynamically
		if ((conn_port = scif_bind(epd_t, l_port + i)) < 0) {
			cerr << "ScifClient::connect bind failed" << endl;
			// error handling omitted for simplicity 
		}

		cerr << "ScifClient::connect successful bind to local port" << conn_port 
								<< " epd: " << epd_t << endl;
		epd.push_back(epd_t);
	}

	// initiate a connection to remote node, when successful returns the
	// peer portID. Re-tries for 20 seconds and exits with error message

	scif_portID portID;

	portID.node = (int)partnerType;
	portID.port = svr_port;

	unsigned int tries;
	for (size_t i = 0; i < epd.size(); ++i) {
		tries = num_try;
		while (tries) {
         		if (scif_connect(epd[i], &portID) < 0) {
				cerr << "ScifClient::connect connection to node " 
						<< portID.node << " failed, trial " << tries << endl;
				tries--;
				sleep(1);
				continue;
			}
         		break;
      		}

		if (tries == 0) {
			cerr << "ScifClient::connect connection failed " << endl;
			// error handling omitted for simplicity 
		}

		barrier(i);
	}

	cout << "all connections successful" << endl;
}


// -----------------------------------------------------------------------
//                                 SCIF Server
// -----------------------------------------------------------------------
ScifServer::ScifServer(const char *node_name, ScifDeviceType type)
: ScifNode(node_name, type)
{
}

ScifServer::~ScifServer()
{
}

// 
// Listen on a port and accept connections as they come in.
//

void ScifServer::accept(const char *client_name, unsigned int local_port, unsigned int n_conn)
{
	strcpy(partner, client_name);
	l_port = local_port;
	
	int conn_port;
	scif_epd_t epd_t;

	// Begin by opening a SCIF endpoint.
	if ((epd_t = scif_open()) < 0) {
		cerr << "ScifServer::accept open failed for server" << endl;
		// error handling omitted for simplicity 
	}

	// bind end pt to specified port 
	if ((conn_port = scif_bind(epd_t, l_port)) < 0) {
		cerr << "ScifServer::accept bind failed for server" << endl;
		// error handling omitted for simplicity 
	}

	cerr << "bind success to port " << conn_port << endl;

	// marks andpoint as listening end pt and queues up a maximum of BACKLOG
	// number of incoming connection requests
	if (scif_listen(epd_t, n_conn) != 0) {
		cerr << "ScifServer::accept listen failed" << endl;
		// error handling omitted for simplicity 
	}

	scif_portID portID;
	scif_epd_t newepd;

	// accepts a conn request by creating a new end pt that connects to peer 
	for (unsigned int i = 0; i < n_conn; ++i) {
		if (scif_accept(epd_t, &portID, &newepd, SCIF_ACCEPT_SYNC) != 0) {
			cerr << "ScifServer::accept accept failed" << endl;
			// error handling omitted for simplicity 
		}
		
		cerr << "ScifServer::accept accepted connection to local port " 
					<< portID.port << " epd: " << newepd <<  endl;

		epd.push_back(newepd);		
		barrier(i);
	}

	cout << "all connections successful" << endl;

	if (scif_close(epd_t) != 0) {
		cerr << "ScifServer::accept scif_close failed" << endl;
		// error handling omitted for simplicity 
	}
}
