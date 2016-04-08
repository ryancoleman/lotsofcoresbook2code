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
// A simple abstraction of a node. Consists of a set of end points through which it communicates 
// with its "peers". The expected interfaces are around connect and accept for client and 
// server respectively.
//

#ifndef __SCIFNODE_H__
#define __SCIFNODE_H__

#include <scif.h>
#include <vector>

#define MAX_NAME_SIZE 32
#define SCIF_NODE_MAX_MSG_SIZE 32

enum ScifDeviceType 
{
	SCIF_DEVICE_HOST,
	SCIF_DEVICE_MIC
};

class ScifNode
{
	public:
		bool sync(scif_epd_t epd, const char *token) const;

		scif_epd_t get_epd(size_t index) const;
		scif_epd_t pop_epd() const;		

		std::vector<scif_epd_t> get_epd() const;
				
		inline unsigned int local_port() const { return l_port;} 
		inline void local_port(unsigned int port) { l_port = port; }

	protected:
		ScifNode(const char *node_name, ScifDeviceType type);
		virtual ~ScifNode();

		void barrier(unsigned int index) const;

	 	char name[MAX_NAME_SIZE];
		char partner[MAX_NAME_SIZE];
	
		ScifDeviceType deviceType;
		ScifDeviceType partnerType;
	
		// base local port
		unsigned int l_port;
		std::vector<scif_epd_t> epd;		

		// index for next available epd
		static unsigned int pop_epd_index;

	private:
		ScifNode(const ScifNode&);
		ScifNode& operator=(const ScifNode&);
};

//
// SCIF Client Node. Actively connects to a server.
//
class ScifClient : public ScifNode
{
	public:
		ScifClient(const char *node_name, ScifDeviceType type);
		~ScifClient();

		inline unsigned int server_port() const { return svr_port; };
		void connect(const char *server_name, unsigned int server_port, unsigned int n_conn = 1);	

	private:
		// Base server port
		unsigned int svr_port; 

		// retry count for connection establishment.
		unsigned int num_try;
};

//
// SCIF Server Node. Accepts incoming connections.
//
class ScifServer : public ScifNode
{
	public:
		ScifServer(const char *node_name, ScifDeviceType device);
		~ScifServer();

		void accept(const char *client_name, unsigned int local_port, unsigned int n_conn = 1);

	private:
};

#endif
