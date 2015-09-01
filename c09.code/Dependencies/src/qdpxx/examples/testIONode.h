#ifndef TEST_IO_NODE_H
#define TEST_IO_NODE_H

// These tests are disabled -- pending what we do with the 
// QIO

#if 0 

#include <qdp.h>
#include <qio.h>
#include "unittest.h"

using namespace Assertions;


//! Single File IONode TestCases
class TestSingleFileIONode : public TestCase {
 public:
  void run(void) {
    // IO Node is Node 0 -- the only node
    assertEquals<int>(0, SingleFileIONode::IONode(0));
    // Master IO Node is Node 0 -- the only node
    assertEquals<int>(0, SingleFileIONode::masterIONode());
  }
};

class TestMultiFileIONode : public TestCase { 
 public:
  void run() {
    // Multifile -- every file writes
    // Check that a node's IO node is itself
    assertEquals<int>(Layout::nodeNumber(), MultiFileIONode::IONode(Layout::nodeNumber()));
    // The master IO Node is whatever QDP Claims it to be
    assertEquals<int>(DML_master_io_node(), MultiFileIONode::masterIONode());
  }
};

class TestPartFileIONode1 : public TestCase { 
 public:
  void run() {
    // Choke Hypercube down to its "origin"
    // In this implementation everyIO node has to be divisible by 2 for
    // every node
    int my_io_node = PartFileIONode::IONode(Layout::nodeNumber());
    multi1d<int> ioNodeCoords = Layout::getLogicalCoordFrom(my_io_node);
    for(int i=0; i < ioNodeCoords.size(); i++) { 
      assertEquals<bool>(true, (ioNodeCoords[i]%2 == 0));
    }
  }
};

class TestPartFileIONode2 : public TestCase { 
 public:
  void run() {
    // Choke Hypercube down to its "origin"
    // In this implementation everyIO node has to be divisible by 2 for
    // every node
    int my_io_node = PartFileIONode::IONode(Layout::nodeNumber());
    multi1d<int> ioNodeCoords = Layout::getLogicalCoordFrom(my_io_node);
    multi1d<int> myCoords = Layout::getLogicalCoordFrom(Layout::nodeNumber());
    for(int i=0; i < ioNodeCoords.size(); i++) { 
      assertEquals<bool>(true, ( (myCoords[i] - ioNodeCoords[i]) >= 0));
      assertEquals<bool>(true, ( (myCoords[i] - ioNodeCoords[i]) < 2));
    }
  }
};

class TestPartFileIONode3 : public TestCase { 
 public:
  void run() {
    assertEquals<int>(DML_master_io_node(),  PartFileIONode::masterIONode());
  }
};

#endif 
#endif
