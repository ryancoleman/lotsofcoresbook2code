/*! \file
 *  \brief A Map Object that works lazily from Disk
 */

#include "qdp_map_obj_disk.h"

#include <errno.h>
#include <sys/stat.h>

namespace QDP 
{ 
  namespace MapObjDiskEnv 
  {
    // Anonymous namespace
    namespace {
      const std::string file_magic="XXXXQDPLazyDiskMapObjFileXXXX";
    };

    // Magic string at start of file.
    std::string getFileMagic() {return file_magic;}

    // Get meta-data
    std::string getMetaData(const std::string& filename)
    {
      BinaryFileReader reader;
      reader.open(filename);

      std::string read_magic;
      reader.readDesc(read_magic);
    
      // Check magic
      if (read_magic != MapObjDiskEnv::file_magic) { 
	QDPIO::cerr << "Magic String Wrong: Expected: " << MapObjDiskEnv::file_magic << " but read: " << read_magic << std::endl;
	QDP_abort(1);
      }

      MapObjDiskEnv::file_version_t read_version;
      read(reader, read_version);
      
      std::string user_data;
      readDesc(reader, user_data);
    
      reader.close();
      return user_data;
    }


    // Check for a file
    bool checkForNewFile(const std::string& file, std::ios_base::openmode mode)
    {
      bool new_file = false;

      if (Layout::primaryNode()) 
      {
	struct stat statbuf;
	if ((mode & std::ios_base::trunc) || (stat(file.c_str(), &statbuf) != 0 && (errno == ENOENT)))
	{
	  if (errno == ENOENT)
	    errno = 0;	/* In case someone looks at errno. */
	  new_file = true;
	}
      }

      QDPInternal::broadcast(new_file);

      return new_file;
    }
  }
    
}
