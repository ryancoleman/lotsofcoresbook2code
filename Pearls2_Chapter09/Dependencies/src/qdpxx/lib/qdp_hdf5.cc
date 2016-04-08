//
// QDP data parallel interface
/*!
* @file
* @brief  HDF5 reading/writing of gauges and propagators
*/

#include "qdp.h"

#ifdef QDP_USE_HDF5

// optional hdf5 code goes here
namespace QDP {
	using std::string;

	// HDF5 classes
	//-------------------------------------------------------------------------------- 
	//--------------------------------------------------------------------------------
	// Base class
	//-------------------------------------------------------------------------------- 
	//-------------------------------------------------------------------------------- 
	//constructor, never called by user:

	HDF5::HDF5(const long int& stripesizee, const long int& maxalignn) :
		error_stack(H5E_DEFAULT), par_init(false), file_comm(H5P_DEFAULT),
		file_id(-1), current_group(-1),
		stripesize(stripesizee),
		maxalign(maxalignn)
	{
		//no profiling
		profile=false;
		isprefetched=false;

		//do not do cleanup after finishing:
		herr_t err=H5dont_atexit();

		//set error handler to custom error handler:
		H5E_auto2_t  oldfunc;
		void *client_data;
    
		//get old error handler
		H5Eget_auto(error_stack, &oldfunc, &client_data);
		
		//get new error handler
		//H5Eset_auto(error_stack,errorHandler,NULL);
		//H5Eset_auto(error_stack,NULL,NULL);
	}

	//checks whether an HDF5 file already exists:
	bool HDF5::check_exists(const std::string& filename)const{
		//on node 0, test if file exists:
		bool exists=false;
		if(Layout::nodeNumber()==0){
			std::ifstream input(filename.c_str(),std::ios_base::binary);
			if(input.good()){
				exists=true;
			}
			input.close();
			
			if(exists){
				htri_t ex=H5Fis_hdf5(filename.c_str());
				if(ex<=0) exists=false;
			}
		}
		QDPInternal::broadcast(exists);
		
		return exists;
	}

	//error handler:
	hid_t HDF5::errorHandler(hid_t errstack, void* unused){
		QDPIO::cout << "Some error occured, but we do not care at the moment!" << std::endl;
		return errstack;
	}

	//lookup routines:
	std::string HDF5::pwd()const{
		return getNameById(current_group);
	}

	std::string HDF5::parentDir()const{
		size_t pos;
		std::string dir=pwd();
		if((pos=dir.find_last_of("/"))==0) return std::string("/");
		return dir.substr(0,pos);
	}
  
	std::string HDF5::getNameById(hid_t id)const{
		ssize_t size=H5Iget_name(id,NULL,0);
		char* val=new char[size+1];
		H5Iget_name(id,val,size+1);
		std::string name(val);
		delete [] val;
		return name;
	}

	std::vector<std::string> HDF5::splitPathname(const std::string& name){
		//first split the string according to the file separators, such that the linking is correct:
		std::vector<std::string> dirlist;
		std::string nm=name;

		//get rid of possible leading file separator:
		if(nm.find_first_of("/")==0) nm=nm.substr(1);

		//no file separator left:
		if(nm.find_first_of("/")==std::string::npos){
			dirlist.push_back(name);
		}
		else{
			size_t pos=0;
			while((pos=nm.find_first_of("/"))!=std::string::npos){
				std::string tmpstr(nm.substr(0,pos));
				dirlist.push_back(tmpstr);
				nm=nm.substr(pos+1);
			}
			if(nm.compare("")!=0) dirlist.push_back(nm);
		}

		return dirlist;
	}

	//check if object exists:
	bool HDF5::objectExists(const ::std::string& name){
		return objectExists(file_id,name);
	}
	
	bool HDF5::objectExists(hid_t loc_id, const std::string& name){
		std::vector<std::string> dirlist=splitPathname(name);

		hid_t start_group;
		if(name.find_first_of("/")==0){
			start_group=file_id;
		}
		else start_group=loc_id;

		//iterate through the tree and check whether everything on the way exists:
		htri_t exists=H5Lexists(start_group,dirlist[0].c_str(),H5P_DEFAULT);
		if(exists!=1) return false;
		exists=H5Oexists_by_name(start_group,dirlist[0].c_str(),H5P_DEFAULT);
		if(exists!=1) return false;
		std::string tmpstring=dirlist[0];
		for(unsigned int i=1; i<dirlist.size(); i++){
			tmpstring+="/"+dirlist[i];
			exists=H5Lexists(start_group,tmpstring.c_str(),H5P_DEFAULT);
			if(exists!=1) return false;
			exists=H5Oexists_by_name(start_group,tmpstring.c_str(),H5P_DEFAULT);
			if(exists!=1) return false;
		}
		return true;
	}

	::std::string HDF5::objectType(const ::std::string& name){
		return objectType(file_id, name);
	}
    
	::std::string HDF5::objectType(hid_t loc_id, const ::std::string& name){
		std::string result="Null";
		if(objectExists(loc_id, name)){
			H5O_info_t objinfo;
			herr_t errhandle=H5Oget_info_by_name(loc_id,name.c_str(),&objinfo,H5P_DEFAULT);
			if(errhandle<0){
				::std::cerr << "HDF5::objectType: error, something went wrong with looking up " << name << "!" << std::endl;
			}
			else{
				switch(objinfo.type){
					case H5O_TYPE_GROUP:
					result="Group";
					break;
                        
					case H5O_TYPE_DATASET:
					result="Dataset";
					break;
                        
					case H5O_TYPE_NAMED_DATATYPE:
					result="NamedDatatype";
					break;
                        
					case H5O_TYPE_NTYPES:
					result="VariousTypes";
					break;
                        
					default:
					result="Unknown";
					break;
				}
			}
		}
		return result;
	}
	

	//navigation routines:
	void HDF5::pop(){
		//only do something if not at root-level:
		if(current_group!=file_id){
			hid_t tmp_group;

			//lookup what current and parent directory is:
			std::string cdir=pwd();
			std::string pdir=parentDir();
      
			//close group to prevent data loss and then go to new group:
			H5Gclose(current_group);
			tmp_group=H5Gopen(file_id,pdir.c_str(),H5P_DEFAULT);
      
			current_group=tmp_group;
		}
	}
  
	void HDF5::cd(const std::string& dirname){
		std::string tmpdirname(dirname);
		if(tmpdirname.compare("..")==0) pop();
		else if(tmpdirname.compare(pwd())!=0){
			hid_t tmp_group;
      
			//open HDF5 group relative to current_group:
			tmp_group=H5Gopen(current_group,tmpdirname.c_str(),H5P_DEFAULT);
      
			if(tmp_group<0){
				QDPIO::cerr << "HDF5::cd: error, the group " << tmpdirname << " does not exist!" << std::endl;
				return;
			}
			//only close current_group if tmp_group is not further below in the tree or completely outside:
			//get the names of the new and current dir:
			std::string ndir=getNameById(tmp_group);
			std::string cdir=pwd();
      
			//perform different tests:
			if(cdir.find(ndir)==0){
				//ndir is parent of cdir: iterate up from cdir and close all dirs on the way:
				while(current_group!=tmp_group){
					pop();
				}
			}
			else if(ndir.find(cdir)==string::npos){
				//ndir is not inside cdir tree, so close cdir tree:
				while(current_group!=file_id){
					pop();
				}       
			}
			current_group=tmp_group;
		}
	}

	//close-routine
	int HDF5::close(){
		if(file_id>0){
			herr_t err=1;
      
			unsigned int types=H5F_OBJ_DATASET | H5F_OBJ_GROUP |
				H5F_OBJ_DATATYPE | H5F_OBJ_ATTR;

			//get number of still open objects:
			ssize_t num_open = H5Fget_obj_count(file_id,types);
			if (num_open > 0) { 
				std::vector<hid_t> open_object_ids(num_open, 0); 
				H5Fget_obj_ids(file_id, types, num_open, &(open_object_ids.front()) ); 
				for(unsigned int i=0; i<num_open; i++){
					err=H5Oclose(open_object_ids[i]);
				}
			}
			err=H5Fclose(file_id);
		}

		//reset file_id:
		file_id=-1;

		return EXIT_SUCCESS;
	}

	//prefetch mapping for CB->lexicographical:
	int HDF5::prefetchLatticeCoordinates(){
		const int mynode = Layout::nodeNumber();
		
		//measure
		reordermap.resize(Layout::sitesOnNode());
		
		unsigned int run=0;
		for(int site=0; site < Layout::vol(); ++site){
			multi1d<int> coord = crtesn(site, Layout::lattSize());
			int node = Layout::nodeNumber(coord);

			if(node==mynode){
				reordermap[run]=Layout::linearSiteIndex(coord);
				run++;
			}
		}
		isprefetched=true;
		
		return EXIT_SUCCESS;
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//DATATYPE HELPERS                                                                                                                   
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//create complex
	hid_t HDF5::createComplexType(const unsigned int& float_size){
		hid_t complex_id = H5Tcreate(H5T_COMPOUND, 2*float_size);
		hid_t base_type;

		if(float_size==4) base_type=H5Tcopy(H5T_NATIVE_FLOAT);
		else base_type=H5Tcopy(H5T_NATIVE_DOUBLE);

		H5Tinsert(complex_id, "r", 0, base_type);
		H5Tinsert(complex_id, "i", float_size, base_type);
		H5Tclose(base_type);

		return complex_id;
	}

	//check complex:
	bool HDF5::checkComplexType(const hid_t& type_id, hid_t& base_type_id){
		//check if type_id is complex class:
		if(H5Tget_class(type_id)!=H5T_COMPOUND){
			return false;
		}

		if(H5Tget_nmembers(type_id)!=2) return false;
		base_type_id=H5Tget_member_type(type_id,0);
		htri_t err=H5Tequal(base_type_id,H5Tget_member_type(type_id,1));
		if(err<=0){
			base_type_id=-1;
			return false;
		}

		//check if data members are named r and i:
		bool error=false;
		char* name=H5Tget_member_name(type_id,0);
		if( strcmp(name,"r") != 0 ) error=true;
		free(name);
		name=H5Tget_member_name(type_id,1);
		if( strcmp(name,"i") != 0 ) error=true;
		free(name);
		if(error){
			base_type_id=-1;
			return false;
		}

		return true;
	}

	//create colormat                                                                            
	hid_t HDF5::createColorMatrixType(const hid_t& complex_id, const unsigned int& rank){
		if(rank>Nc || rank==0){
			HDF5_error_exit("HDF::createColorMatrixType: error, the rank should be between 1 and Nc");
		}
		hsize_t dims[]={rank,rank};
		hid_t colmat_id = H5Tarray_create(complex_id,2,dims);

		return colmat_id;
	}

	//check colormat:                                                                             
	bool HDF5::checkColorMatrixType(const hid_t& type_id, const unsigned int& rank, hid_t& base_type_id){
		///datatype should be array and have dimension 2:                                                                                                                                                                                                                                                                                                                     
		if( (H5Tget_class(type_id)!=H5T_ARRAY) || (H5Tget_array_ndims(type_id)!=2) ){
			return false;
		}

		//each dimension should have rank rank:
		hsize_t dims[2];
		H5Tget_array_dims(type_id,dims);
		if( (dims[0]!=dims[1]) || (dims[0]!=rank) ){
			return false;
		}

		//get base datatype and see whether it is the same for both dimensions:
		hid_t member_type=H5Tget_super(type_id);
		if(!checkComplexType(member_type,base_type_id)){
			return false;
		}
		H5Tclose(member_type);

		return true;
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//READING ATTRIBUTES                                                                                                                 
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//single datum
	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, short& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, unsigned short& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,false);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, int& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, unsigned int& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,false);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, unsigned long long& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,false);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, float& datum){
		rdAtt(obj_name,attr_name,datum,H5T_FLOAT,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, double& datum){
		rdAtt(obj_name,attr_name,datum,H5T_FLOAT,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, std::string& datum){
		std::string oname(obj_name), aname(attr_name);
		bool exists=objectExists(current_group,oname);
		if(!exists){
			HDF5_error_exit("HDF5::readAttribute: error, object "+oname+" you try to read attribute from does not exists!");
		}

		//do sanity checks and get datatype                                                                                                                                                                   
		hid_t ex=H5Aexists_by_name(current_group,oname.c_str(),aname.c_str(),H5P_DEFAULT);
		if(ex!=1){
			HDF5_error_exit("HDF5::readAttribute: error, the attribute "+aname+" you try to read does not exists!");
		}
		hid_t attr_id=H5Aopen_by_name(current_group,oname.c_str(),aname.c_str(),H5P_DEFAULT,H5P_DEFAULT);
		if(attr_id<0){
			HDF5_error_exit("HDF5::readAttribute: error, cannot open attribute "+aname+" attached to "+oname+"!");
		}
		hid_t type_id=H5Aget_type(attr_id);
		if(H5Tget_class(type_id)!=H5T_STRING){
			HDF5_error_exit("HDF5::readAttribute: error, datatype mismatch in attribute "+aname+"!");
		}
   
		//memory size
		hsize_t size=H5Tget_size(type_id);
		char* datumcpy=new char[size];
		hid_t nat_type_id=H5Tget_native_type(type_id,H5T_DIR_ASCEND);
		H5Aread(attr_id,nat_type_id,reinterpret_cast<void*>(datumcpy));
		datum=std::string(datumcpy);
		delete [] datumcpy;
		H5Aclose(attr_id);
		H5Tclose(nat_type_id);
		H5Tclose(type_id);
	}

	//array
	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<short>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<unsigned short>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,false);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<int>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<unsigned int>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,false);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<unsigned long long>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_INTEGER,false);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<float>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_FLOAT,true);
	}

	void HDF5::readAttribute(const std::string& obj_name, const std::string& attr_name, multi1d<double>& datum){
		rdAtt(obj_name,attr_name,datum,H5T_FLOAT,true);
	}
  
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//READING DATASETS                                                                                                                   
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//single datum
	void HDF5::read(const std::string& obj_name, short& datum){
		rd(obj_name,datum,H5T_INTEGER,true);
	}

	void HDF5::read(const std::string& obj_name, unsigned short& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, int& datum){
		rd(obj_name,datum,H5T_INTEGER,true);
	}

	void HDF5::read(const std::string& obj_name, unsigned int& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, unsigned long long& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, float& datum){
		rd(obj_name,datum,H5T_FLOAT,true);
	}

	void HDF5::read(const std::string& obj_name, double& datum){
		rd(obj_name,datum,H5T_FLOAT,true);
	}

	void HDF5::read(const std::string& dataname, std::string& datum){
		std::string dname(dataname);

		bool exists=objectExists(current_group,dname);
		if(!exists){
			HDF5_error_exit("HDF5::read: error, dataset does not exists!");
		}
		H5O_info_t objinfo;
		herr_t errhandle=H5Oget_info_by_name(current_group,dname.c_str(),&objinfo,H5P_DEFAULT);
		if(objinfo.type!=H5O_TYPE_DATASET){
			HDF5_error_exit("HDF5::read: error, "+dname+" exists but it is not a dataset!");
		}
		hid_t dset_id=H5Dopen(current_group,dname.c_str(),H5P_DEFAULT);
		if(dset_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open dataset!");
		}
		hid_t type_id=H5Dget_type(dset_id);
		if(H5Tget_class(type_id)!=H5T_STRING){
			HDF5_error_exit("HDF5::read: error, datatype mismatch in dataset "+dataname+"!");
		}

		//memory size
		hsize_t size=H5Tget_size(type_id);
		char* datumcpy=new char[size];
		hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		hid_t nat_type_id=H5Tget_native_type(type_id,H5T_DIR_ASCEND);
		H5Dread(dset_id,nat_type_id,H5S_ALL,H5S_ALL,plist_id,reinterpret_cast<void*>(datumcpy));
		datum=std::string(datumcpy);
		delete [] datumcpy;
		H5Pclose(plist_id);
		H5Dclose(dset_id);
		H5Tclose(nat_type_id);
		H5Tclose(type_id);
	}

	//array:
	//1D
	void HDF5::read(const std::string& obj_name, multi1d<short>& datum){
		rd(obj_name,datum,H5T_INTEGER,true);
	}

	void HDF5::read(const std::string& obj_name, multi1d<unsigned short>& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, multi1d<int>& datum){
		rd(obj_name,datum,H5T_INTEGER,true);
	}

	void HDF5::read(const std::string& obj_name, multi1d<unsigned int>& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, multi1d<unsigned long long>& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, multi1d<float>& datum){
		rd(obj_name,datum,H5T_FLOAT,true);
	}

	void HDF5::read(const std::string& obj_name, multi1d<double>& datum){
		rd(obj_name,datum,H5T_FLOAT,true);
	}
	//2D
	void HDF5::read(const std::string& obj_name, multi2d<short>& datum){
		rd(obj_name,datum,H5T_INTEGER,true);
	}

	void HDF5::read(const std::string& obj_name, multi2d<unsigned short>& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, multi2d<int>& datum){
		rd(obj_name,datum,H5T_INTEGER,true);
	}

	void HDF5::read(const std::string& obj_name, multi2d<unsigned int>& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, multi2d<unsigned long long>& datum){
		rd(obj_name,datum,H5T_INTEGER,false);
	}

	void HDF5::read(const std::string& obj_name, multi2d<float>& datum){
		rd(obj_name,datum,H5T_FLOAT,true);
	}

	void HDF5::read(const std::string& obj_name, multi2d<double>& datum){
		rd(obj_name,datum,H5T_FLOAT,true);
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//READING Compound types:                                                                                                            
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	void HDF5::readPrepare(const std::string& name, hid_t& type_id){
		//determine whether there is a dataset with the specified name:
		bool exists=objectExists(current_group,name);
		if(!exists){
			HDF5_error_exit("HDF5::read: error, dataset does not exists!");
		}
		H5O_info_t objinfo;
		herr_t errhandle=H5Oget_info_by_name(current_group,name.c_str(),&objinfo,H5P_DEFAULT);
		if(objinfo.type!=H5O_TYPE_DATASET){
			HDF5_error_exit("HDF5::read: error, "+name+" exists but it is not a dataset!");
		}

		//open dataset
		hid_t dset_id=H5Dopen(current_group,name.c_str(),H5P_DEFAULT);
		if(dset_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open dataset!");
		}
		type_id=H5Dget_type(dset_id);
		H5Dclose(dset_id);		
	}

	//complex types:
	//Single element:
	template<>void HDF5::read< PScalar< PScalar< RComplex<float> > > >(const std::string& dataname, ComplexF& datum){
		//get type:
		hid_t type_id;
		readPrepare(dataname,type_id);
		if(type_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open datatype!");
		}
		hid_t base_type_id;
		if(!checkComplexType(type_id,base_type_id)){
			HDF5_error_exit("HDF5::read: error, datatype mismatch!");
		}
		unsigned int float_size=H5Tget_size(base_type_id);
		if(float_size!=4){
			HDF5_error_exit("HDF5:read: error, datatype size mismatch!");
		}
		H5Tclose(type_id);
		H5Tclose(base_type_id);
    
		//perform the actual read:
		rd(dataname,datum,H5T_COMPOUND,true,false);
	}

	template<>void HDF5::read< PScalar< PScalar< RComplex<double> > > >(const std::string& dataname, ComplexD& datum){
		//get type:
		hid_t type_id;
		readPrepare(dataname,type_id);
		if(type_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open datatype!");
		}
		hid_t base_type_id;
		if(!checkComplexType(type_id,base_type_id)){
			HDF5_error_exit("HDF5::read: error, datatype mismatch!");
		}
		unsigned int float_size=H5Tget_size(base_type_id);
		if(float_size!=8){
			HDF5_error_exit("HDF5:read: error, datatype size mismatch!");
		}
		H5Tclose(type_id);
		H5Tclose(base_type_id);

		//perform the actual read:
		rd(dataname,datum,H5T_COMPOUND,true,false);
	}

	//array types
	template<>void HDF5::read<  PScalar< PScalar< RComplex<float> > > >(const std::string& dataname, multi1d<ComplexF>& datum){
		//get type:
		hid_t type_id;
		readPrepare(dataname,type_id);
		if(type_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open datatype!");
		}
		hid_t base_type_id;
		if(!checkComplexType(type_id,base_type_id)){
			HDF5_error_exit("HDF5::read: error, datatype mismatch!");
		}
		unsigned int float_size=H5Tget_size(base_type_id);
		if(float_size!=4){
			HDF5_error_exit("HDF5:read: error, datatype size mismatch!");
		}
		H5Tclose(type_id);
		H5Tclose(base_type_id);

		//perform the actual read:
		rd(dataname,datum,H5T_COMPOUND,true,false);
	}

	template<>void HDF5::read< PScalar< PScalar< RComplex<double> > > >(const std::string& dataname, multi1d<ComplexD>& datum){
		//get type:
		hid_t type_id;
		readPrepare(dataname,type_id);
		if(type_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open datatype!");
		}
		hid_t base_type_id;
		if(!checkComplexType(type_id,base_type_id)){
			HDF5_error_exit("HDF5::read: error, datatype mismatch!");
		}
		unsigned int float_size=H5Tget_size(base_type_id);
		if(float_size!=8){
			HDF5_error_exit("HDF5:read: error, datatype size mismatch!");
		}
		H5Tclose(type_id);
		H5Tclose(base_type_id);

		//perform the actual read:
		rd(dataname,datum,H5T_COMPOUND,true,false);
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//READING Lattice Types:                                                                                                            
	//***********************************************************************************************************************************
	//*********************************************************************************************************************************** 
	void HDF5::readPrepareLattice(const std::string& name, hid_t& type_id, multi1d<ullong>& sizes){
		//determine whether there is a lattice with the specified name:
		bool exists=objectExists(current_group,name);
		if(!exists){
			HDF5_error_exit("HDF5::read: error, dataset does not exists!");
		}
		H5O_info_t objinfo;
		herr_t errhandle=H5Oget_info_by_name(current_group,name.c_str(),&objinfo,H5P_DEFAULT);
		if(objinfo.type!=H5O_TYPE_DATASET){
			HDF5_error_exit("HDF5::read: error, "+name+" exists but it is not a dataset!");
		}

		//determine dimension of dataset:                                                               
		//open dataset and filespace:                                                                   
		hid_t dset_id=H5Dopen(current_group,name.c_str(),H5P_DEFAULT);
		if(dset_id<0){
			HDF5_error_exit("HDF5::read: error, cannot open dataset!");
		}
		type_id=H5Dget_type(dset_id);

		//get extent of dataset:                                                                        
		hid_t filespace=H5Dget_space(dset_id);
		int Ndims=H5Sget_simple_extent_ndims(filespace);
		sizes.resize(Ndims);
		hsize_t* inttoken=new hsize_t[Ndims];
		errhandle=H5Sget_simple_extent_dims(filespace,inttoken, NULL);
		for(unsigned int dd=0; dd<Ndims; dd++){
			sizes[dd]=inttoken[dd];
		}
		delete [] inttoken;

		H5Sclose(filespace);
		H5Dclose(dset_id);
		
		//prefetch for faster I/O:
		if(!isprefetched && sizes.size()>1) prefetchLatticeCoordinates();
	}

	void HDF5::readLattice(const std::string& name, const hid_t& type_id, const hid_t& base_type_id, const ullong& obj_size, const ullong& tot_size, REAL* buf) {
		// determine local sizes
		const int mynode = Layout::nodeNumber();
		const int nodeSites = Layout::sitesOnNode();

		// set up and create hyperslap:
		unsigned int dimensions;
		if(obj_size>1) dimensions=Nd+1;
		else dimensions=Nd;
		hsize_t* node_offset = new hsize_t[dimensions];
		hsize_t* offset = new hsize_t[dimensions];
		hsize_t* total_count = new hsize_t[dimensions];
		hsize_t* dim_size = new hsize_t[dimensions];

		// reorder such that x is fastest:
		for(unsigned int i = 0; i < Nd; i ++) {
			dim_size[i] = total_count[i] = Layout::subgridLattSize()[(Nd - 1) - i];
			offset[i] = node_offset[i] = Layout::nodeCoord()[(Nd - 1) - i] * total_count[i];
		} // for                                                                                                                                                                                                                                                                                                                                                              
		if(obj_size>1){
			dim_size[Nd] = total_count[Nd] = obj_size;
			offset[Nd] = node_offset[Nd] = 0;
		}

		//create property list for readin:
		hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		// open dataset and filespace and create hyperslab:
		hid_t dset_id = H5Dopen(current_group, name.c_str(), H5P_DEFAULT);
		hid_t filespace = H5Dget_space(dset_id);

		unsigned int float_size=H5Tget_size(base_type_id);
		if( (float_size!=4) && (float_size!=8) ){
			HDF5_error_exit("HDF5Reader::read: error, invalid base datatype while trying to read lattice!");
		}
		if(H5Tget_class(base_type_id)!=H5T_FLOAT){
			HDF5_error_exit("HDF5Reader::read: error, invalid base datatype while trying to read lattice!");
		}

		hid_t err;
		REAL32* buf32 = NULL;
		REAL64* buf64 = NULL;

		size_t two_gb = (size_t) 2 * 1024 * 1024 * 1024;
		size_t total_size = tot_size;
		unsigned int blocks = 1;
		while(total_size * float_size > two_gb) {
			dim_size[0] = dim_size[0] >> 1;
			total_size = total_size >> 1;
			blocks = blocks << 1;
		} // while  

		if(float_size == 4) {
			buf32 = new (std::nothrow) REAL32[total_size];
			if(buf32 == 0x0) {
				HDF5_error_exit("Unable to allocate buf\n");
			} // if                                                                                                                                                                                                                                                                                                                                                             
		} else {
			buf64 = new (std::nothrow) REAL64[total_size];
			if(buf64 == 0x0) {
				HDF5_error_exit("Unable to allocate buf\n");
			} // if                                                                                                                                                                                                                                                                                                                                                             
		} // if-else                                                                                                                                                                                                                                                                                                                                                          

		hsize_t rank = static_cast<hsize_t>(dimensions);
		hid_t memspace = H5Screate_simple(rank, dim_size, NULL);
		hid_t nat_type_id=H5Tget_native_type(type_id,H5T_DIR_ASCEND);

		// read:
		for(int i = 0; i < blocks; ++ i) {
			offset[0] = node_offset[0] + i * dim_size[0];
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, const_cast<const hsize_t*>(offset),
			NULL, const_cast<const hsize_t*>(dim_size), NULL);
			if(float_size == 4) {
				err = H5Dread(dset_id, nat_type_id, memspace, filespace, plist_id, static_cast<void*>(buf32));
				for(ullong j = 0; j < total_size; j ++)
					(buf + i * total_size)[j] = static_cast<REAL>(buf32[j]);
			} else {
				err = H5Dread(dset_id, nat_type_id, memspace, filespace, plist_id, static_cast<void*>(buf64));
				for(ullong j = 0; j < total_size; j ++)
					(buf + i * total_size)[j] = static_cast<REAL>(buf64[j]);
			} // if-else                                                                                                                                                                                                                                                                                                                                                        
		} // for                                                                                                                                                                                                                                                                                                                                                              

		if(buf32 != NULL) delete [] buf32;
		if(buf64 != NULL) delete [] buf64;
		delete [] dim_size;
		delete [] total_count;
		delete [] offset;
		delete [] node_offset;
		H5Pclose(plist_id);
		H5Sclose(memspace);
		H5Sclose(filespace);
		H5Tclose(nat_type_id);
		H5Dclose(dset_id);
	} // readLattice()  


	//read LatticeColorMatrix
	template<>void HDF5::read< PScalar< PColorMatrix< RComplex<REAL64>, 3> > >(const std::string& name, LatticeColorMatrixD3& field){
		StopWatch swatch_prepare, swatch_datatypes, swatch_reorder, swatch_read;
		
		//read dataset extents:
		if(profile) swatch_prepare.start();                                                                        
		multi1d<ullong> sizes;
		hid_t type_id;
		readPrepareLattice(name,type_id,sizes);
		if(profile) swatch_prepare.stop();

		//do some datatype sanity checks and get the basis precision:
		if(profile) swatch_datatypes.start();
		hid_t base_type_id;
		if( !checkColorMatrixType(type_id,Nc,base_type_id) ){
			HDF5_error_exit("HDF5::read: object is not a color matrix!");
		}

		unsigned int float_size=H5Tget_size(base_type_id);
		if(sizes.size()!=Nd){
			HDF5_error_exit("HDF5::read: error, wrong dimensionality!");
			return;
		}
		for(unsigned int dd=0; dd<Nd; dd++){
			if(sizes[Nd-dd-1]!=Layout::lattSize()[dd]){
				HDF5_error_exit("HDF5::read: mismatching lattice extents.");
			}
		}
		if(profile) swatch_datatypes.stop();

		//determine local sizes, allocate memory and read
		if(profile) swatch_read.start();
		const int mynode=Layout::nodeNumber();
		const int nodeSites = Layout::sitesOnNode();
		size_t obj_size=sizeof(ColorMatrixD3)/float_size;
		size_t tot_size = nodeSites*obj_size;
		REAL* buf = new(std::nothrow) REAL[tot_size];
		if( buf == 0x0 ) {
			HDF5_error_exit("Unable to allocate buf\n");
		}

		readLattice(name,type_id,base_type_id,1,tot_size,buf);
		H5Tclose(type_id);
		H5Tclose(base_type_id);
		if(profile) swatch_read.stop();

		//put lattice into u-field and reconstruct as well as reorder them on the fly:
		// Reconstruct the gauge field
		if(profile) swatch_reorder.start();
		/*#pragma omp parallel for firstprivate(nodeSites,obj_size,float_size) shared(buf,field)
		for(unsigned int run=0; run<nodeSites; run++){
			memcpy(&(field.elem(reordermap[run])),reinterpret_cast<char*>(buf+run*obj_size),float_size*obj_size);
		}*/
		CvtToLayout(field,reinterpret_cast<void*>(buf),nodeSites,float_size*obj_size);
		delete [] buf;
		if(profile) swatch_reorder.stop();
		
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics. Read:" << std::endl;
			QDPIO::cout << "\t preparing: " << swatch_prepare.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t datatype-handling: " << swatch_datatypes.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t reordering: " << swatch_reorder.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t read: " << swatch_read.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB read: " << Layout::vol()*sizeof(ColorMatrixD3)/1024/1024 << std::endl;
		}
	}

	//QDP Lattice IO:
	template<>void HDF5::read< PScalar< PColorMatrix< RComplex<REAL64>, 3> > >(const std::string& name, multi1d<LatticeColorMatrixD3>& field){
		StopWatch swatch_prepare, swatch_datatypes, swatch_reorder, swatch_read;
		
		//read dataset extents: 
		if(profile) swatch_prepare.start();                                                                                                                                                               
		multi1d<ullong> sizes;
		hid_t type_id;
		readPrepareLattice(name,type_id,sizes);
		if(profile) swatch_prepare.stop();  

		//do some datatype sanity checks and get the basis precision:
		if(profile) swatch_datatypes.start();
		hid_t base_type_id;
		if( !checkColorMatrixType(type_id,Nc,base_type_id) ){
			HDF5_error_exit("HDF5::read: object is not a color matrix!");
		}

		unsigned int float_size=H5Tget_size(base_type_id);
		if(sizes.size()!=(Nd+1)){
			HDF5_error_exit("HDF5::read: error, wrong dimensionality!");
			return;
		}
		for(unsigned int dd=0; dd<Nd; dd++){
			if(sizes[Nd-dd-1]!=Layout::lattSize()[dd]){
				HDF5_error_exit("HDF5::read: mismatching lattice extents.");
			}
		}
		field.resize(sizes[Nd]);
		if(profile) swatch_datatypes.stop();

		//determine local sizes, allocate memory and read
		if(profile) swatch_read.start();
		const int mynode=Layout::nodeNumber();
		const int nodeSites = Layout::sitesOnNode();
		size_t obj_size=sizeof(ColorMatrixD3)/float_size;
		size_t tot_size = nodeSites*sizes[Nd]*obj_size;
		REAL* buf = new(std::nothrow) REAL[tot_size];
		if( buf == 0x0 ) {
			HDF5_error_exit("Unable to allocate buf\n");
		}

		readLattice(name,type_id,base_type_id,sizes[Nd],tot_size,buf);
		H5Tclose(type_id);
		H5Tclose(base_type_id);
		if(profile) swatch_read.stop();

		//put lattice into u-field and reconstruct as well as reorder them on the fly:
		// Reconstruct the gauge field
		if(profile) swatch_reorder.start();
		unsigned int arr_size=sizes[Nd];
		/*#pragma omp parallel for firstprivate(nodeSites,arr_size,obj_size,float_size) shared(buf,field)
		for(unsigned int run=0; run<nodeSites; run++){
			for(unsigned int dd=0; dd<arr_size; dd++){
				memcpy(&(field[dd].elem(reordermap[run])),reinterpret_cast<char*>(buf+(dd+arr_size*run)*obj_size),float_size*obj_size);
			}
		}*/
		CvtToLayout(field,reinterpret_cast<void*>(buf),nodeSites,arr_size,float_size*obj_size);
		delete [] buf;
		if(profile) swatch_reorder.stop();
		
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics. Read:" << std::endl;
			QDPIO::cout << "\t preparing: " << swatch_prepare.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t datatype-handling: " << swatch_datatypes.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t reordering: " << swatch_reorder.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t read: " << swatch_read.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB read: " << Layout::vol()*field.size()*sizeof(ColorMatrixD3)/1024/1024 << std::endl;
		}
	}

	//special lattice IO:
	void HDF5::readFUEL(const std::string& name, multi1d<LatticeColorMatrixD3>& field){
		StopWatch swatch_complete;
		
		QDPIO::cout << "Reading FUEL config ..." << std::flush;

		//do some checks
		if(profile) swatch_complete.start();
		if(!objectExists(current_group,name)){
			HDF5_error_exit("HDF5::readFUEL: error, configuration "+name+"does not exist!");
		}
		//check if the object is a group and if it contains the datasets 0 to Nd-1 but not more:
		H5O_info_t objinfo;
		herr_t errhandle=H5Oget_info_by_name(current_group,name.c_str(),&objinfo,H5P_DEFAULT);
		if(objinfo.type!=H5O_TYPE_GROUP){
			HDF5_error_exit("HDF5::readFUEL: error, "+name+" exists but it is not a FUEL config!");
		}
		//check if datasets from 0 to Nd-1 exist:
		for(unsigned int dd=0; dd<Nd; dd++){
			std::stringstream stream;
			stream << name << "/" << dd;
			std::string dname=stream.str();
			if(!objectExists(current_group,dname)){
				HDF5_error_exit("HDF5::readFUEL: error, "+name+" exists but it is not a FUEL config! Dataset "+dname+" was not found!");
			}
		}
		//now, check if there is another entry and whether all these entries are LatticeColorMatrices:
		std::stringstream stream;
		stream << name << "/" << Nd;
		std::string dname=stream.str();
		if(objectExists(current_group,dname)){
			HDF5_error_exit("HDF5::readFUEL: error, "+name+" exists but it is not a FUEL config!");
		}

		//read LatticeColorMatrices:
		field.resize(Nd);
		for(unsigned int dd=0; dd<Nd; dd++){
			std::stringstream stream;
			stream << name << "/" << dd;
			std::string dname=stream.str();
			read(dname,field[dd]);
		}
		if(profile) swatch_complete.stop();
		QDPIO::cout << "done!" << std::endl;
		
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics for FUEL-read: " << std::endl;
			QDPIO::cout << "\t total: " << swatch_complete.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB read: " << Layout::vol()*field.size()*sizeof(ColorMatrixD3)/1024/1024 << std::endl;
		}
	}

	//-------------------------------------------------------------------------------- 
	//--------------------------------------------------------------------------------
	// Reader class
	//-------------------------------------------------------------------------------- 
	//-------------------------------------------------------------------------------- 
	//constructors:
	HDF5Reader::HDF5Reader() : HDF5() {};

	HDF5Reader::HDF5Reader(const long int& blocksize, const long int& maxalign) : HDF5(blocksize, maxalign) {};
  
	HDF5Reader::HDF5Reader(const std::string& filename): HDF5(){
		open(filename);
	};

	//destructors:
	HDF5Reader::~HDF5Reader(){
		close();
	};
  
	void HDF5Reader::open(const std::string& filename){
		//close actual file:
		close();
		
		//check if file exists first:
		bool exists=check_exists(filename);
		if(!exists){
			HDF5_error_exit("HDF5Reader::open: error, file does not exist!");
		}

		//create file access property control list:
		hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);

		//switch on LUSTRE optimizations:
		if(stripesize>0){
			//memory alignment:
			H5Pset_alignment(fapl_id,maxalign,stripesize);
		}

		//create MPI_IO accessor:
		MPI_Info info  = MPI_INFO_NULL;
		QMP_get_hidden_comm(QMP_comm_get_default(),reinterpret_cast<void**>(&mpicomm));
		H5Pset_fapl_mpio(fapl_id,*mpicomm, info);
   
		
		file_id=H5Fopen(filename.c_str(),H5F_ACC_RDONLY,fapl_id);
		H5Pclose(fapl_id);
      
		if(file_id<0){
			HDF5_error_exit("HDF5Reader::open: could not open file "+filename+" for reading. Some error occured!");
		}
      
		current_group=file_id;
	}


	//-------------------------------------------------------------------------------- 
	//--------------------------------------------------------------------------------
	// Writer class
	//-------------------------------------------------------------------------------- 
	//-------------------------------------------------------------------------------- 
	//constructors
	HDF5Writer::HDF5Writer() : HDF5() {};
  
	HDF5Writer::HDF5Writer(const long int& stripesizee, const long int& maxalignn) : HDF5(stripesizee,maxalignn) {};

	HDF5Writer::HDF5Writer(const std::string& filename, const HDF5Base::writemode& mode) : HDF5(){
		open(filename,mode);
	};

	//destructors
	HDF5Writer::~HDF5Writer(){
		close();
	};
  
	//member functions
	void HDF5Writer::open(const std::string& filename, const HDF5Base::writemode& mode){
		//close actual file if open:
		close();

		//try to create file. First create file access control list:
		hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
		hid_t fcpl_id = H5Pcreate(H5P_FILE_CREATE);

		//switch on LUSTRE optimizations:
		if(stripesize>0){ 
			//memory alignment:
			H5Pset_alignment(fapl_id,maxalign,stripesize);
			//binary-tree optimization:
			int btree_ik = ceil((stripesize - 4096) / 96);
			H5Pset_istore_k(fcpl_id,btree_ik);
		}

		//activate parallel IO:
		//this is MPI stuff and breaks the whole USQCD software design paradigm, however, at the
		//moment it is not possible to get MPI objects from QMP routines. This will be implemented
		//as soon as this hack works:  
		MPI_Info info  = MPI_INFO_NULL;
		QMP_get_hidden_comm(QMP_comm_get_default(),reinterpret_cast<void**>(&mpicomm));
		H5Pset_fapl_mpio(fapl_id,*mpicomm,info);

		//collective file creation/opening:
		//on node 0, test if file exists:
		bool exists=false;
		if(Layout::nodeNumber()==0){
			std::ifstream input(filename.c_str(),std::ios_base::binary);
			if(input.good()){
				exists=true;
			}
			input.close();
		}
		QDPInternal::broadcast(exists);

		if(exists){
			//file exists, check if it is an HDF5 file:
			htri_t ex=H5Fis_hdf5(filename.c_str());
			if(ex>0){
				//file is existing HDF5 file:
				if(mode&HDF5Base::trunc){
					//it exists and should be overwritten, so truncate it:
					file_id=H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,fcpl_id,fapl_id);
				}
				else{
					//it exists and should not be overwritten, so open it rw:
					file_id=H5Fopen(filename.c_str(),H5F_ACC_RDWR,fapl_id);
				}
			}
			else{
				//file exists but is not HDF5 file:
				if(mode&HDF5Base::trunc){
					//it should be overwritten, so truncate it:
					remove(filename.c_str());
					file_id=H5Fcreate(filename.c_str(),H5F_ACC_EXCL,fcpl_id,fapl_id);
				}
				else{
					HDF5_error_exit("HDF5Writer::open: error, file "+filename+" already exists and is not HDF5 file! Please use overwrite=true to create anew file!");
				}
			}
		}
		else{
			//it does not exists, so simply create a new file:
			file_id=H5Fcreate(filename.c_str(),H5F_ACC_EXCL,fcpl_id,fapl_id);
		}
		H5Pclose(fapl_id);
		H5Pclose(fcpl_id);

		if(file_id<0){
			HDF5_error_exit("HDF5Writer::open: could not open file "+filename+" for writing.");
		}

		current_group=file_id;
	}

	//create a new group inside current one w/o steping into it:
	void HDF5Writer::mkdir(const ::std::string& name){
		std::string cwd=pwd();
		push(name);
		cd(cwd);
	}

	//create a new group inside current one and step into it:
	void HDF5Writer::push(const std::string& name){
		std::vector<std::string> dirlist=splitPathname(name);
    
		//if this is an absolute path, go to top dir first:
		if(name.find_first_of("/")==0){
			cd("/");
		}

		//create groups iteratively:
		hid_t last_group;
		for(unsigned int i=0; i<static_cast<unsigned int>(dirlist.size()); i++){
			last_group=current_group;

			//check if group exists:
			htri_t ex=H5Lexists(last_group,dirlist[i].c_str(),H5P_DEFAULT);
			if(ex==0){
				//does not exists, create a link:
				current_group=H5Gcreate(last_group,dirlist[i].c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
			}
			else{
				//the link exists, check if the group exists as well:
				ex=H5Oexists_by_name(last_group,dirlist[i].c_str(),H5P_DEFAULT);
				if(ex==0){
					QDPIO::cout << "HDF5Writer::push: path traversal error! You seem to have dangling links, I will try to clean up" << std::endl;
					//delete old link:
					H5Ldelete(last_group,dirlist[i].c_str(),H5P_DEFAULT);
					//create new group:
					current_group=H5Gcreate(last_group,dirlist[i].c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
				}
				else{
					//check if object is a group: if yes, step into it, if no, cast an error:
					H5O_info_t objinfo;
					herr_t errhandle=H5Oget_info_by_name(last_group,dirlist[i].c_str(),&objinfo,H5P_DEFAULT);
					if(objinfo.type!=H5O_TYPE_GROUP){
						QDPIO::cout << "HDF5Writer::push: error, the object " << dirlist[i] << " in " << getNameById(last_group) << " already exists and is not a group!"  << std::endl;
						current_group=last_group;
						break;
					}
					else{
						current_group=H5Gopen(last_group,dirlist[i].c_str(),H5P_DEFAULT);
					}
				}
			}
      
			if(current_group<0){
				HDF5_error_exit("HDF5Writer::push: something went wrong, aborting!");
			}
		}
	}
  
	//attribute handling:
	herr_t HDF5Writer::rmAtt(hid_t location_id, const char *attr_name, const H5A_info_t* attrinfo, void* opdata){
		return H5Adelete(location_id,attr_name);
	}

	void HDF5Writer::deleteAllAttributes(const std::string& obj_name){
		H5Aiterate_by_name(current_group,obj_name.c_str(),H5_INDEX_NAME,H5_ITER_NATIVE,NULL,rmAtt,NULL,H5P_DEFAULT);
	}

	void HDF5Writer::deleteAttribute(const std::string& obj_name, const std::string& attr_name){
		std::string oname(obj_name), aname(attr_name);
		htri_t exists=H5Aexists_by_name(current_group,oname.c_str(),aname.c_str(),H5P_DEFAULT);
		if(exists!=1){
			QDPIO::cout << "HDF5Writer::deleteAttribute: error, attribute does not exists!" << std::endl; \
				return;
		} 
		herr_t errhandle=H5Adelete_by_name(current_group,oname.c_str(),aname.c_str(),H5P_DEFAULT);
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//WRITING ATTRIBUTES
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************

	template<typename ctype>
	bool get_global(ctype& global, const ctype& local)
	{
		global = local;
		QDPInternal::broadcast(global);

		int chcksum = global != local ? 1 : 0;
		QDPInternal::globalSum(chcksum);

		return chcksum==0;
	}

	template<typename ctype>
	bool get_global(multi1d<ctype>& global, const multi1d<ctype>& local)
	{
		int s0;
		int chcksum = not get_global(s0, local.size());

		if (Layout::nodeNumber()==0)
			global = local;
		else
			global.resize(s0);

		QDPInternal::broadcast(&(global[0]), sizeof(ctype)*s0);

		if (chcksum!=0) // this is and MUST be global, otherwise a possible deadlock
			return false;
           
		for (int j=0; j<s0; ++j)
			chcksum += global[j] != local[j] ? 1 : 0;

		QDPInternal::globalSum(chcksum);

		return chcksum==0;
	}

	template bool get_global<short>(                multi1d<short>& global,              const multi1d<short>& local);
	template bool get_global<int>(                  multi1d<int>& global,                const multi1d<int>& local);
	template bool get_global<unsigned short>(       multi1d<unsigned short>& global,     const multi1d<unsigned short>& local);
	template bool get_global<unsigned int>(         multi1d<unsigned int>& global,       const multi1d<unsigned int>& local);
	template bool get_global<unsigned long>(        multi1d<unsigned long>& global,      const multi1d<unsigned long>& local);
	template bool get_global<unsigned long long>(   multi1d<unsigned long long>& global, const multi1d<unsigned long long>& local);
	template bool get_global<float>(                multi1d<float>& global,              const multi1d<float>& local);
	template bool get_global<double>(               multi1d<double>& global,             const multi1d<double>& local);

	//single datum
	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const short& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_SHORT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const unsigned short& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_USHORT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const int& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_INT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const unsigned int& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_UINT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const unsigned long long& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_ULLONG,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const float& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_FLOAT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const double& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_DOUBLE,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const std::string& datum, const HDF5Base::writemode& mode){
		std::string oname(obj_name), aname(attr_name);

                std::string datum_0;
                if (not get_global(datum_0,datum)) {
                    QDPIO::cerr << "HDF5Writer::writeAttribute() warning: " << obj_name
                        << ".attrib(" << attr_name << ") was NOT global. Using node=0 value now." << std::endl; 
                }

		bool exists=objectExists(current_group,oname);
		if(!exists){
			HDF5_error_exit("HDF5Writer::writeAttribute: error, object "+oname+" you try to write attribute to does not exists!");
		}

		hid_t ex=H5Aexists_by_name(current_group,oname.c_str(),aname.c_str(),H5P_DEFAULT);
		if(ex==1){
			if(!(mode&HDF5Base::trunc)){
				HDF5_error_exit("HDF5Writer::writeAttribute: error, attribute "+aname+" already exists!");
			}
			herr_t errhandle=H5Adelete_by_name(current_group,oname.c_str(),aname.c_str(),H5P_DEFAULT);
		}

                if (datum_0.length()+1>64*1024) {
                    QDPIO::cerr << "HDF5Writer::writeAttribute() error: " << obj_name
                        << ".attrib(" << attr_name
                        << ") exceeds the maximum hdf5 attrib size (64kB)." << std::endl;

                    HDF5_error_exit("bad string attrib write");
                }

		//create string datatytpe and set encoding to UTF-8:
		hid_t typid=H5Tcreate(H5T_STRING,datum_0.length()+1);
		H5Tset_cset(typid,H5T_CSET_UTF8);
		//create space:
		hid_t attr_space_id=H5Screate(H5S_SCALAR);
		hid_t attr_id=H5Acreate_by_name(current_group,oname.c_str(),aname.c_str(),typid,attr_space_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Sclose(attr_space_id);

		//write:
		H5Awrite(attr_id,typid,reinterpret_cast<void*>(const_cast<char*>(datum_0.c_str())));
		H5Aclose(attr_id);
		H5Tclose(typid);
	}
  
	//array.
	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<short>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_SHORT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<unsigned short>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_USHORT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<int>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_INT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<unsigned int>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_UINT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<unsigned long long>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_ULLONG,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<float>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_FLOAT,mode);
	}

	void HDF5Writer::writeAttribute(const std::string& obj_name, const std::string& attr_name, const multi1d<double>& datum, const HDF5Base::writemode& mode){
		wtAtt(obj_name,attr_name,datum,H5T_NATIVE_DOUBLE,mode);
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//WRITING DATASETS
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//single datum
	void HDF5Writer::write(const std::string& obj_name, const short& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_SHORT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const unsigned short& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_USHORT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const int& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_INT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const unsigned int& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_UINT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const unsigned long long& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_ULLONG,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const float& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_FLOAT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const double& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_DOUBLE,mode);
	}

	void HDF5Writer::write(const std::string& dataname, const std::string& datum, const HDF5Base::writemode& mode){
		std::string dname(dataname);
                std::string datum_0;
                if (not get_global(datum_0,datum)) {
                    QDPIO::cerr << "HDF5Writer::write() warning: " << dataname
                        << " was NOT global. Using node=0 value now." << std::endl; 
                }

		bool exists=objectExists(current_group,dname);
		if(exists){
			if(!(mode&HDF5Base::trunc)){
				QDPIO::cout << "HDF5Writer::write: error, object named " << dname << " already exists!" << std::endl;
				return;
			}
			H5O_info_t objinfo;
			hid_t errhandle=H5Oget_info_by_name(current_group,dname.c_str(),&objinfo,H5P_DEFAULT);
			if(objinfo.type!=H5O_TYPE_DATASET){
				QDPIO::cout << "HDF5Writer::write: error, object you try to write does already exist and is of different type!" << std::endl;
				return;
			}
			errhandle=H5Ldelete(current_group,dname.c_str(),H5P_DEFAULT);
		}

		//create string datatytpe and set encoding to UTF-8:                                                                                                                                                  
		hid_t dataid, spaceid, typid=H5Tcreate(H5T_STRING,datum_0.length()+1);
		H5Tset_cset(typid,H5T_CSET_UTF8);
		//create space:                                                                                                                                                                                       
		spaceid=H5Screate(H5S_SCALAR);
		dataid=H5Dcreate(current_group,dname.c_str(),typid,spaceid,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Sclose(spaceid);

		hid_t plist_id = H5Pcreate (H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		H5Dwrite(dataid,typid,H5S_ALL,H5S_ALL,plist_id,reinterpret_cast<void*>(const_cast<char*>(datum_0.c_str())));
		H5Pclose(plist_id);
		H5Tclose(typid);
		H5Dclose(dataid);
	}

	//array:
	//1D
	void HDF5Writer::write(const std::string& obj_name, const multi1d<short>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_SHORT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi1d<unsigned short>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_USHORT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi1d<int>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_INT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi1d<unsigned int>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_UINT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi1d<unsigned long long>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_ULLONG,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi1d<float>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_FLOAT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi1d<double>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_DOUBLE,mode);
	}
	//2D
	void HDF5Writer::write(const std::string& obj_name, const multi2d<short>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_SHORT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi2d<unsigned short>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_USHORT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi2d<int>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_INT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi2d<unsigned int>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_UINT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi2d<unsigned long long>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_ULLONG,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi2d<float>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_FLOAT,mode);
	}

	void HDF5Writer::write(const std::string& obj_name, const multi2d<double>& datum, const HDF5Base::writemode& mode){
		wt(obj_name,datum,H5T_NATIVE_DOUBLE,mode);
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//WRITING Compound types:                                                                                                                   
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//complex types:
	//Single element:
	template<>
	void HDF5Writer::write< PScalar< PScalar< RComplex<float> > > >(const std::string& dataname, const ComplexF& datum, const HDF5Base::writemode& mode){
		hid_t type_id;
		bool exists=objectExists(file_id,".ComplexFloat");
		if(!exists){
			type_id=createComplexType(sizeof(REAL32));
			commitType(".ComplexFloat",type_id);
		}
		else{
			type_id=H5Topen(file_id,".ComplexFloat",H5P_DEFAULT);
			if(type_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		
		//perform the actual write:
		wt(dataname,datum,type_id,mode);
		H5Tclose(type_id);
	}

	template<>
	void HDF5Writer::write< PScalar< PScalar< RComplex<double> > > >(const std::string& dataname, const ComplexD& datum, const HDF5Base::writemode& mode){
		hid_t type_id;
		bool exists=objectExists(file_id,".ComplexDouble");
		if(!exists){
			type_id=createComplexType(sizeof(REAL64));
			commitType(".ComplexDouble",type_id);
		}
		else{
			type_id=H5Topen(file_id,".ComplexDouble",H5P_DEFAULT);
			if(type_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		//perform the actual write:
		wt(dataname,datum,type_id,mode);
		H5Tclose(type_id);
	}

	//array types
	template<>
	void HDF5Writer::write< PScalar< PScalar< RComplex<float> > > >(const std::string& dataname, const multi1d<ComplexF>& datum, const HDF5Base::writemode& mode){
		hid_t type_id;
		bool exists=objectExists(file_id,".ComplexFloat");
		if(!exists){
			type_id=createComplexType(sizeof(REAL32));
			commitType(".ComplexFloat",type_id);
		}
		else{
			type_id=H5Topen(file_id,".ComplexFloat",H5P_DEFAULT);
			if(type_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		//perform the actual write:
		wt(dataname,datum,type_id,mode);
		H5Tclose(type_id);
	}

	template<>
	void HDF5Writer::write< PScalar< PScalar< RComplex<double> > > >(const std::string& dataname, const multi1d<ComplexD>& datum, const HDF5Base::writemode& mode){
		hid_t type_id;
		bool exists=objectExists(file_id,".ComplexDouble");
		if(!exists){
			type_id=createComplexType(sizeof(REAL64));
			commitType(".ComplexDouble",type_id);
		}
		else{
			type_id=H5Topen(file_id,".ComplexDouble",H5P_DEFAULT);
			if(type_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		//perform the actual write:                                                                                                                                                                           
		wt(dataname,datum,type_id,mode);
		H5Tclose(type_id);
	}

	//ColorMatrix:
	template<>
	void HDF5Writer::write< PScalar< PColorMatrix< RComplex<REAL32>, 3> > >(const std::string& dataname, const ColorMatrixF3& datum, const HDF5Base::writemode& mode){
		//first get complex type:
		hid_t complex_id, colmat_id;
		bool exists=objectExists(file_id,".ComplexFloat");
		if(!exists){
			complex_id=createComplexType(sizeof(REAL32));
			commitType(".ComplexFloat",complex_id);
		}
		else{
			complex_id=H5Topen(file_id,".ComplexFloat",H5P_DEFAULT);
			if(complex_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		//get color matrix type:
		exists=objectExists(file_id,".ColorMatrixFloat3");
		if(!exists){
			colmat_id=createColorMatrixType(complex_id,Nc);
			commitType(".ColorMatrixFloat3",colmat_id);
		}
		else{
			colmat_id=H5Topen(file_id,".ColorMatrixFloat3",H5P_DEFAULT);
			if(colmat_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		wt(dataname,datum,colmat_id,mode);
		H5Tclose(complex_id);
		H5Tclose(colmat_id);
	}

	template<>
	void HDF5Writer::write< PScalar< PColorMatrix< RComplex<REAL64>, 3> > >(const std::string& dataname, const ColorMatrixD3& datum, const HDF5Base::writemode& mode){
		//first get complex type:
		hid_t complex_id, colmat_id;
		bool exists=objectExists(file_id,".ComplexDouble");
		if(!exists){
			complex_id=createComplexType(sizeof(REAL64));
			commitType(".ComplexDouble",complex_id);
		}
		else{
			complex_id=H5Topen(file_id,".ComplexDouble",H5P_DEFAULT);
			if(complex_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		//get color matrix type:
		exists=objectExists(file_id,".ColorMatrixDouble3");
		if(!exists){
			colmat_id=createColorMatrixType(complex_id,Nc);
			commitType(".ColorMatrixDouble3",colmat_id);
		}
		else{
			colmat_id=H5Topen(file_id,".ColorMatrixDouble3",H5P_DEFAULT);
			if(colmat_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}

		wt(dataname,datum,colmat_id,mode);
		H5Tclose(complex_id);
		H5Tclose(colmat_id);
	}

	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//Lattice IO
	//***********************************************************************************************************************************
	//***********************************************************************************************************************************
	//helper routines for Lattice field I/O:
	void HDF5Writer::writePrepare(const std::string& name, const HDF5Base::writemode& mode){
		//before writing is performed, check if dataset exists:
		herr_t errhandle;

		bool exists=objectExists(current_group,name);
		if(exists){
			if(!(mode&HDF5Base::trunc)){
				HDF5_error_exit("HDF5Writer::write: error, dataset already exists and you specified not to overwrite!");
			}
			H5O_info_t objinfo;
			errhandle=H5Oget_info_by_name(current_group,name.c_str(),&objinfo,H5P_DEFAULT);
			if(objinfo.type!=H5O_TYPE_DATASET){
				HDF5_error_exit("HDF5Writer::write: error, the object you try to write does already exist and is of different type!");
			}
			//delete attributes (if present) and unlink storage:
			deleteAllAttributes(name);
			errhandle=H5Ldelete(current_group,name.c_str(),H5P_DEFAULT);
		}
		
		//prefetch for faster I/O:
		if(!isprefetched) prefetchLatticeCoordinates();
	}

	void HDF5Writer::writeLattice(const std::string& name, const hid_t& datatype, const ullong& obj_size, char* buf){
		//writing out reordered data array:
		//determine the dimension of the array: this is useful since not for all classes a write routine will be implemented. In that case,
		//it falls back to a floating point array. In other cases, more sophisticated datatypes will be written and stored
		unsigned int dimension;
		if(obj_size>1) dimension=Nd+1;
		else dimension=Nd;

		//create dataspace and dataset: 
		hsize_t rank = static_cast<hsize_t>(dimension);
		hsize_t* spacesize = new hsize_t[dimension];
		//reorder such that x is fastest: 
		for(unsigned int i = 0; i < Nd; i ++) {
			spacesize[i] = Layout::lattSize()[(Nd - 1) - i];
		} // for                                                                                                                                                                                                                                                                                                                                                              
		if(obj_size>1) spacesize[Nd] = obj_size;
		hid_t filespace = H5Screate_simple(rank, const_cast<const hsize_t*>(spacesize), NULL);
		hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);

		//hyperslab creation and selection:
		hsize_t* node_offset = new hsize_t[dimension];
		hsize_t* offset = new hsize_t[dimension];
		hsize_t* total_count = new hsize_t[dimension];
		hsize_t* dim_size = new hsize_t[dimension];
		//reorder such that x is fastest:                                                                                                                                                                                                                                                                                                                                     
		for(unsigned int i = 0; i < Nd; ++ i){
			dim_size[i] = total_count[i] = Layout::subgridLattSize()[(Nd - 1) - i];
			offset[i] = node_offset[i] = Layout::nodeCoord()[(Nd - 1) - i] * total_count[i];
		} // for                                                                                                                                                                                                                                                                                                                                                              
		if(obj_size>1){
			dim_size[Nd] = total_count[Nd] = obj_size;
			offset[Nd] = node_offset[Nd] = 0;
		}

		//LUSTRE optimization:
		if(stripesize > 0) H5Pset_chunk(dcpl_id,dimension, total_count);

		//create dataset:
		hid_t dset_id = H5Dcreate(current_group, name.c_str(), datatype, filespace,
		H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
		H5Sclose(filespace);
		H5Pclose(dcpl_id);
		delete [] spacesize;

		//create property list for writeout:
		hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

		size_t typesize=H5Tget_size(datatype);
		size_t total_size = typesize * obj_size * Layout::sitesOnNode();
		size_t two_gb = (size_t) 2 * 1024 * 1024 * 1024;

		int power = 0;
		while(total_size * sizeof(REAL) > two_gb) {
			dim_size[0] = dim_size[0] >> 1;         // assuming there is enough room in this dim                                                                                                                                                                                                                                                                                
			total_size = total_size >> 1;
			++ power;
		} // while                                                                                                                                                                                                                                                                                                                                                            
		int blocks = 1 << power;  // number of blocks to write

		// create memspace
		hid_t memspace = H5Screate_simple(rank, dim_size, NULL);
		filespace = H5Dget_space(dset_id);

		for(int i = 0; i < blocks; ++ i) {
			offset[0] = node_offset[0] + i * dim_size[0];
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, const_cast<const hsize_t*>(offset), NULL,
			const_cast<const hsize_t*>(dim_size), NULL);
			// write
			H5Dwrite(dset_id, datatype, memspace, filespace, plist_id, buf + i * total_size);
		} // for

		// cleaning up
		delete [] dim_size;
		delete [] total_count;
		delete [] offset;
		delete [] node_offset;

		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Sclose(memspace);
		H5Pclose(plist_id);
	} // writeLattice() 

	//float lattice color matrix:
	template<>
	void HDF5Writer::write(const std::string& name, const LatticeColorMatrixF3& field, const HDF5Base::writemode& mode){
		StopWatch swatch_prepare, swatch_reorder, swatch_write, swatch_datatypes;
	  
		//before writing is performed, check if dataset exists:
		if(profile) swatch_prepare.start();
		writePrepare(name,mode);
		if(profile) swatch_prepare.stop();

		//color matrix datatype:
		if(profile) swatch_datatypes.start();
		hid_t complex_id, colmat_id;
		bool exists=objectExists(file_id,".ComplexFloat");
		if(!exists){
			complex_id=createComplexType(sizeof(REAL32));
			commitType(".ComplexFloat",complex_id);
		}
		else{
			complex_id=H5Topen(file_id,".ComplexFloat",H5P_DEFAULT);
			if(complex_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		exists=objectExists(file_id,".ColorMatrixFloat3");
		if(!exists){
			colmat_id=createColorMatrixType(complex_id,Nc);
			commitType(".ColorMatrixFloat3",colmat_id);
		}
		else{
			colmat_id=H5Topen(file_id,".ColorMatrixFloat3",H5P_DEFAULT);
			if(colmat_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		if(profile) swatch_datatypes.stop();

		//get node information:
		if(profile) swatch_reorder.start();
		const int mynode=Layout::nodeNumber();
		const int nodeSites = Layout::sitesOnNode();

		//copy buffer into data
		size_t float_size=sizeof(REAL32);
		size_t obj_size=sizeof(ColorMatrixF3)/float_size;
		REAL32* buf=new REAL32[nodeSites*obj_size];
		/*#pragma omp parallel for firstprivate(nodeSites,obj_size,float_size) shared(buf,field)
		for(unsigned int run=0; run<nodeSites; run++){
			memcpy(reinterpret_cast<char*>(buf+run*obj_size),&(field.elem(reordermap[run])),float_size*obj_size);
		}*/
		CvtToHost(reinterpret_cast<void*>(buf),field,nodeSites,float_size*obj_size);
		if(profile) swatch_reorder.stop();

		//write out the stuff:
		if(profile) swatch_write.start();
		writeLattice(name,colmat_id,1,reinterpret_cast<char*>(buf));

		//clean up
		H5Tclose(colmat_id);
		H5Tclose(complex_id);
		delete [] buf;
		if(profile) swatch_write.stop();
	  
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics. Write:" << std::endl;
			QDPIO::cout << "\t preparing: " << swatch_prepare.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t datatype-handling: " << swatch_datatypes.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t reordering: " << swatch_reorder.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t write: " << swatch_write.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB written: " << Layout::vol()*sizeof(ColorMatrixF3)/1024/1024 << std::endl;
		}
	}

	//double lattice color matrix:
	template<>
	void HDF5Writer::write< PScalar< PColorMatrix< RComplex<REAL64>, 3> > >(const std::string& name, const LatticeColorMatrixD3& field, const HDF5Base::writemode& mode)
	{
		StopWatch swatch_prepare, swatch_reorder, swatch_write, swatch_datatypes;
	  
		//before writing is performed, check if dataset exists:
		if(profile) swatch_prepare.start();
		writePrepare(name,mode);
		if(profile) swatch_prepare.stop();

		//color matrix datatype:
		if(profile) swatch_datatypes.start();
		hid_t complex_id, colmat_id;
		bool exists=objectExists(file_id,".ComplexDouble");
		if(!exists){
			complex_id=createComplexType(sizeof(REAL64));
			commitType(".ComplexDouble",complex_id);
		}
		else{
			complex_id=H5Topen(file_id,".ComplexDouble",H5P_DEFAULT);
			if(complex_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		exists=objectExists(file_id,".ColorMatrixDouble3");
		if(!exists){
			colmat_id=createColorMatrixType(complex_id,Nc);
			commitType(".ColorMatrixDouble3",colmat_id);
		}
		else{
			colmat_id=H5Topen(file_id,".ColorMatrixDouble3",H5P_DEFAULT);
			if(colmat_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		if(profile) swatch_datatypes.stop();

		//get node information:
		if(profile) swatch_reorder.start();
		const int mynode=Layout::nodeNumber();
		const int nodeSites = Layout::sitesOnNode();

		//copy buffer into data
		size_t float_size=sizeof(REAL64);
		size_t obj_size=sizeof(ColorMatrixD3)/float_size;
		REAL64* buf=new REAL64[nodeSites*obj_size];
		/*#pragma omp parallel for firstprivate(nodeSites,obj_size,float_size) shared(buf,field)
		for(unsigned int run=0; run<nodeSites; run++){
			memcpy(reinterpret_cast<char*>(buf+run*obj_size),&(field.elem(reordermap[run])),float_size*obj_size);
		}*/
		CvtToHost(reinterpret_cast<void*>(buf),field,nodeSites,float_size*obj_size);
		if(profile) swatch_reorder.stop();
	
		//write out the stuff:
		if(profile) swatch_write.start();
		writeLattice(name,colmat_id,1,reinterpret_cast<char*>(buf));

		//clean up
		H5Tclose(colmat_id);
		H5Tclose(complex_id);
		delete [] buf;
		if(profile) swatch_write.stop();
	  
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics. Write:" << std::endl;
			QDPIO::cout << "\t preparing: " << swatch_prepare.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t datatype-handling: " << swatch_datatypes.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t reordering: " << swatch_reorder.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t write: " << swatch_write.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB written: " << Layout::vol()*sizeof(ColorMatrixD3)/1024/1024 << std::endl;
		}
	}

	//write chroma configuration:
	template<>
	void HDF5Writer::write< PScalar< PColorMatrix< RComplex<REAL64>, 3> > >(const std::string& name, const multi1d<LatticeColorMatrixD3>& field, const HDF5Base::writemode& mode)
	{
		StopWatch swatch_prepare, swatch_reorder, swatch_write, swatch_datatypes;
	  
		//before writing is performed, check if dataset exists:
		if(profile) swatch_prepare.start();
		writePrepare(name,mode);
		if(profile) swatch_prepare.stop();

		//color matrix datatype:
		if(profile) swatch_datatypes.start();
		hid_t complex_id, colmat_id;
		bool exists=objectExists(file_id,".ComplexDouble");
		if(!exists){
			complex_id=createComplexType(sizeof(REAL64));
			commitType(".ComplexDouble",complex_id);
		}
		else{
			complex_id=H5Topen(file_id,".ComplexDouble",H5P_DEFAULT);
			if(complex_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		exists=objectExists(file_id,".ColorMatrixDouble3");
		if(!exists){
			colmat_id=createColorMatrixType(complex_id,Nc);
			commitType(".ColorMatrixDouble3",colmat_id);
		}
		else{
			colmat_id=H5Topen(file_id,".ColorMatrixDouble3",H5P_DEFAULT);
			if(colmat_id<0){
				HDF5_error_exit("HDF5Writer::write: error, cannot open committed Datatype!");
			}
		}
		if(profile) swatch_datatypes.stop();

		//get node information:
		if(profile) swatch_reorder.start();
		const int mynode=Layout::nodeNumber();
		const int nodeSites = Layout::sitesOnNode();

		//copy buffer into data
		size_t float_size=sizeof(REAL64);
		size_t obj_size=sizeof(ColorMatrixD3)/float_size;
		size_t tot_size = nodeSites*field.size()*obj_size;
		REAL64* buf=new REAL64[tot_size];
		unsigned int fsize=field.size();
		/*#pragma omp parallel for firstprivate(nodeSites,fsize,obj_size,float_size) shared(buf,field)
		for(unsigned int run=0; run<nodeSites; run++){
			for(unsigned int dd=0; dd<fsize; dd++){
				memcpy(reinterpret_cast<char*>(buf+(dd+fsize*run)*obj_size),&(field[dd].elem(reordermap[run])),float_size*obj_size);
			}
		}*/
		CvtToHost(reinterpret_cast<void*>(buf),field,nodeSites,fsize,float_size*obj_size);
		if(profile) swatch_reorder.stop();
    
		//write out the stuff:
		if(profile) swatch_write.start();
		writeLattice(name,colmat_id,field.size(),reinterpret_cast<char*>(buf));

		//clean up
		H5Tclose(colmat_id);
		H5Tclose(complex_id);
		delete [] buf;
		if(profile) swatch_write.stop();
	  
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics. Write:" << std::endl;
			QDPIO::cout << "\t preparing: " << swatch_prepare.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t datatype-handling: " << swatch_datatypes.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t reordering: " << swatch_reorder.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t write: " << swatch_write.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB written: " << Layout::vol()*field.size()*sizeof(ColorMatrixD3)/1024/1024 << std::endl;
		}
	}

	//write FUEL configuration
	void HDF5Writer::writeFUEL(const std::string& name, const multi1d<LatticeColorMatrixD3>& field, const HDF5Base::writemode& mode){
		StopWatch swatch_complete;
	  
		QDPIO::cout<< "Writing FUEL config..." << std::flush;
		if(profile) swatch_complete.start();
		if(field.size()!=Nd) HDF5_error_exit("HDF5Writer::writeFUEL: passed vector is not a gauge field!");
		if(objectExists(current_group,name)){
			if(!(mode&HDF5Base::trunc)){
				HDF5_error_exit("HDF5Writer::writeFUEL: error, object "+name+" does already exist!");
			}
			else{
				//check if the object is a group and if it contains the datasets 0 to Nd-1:
				H5O_info_t objinfo;
				herr_t errhandle=H5Oget_info_by_name(current_group,name.c_str(),&objinfo,H5P_DEFAULT);
				if(objinfo.type!=H5O_TYPE_GROUP){
					HDF5_error_exit("HDF5Writer::writeFUEL: error, "+name+" exists but it is not a FUEL config!");
				}
				//check if datasets from 0 to Nd-1 exist:
				for(unsigned int i=0; i<Nd; i++){
					std::stringstream stream;
					stream << name << "/" << i;
					std::string dname=stream.str();
					if(!objectExists(current_group,dname)){
						HDF5_error_exit("HDF5Writer::writeFUEL: error, "+name+" exists but it is not a FUEL config! Dataset "+dname+" was not found!");
					}
				}
				deleteAllAttributes(name);
				H5Ldelete(current_group,name.c_str(),H5P_DEFAULT);
			}
		}

		push(name);
		for(unsigned int i=0; i<Nd; i++){
			std::stringstream stream;
			stream << i;
			std::string dname(stream.str());
			write(dname,field[i],HDF5Base::ate);
			writeAttribute(dname,".kind","LatticeColorMatrix");
		}
		pop();
		if(profile) swatch_complete.stop();
		QDPIO::cout<< "done!" << std::endl;
	  
		if(profile){
			QDPIO::cout << "HDF5-I/O statistics for FUEL-write: " << std::endl;
			QDPIO::cout << "\t total: " << swatch_complete.getTimeInSeconds() << " s." << std::endl;
			QDPIO::cout << "\t MB written: " << Layout::vol()*field.size()*sizeof(ColorMatrixD3)/1024/1024 << std::endl;
		}
	}
}
#endif

