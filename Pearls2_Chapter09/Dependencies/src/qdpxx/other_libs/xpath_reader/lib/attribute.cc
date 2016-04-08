#include <xml_attribute.h>
#include <string>
#include <list>
#include <sstream>

using namespace std;
using namespace  XMLWriterAPI;

Attribute::Attribute(const string& _name, const string& _value)
{
  initString(_name, _value);
}

Attribute::Attribute(const string& _name, const int& _value)
{
  init<int>(_name, _value);
}  

Attribute::Attribute(const string& _name, const unsigned int& _value)
{
  init<unsigned int>(_name, _value);
}

Attribute::Attribute(const string& _name, const long int& _value)
{
  init<long int>(_name, _value);
}  

Attribute::Attribute(const string& _name, const unsigned long int& _value)
{
  init<unsigned long int>(_name, _value);
}

Attribute::Attribute(const string& _name, const short int& _value)
{
  init<short int>(_name, _value);
}  

Attribute::Attribute(const string& _name, const unsigned short int& _value)
{
  init<unsigned short int>(_name, _value);
}

Attribute::Attribute(const string& _name, const float& _value)
{
  init<float>(_name, _value);
}  

Attribute::Attribute(const string& _name, const double& _value)
{
  init<double>(_name, _value);
}

Attribute::Attribute(const string& _name, const bool& _value)
{
  init<bool>(_name, _value);
}

Attribute::Attribute(const string& nsprefix, 
		     const string& _name, const string& _value)
{
  initString(nsprefix, _name, _value);
}

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const int& _value)
{
  init<int>(nsprefix, _name, _value);
}  

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const unsigned int& _value)
{
  init<unsigned int>(nsprefix, _name, _value);
}

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const long int& _value)
{
  init<long int>(nsprefix, _name, _value);
}  

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const unsigned long int& _value)
{
  init<unsigned long int>(nsprefix, _name, _value);
}

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const short int& _value)
{
  init<short int>(nsprefix, _name, _value);
}  

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const unsigned short int& _value)
{
  init<unsigned short int>(nsprefix, _name, _value);
}

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const float& _value)
{
  init<float>(nsprefix, _name, _value);
}  

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const double& _value)
{
  init<double>(nsprefix, _name, _value);
}

Attribute::Attribute(const string& nsprefix,
		     const string& _name, const bool& _value)
{
  init<bool>(nsprefix, _name, _value);
}

Attribute::Attribute() {};
Attribute::~Attribute() {};
Attribute::Attribute(const Attribute& a) 
{
  name = a.name;
  value = a.value;
}

const Attribute& 
Attribute::operator=(const Attribute& a) 
{
  this->name = a.name;
  this->value = a.value;
  return *this;
}

string& 
Attribute::getName(void) 
{
  return name;
} 

string& 
Attribute::getValue(void) {
  return value;
}

// Check that the attribute has both a nonzero
// name and a nonzero value.
bool 
Attribute::isEmpty(void) { 
  bool ret_val;
  // Attribute is empty if either its name or value are zero length
  if ( name.size() == 0 || value.size() == 0 ) {
    ret_val = true;
  }
  else { 
    ret_val = false;
  }
  
  return ret_val;
}

template < class T > 
void 
Attribute::init(const string& _name, const T& _value)
{
  name = _name;
  ostringstream printval;
  printval << boolalpha << _value;
  value = printval.str();
}

template < class T > 
void 
Attribute::init(const string& nsprefix, const string& _name, const T& _value)
{
  ostringstream printval;
  
  if( nsprefix.size() != 0 ) { 
    name = nsprefix + ":" + _name;
  } 
  else {
    name = _name;
  }
  
  printval << boolalpha <<  _value;
  value = printval.str();
}

void Attribute::initString(const string& _name, const string& _value)
{
  name = _name;
  value = _value;
}

void 
Attribute::initString(const string& nsprefix, const string& _name, 
		  const string& _value)
{
  
  if( nsprefix.size() != 0 ) { 
    name = nsprefix + ":" + _name;
  } 
  else {
    name = _name;
  }
  
  value = _value;
}
