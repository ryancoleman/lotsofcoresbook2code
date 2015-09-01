#ifndef QIOXML_USQCD_GAUGE_H
#define QIOXML_USQCD_GAUGE_H

#ifdef __cplusplus
extern "C"
{
#endif

/*********************************************************************/
/* Top level wrapper for USQCD lattice record XML 

   tag           member           description          
   ------------------------------------------------------------
   
*/

typedef struct {
  QIO_TagCharValue usqcdlatticeinfo_tags;
} QIO_USQCDLatticeInfoWrapper;

#define QIO_USQCD_LATTICE_INFO_WRAPPER {\
  {"usqcdInfo", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD lattice record XML

   tag           member           description          
   ------------------------------------------------------------
   version       version        lattice record version number   1.0
   plaq          plaq           Re Tr U_P/3    average plaquette 
   linktr        linktr         Re Tr U_mu/3   average trace of all gauge links
   info          info           XML string     collaboration option
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagCharValue plaq;
  QIO_TagCharValue linktr;
  QIO_TagCharValue info;
} QIO_USQCDLatticeInfo;


#define QIO_USQCDLATTICEFORMATVERSION "1.0"

#define QIO_USQCD_LATTICE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"plaq"   , "", "", 0}, \
   {"linktr" , "", "", 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/

void QIO_encode_usqcd_lattice_info(QIO_String *record_string, 
				     QIO_USQCDLatticeInfo *record_info);
int QIO_decode_usqcd_lattice_info(QIO_USQCDLatticeInfo *record_info,
				    QIO_String *record_string);

int QIO_insert_usqcdlattice_version(QIO_USQCDLatticeInfo *record_info, char *version);
int QIO_insert_usqcdlattice_plaq( QIO_USQCDLatticeInfo *record_info, char *plaq);
int QIO_insert_usqcdlattice_linktr( QIO_USQCDLatticeInfo *record_info, char *linktr);
int QIO_insert_usqcdlattice_info( QIO_USQCDLatticeInfo *record_info, char *info);
int QIO_insert_usqcdlattice_tag_string(QIO_USQCDLatticeInfoWrapper *wrapper,
					 char *recordinfo_tags);


char *QIO_get_usqcd_lattice_info_tag_string(QIO_USQCDLatticeInfoWrapper *wrapper);
char *QIO_get_plaq(QIO_USQCDLatticeInfo *record_info);
char *QIO_get_linktr(QIO_USQCDLatticeInfo *record_info);
char *QIO_get_info(QIO_USQCDLatticeInfo *record_info);
int QIO_defined_plaq(QIO_USQCDLatticeInfo *record_info);
int QIO_defined_linktr(QIO_USQCDLatticeInfo *record_info);
int QIO_defined_info(QIO_USQCDLatticeInfo *record_info);

QIO_USQCDLatticeInfo *QIO_create_usqcd_lattice_info(char *plaq, char *linktr, char *info);
void QIO_destroy_usqcd_lattice_info(QIO_USQCDLatticeInfo *record_info);

#ifdef __cplusplus
}
#endif

#endif /* QIOXML_USQCD_GAUGE_H */
