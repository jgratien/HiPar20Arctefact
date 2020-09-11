/*
 * HDF5Tools.h
 *
 *  Created on: Oct 13, 2015
 *      Author: gratienj
 */

#ifndef SRC_MCGS_COMMON_UTILS_HDF5TOOLS_H_
#define SRC_MCGS_COMMON_UTILS_HDF5TOOLS_H_

#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef HTS_USE_HDF5
#include "hdf5.h"
#endif

#ifdef HTS_USE_XML2
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#endif

namespace HartsSolver
{
  struct Real3
  {
    Real3(double _x=0.,double _y=0., double _z=0.)
    :x(_x),y(_y),z(_z){}
    double x,y,z;
    bool operator()(const Real3& P1,const Real3& P2) {
      if(P1.x==P2.x)
      {
        if(P1.y==P2.y)
        {
          return P1.z<P2.z ;
        }
        else
          return P1.y<P2.y ;
      }
      else
        return P1.x<P2.x ;
    }
  };

  struct HDF5Base
  {
    typedef enum {
      ASCII,
      HDF5
    } eFormatType  ;

    struct FileNode
    {
      FileNode(int _level=0)
      : name("")
      , path_name("")
      , level(_level)
#ifdef HTS_USE_HDF5
      , h_id(-1)
#endif
#ifdef HTS_USE_XML2
      , x_id(NULL)
#endif
      {}

      std::string name ;
      std::string path_name;
      int         level ;


#ifdef HTS_USE_HDF5
      hid_t       h_id ;
#endif
#ifdef HTS_USE_XML2
      xmlNodePtr  x_id ;
#endif
    };

    class StandardTypes
    {
     public:
      StandardTypes()
     {
#ifdef HTS_USE_HDF5
        hid_t type_id = H5Tcreate(H5T_COMPOUND,sizeof(Real3));
        H5Tinsert(type_id,"X",HOFFSET(Real3,x),H5T_NATIVE_DOUBLE);
        H5Tinsert(type_id,"Y",HOFFSET(Real3,y),H5T_NATIVE_DOUBLE);
        H5Tinsert(type_id,"Z",HOFFSET(Real3,z),H5T_NATIVE_DOUBLE);
        m_real3_id = type_id ;
#endif
      }
      ~StandardTypes(){}
     public:
#ifdef HTS_USE_HDF5
      hid_t nativeType(float) const { return H5T_NATIVE_FLOAT; }
      hid_t nativeType(double) const { return H5T_NATIVE_DOUBLE; }
      hid_t nativeType(Real3) const { return m_real3_id ; }
      hid_t nativeType(long double) const { return H5T_NATIVE_LDOUBLE; }
      hid_t nativeType(unsigned int) const { return H5T_NATIVE_UINT; }
      hid_t nativeType(unsigned long) const { return H5T_NATIVE_ULONG; }
      hid_t nativeType(unsigned long long) const { return H5T_NATIVE_ULLONG; }
      hid_t nativeType(int) const { return H5T_NATIVE_INT; }
      hid_t nativeType(long long) const { return H5T_NATIVE_LLONG; }
      hid_t nativeType(long) const { return H5T_NATIVE_LONG; }
      hid_t nativeType(char) const { return H5T_NATIVE_CHAR; }
      hid_t nativeType(unsigned char) const { return H5T_NATIVE_UCHAR; }
      hid_t nativeType(signed char) const { return H5T_NATIVE_SCHAR; }
      hid_t nativeType(unsigned short) const { return H5T_NATIVE_USHORT; }
      hid_t nativeType(short) const { return H5T_NATIVE_SHORT; }
      hid_t nativeType(std::string) const { return H5T_C_S1; }
      //hid_t nativeType(eDataType sd) const;
#endif
      std::string type(double) const { return "real" ; }
      std::string type(Real3) const { return "real3" ; }
      std::string type(int) const { return "int32" ; }
      //std::string type(Int64) const { return "int64" ; }
      std::string type(std::string) const { return "string" ; }
      std::string type(char) const { return "char" ; }


     public:
#ifdef HTS_USE_HDF5
      hid_t m_char_id; //!< Identifiant HDF des entiers sign�s
      hid_t m_uchar_id; //!< Identifiant HDF des caract�res non-sign�s
      hid_t m_int_id; //!< Identifiant HDF des entiers sign�s
      hid_t m_long_id; //!< Identifiant HDF des entiers long sign�s
      hid_t m_uint_id; //!< Identifiant HDF des entiers non sign�s
      hid_t m_ulong_id; //!< Identifiant HDF des entiers long non sign�s
      hid_t m_real_id; //!< Identifiant HDF des r�els
      hid_t m_real3_id; //!< Identifiant HDF pour les Real3
#endif
    };

    bool is_ascii_truncated(double) const {return true ; }
    bool is_ascii_truncated(Real3) const {return true ; }
    bool is_ascii_truncated(int) const {return false ; }
    //bool is_ascii_truncated(Int64) const {return false ; }
    bool is_ascii_truncated(std::string) const {return false ; }

    HDF5Base(std::string const& name)
    : name(name)
    , xfile_name(name+".xml")
    , hfile_name(name+".h5")
    , format("ascii")
    , type(HDF5)
    {}

    std::string name ;
    std::string xfile_name ;
    std::string hfile_name ;
    std::string format ;
    eFormatType type ;
    StandardTypes m_types ;
#ifdef HTS_USE_XML2
    xmlDocPtr doc ;
#endif
#ifdef HTS_USE_HDF5
    hid_t         hfile ;
#endif
    std::vector<double> rbuffer ;
    std::vector<Real3>  r3buffer ;
    //std::vector<Int64>  i64buffer ;
    std::vector<int>    i32buffer ;

  };

  struct Exporter : public HDF5Base
  {

    Exporter(std::string const& name, std::string const& out_format,int prec)
    : HDF5Base(name)
    , fout(xfile_name.c_str())
    {
      fout<<std::fixed<<std::setprecision(prec) ;
      if(out_format.compare("ascii")==0)
      {
        format = "xml" ;
        type = ASCII ;
      }
      if(out_format.compare("hdf5")==0)
      {
        format = "hdf" ;
        type = HDF5 ;
#ifdef HTS_USE_HDF5
        hfile = H5Fcreate(hfile_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
      }
    }

    FileNode createFileRootNode(int level=0) {
      FileNode node(level) ;
#ifdef HTS_USE_HDF5
      node.h_id = hfile ;
#endif
      return node ;
    }

    FileNode createFileNode(FileNode const& parent, std::string name)
    {
      FileNode node ;
      node.name = name ;
      node.path_name = parent.path_name+"/"+name ;
      node.level = parent.level + 1 ;
      for(int i=0;i<parent.level;++i)
        fout<<"\t";
      fout<<"<"<<name<<">"<<std::endl;
      if(type==HDF5)
      {
#ifdef HTS_USE_HDF5
          node.h_id = H5Gcreate2(parent.h_id,node.name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      }
      return node ;
    }

    FileNode createFileNode(FileNode const& parent, std::string name, std::string group_name)
    {
      FileNode node ;
      node.name = name ;
      node.path_name = parent.path_name+"/"+name ;
      node.level = parent.level + 1 ;
      for(int i=0;i<parent.level;++i)
        fout<<"\t";
      fout<<"<"<<name<<" group-name=\""<<group_name<<"\">"<<std::endl;
      if(type==HDF5)
      {
#ifdef HTS_USE_HDF5
          node.h_id = H5Gcreate2(parent.h_id,node.name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      }
      return node ;
    }

    FileNode createFileNode(FileNode const& parent, std::string name, std::string att_name, std::string att_kind,bool is_group=false)
    {
      FileNode node ;
      node.name = name ;
      if(is_group)
        node.path_name = parent.path_name+"/"+name+"_"+att_kind+"_"+att_name ;
      else
        node.path_name = parent.path_name+"/"+att_name ;
      node.level = parent.level + 1 ;
      for(int i=0;i<parent.level;++i)
        fout<<"\t";
      fout<<"<"<<name<<" name=\""<<att_name<<"\" kind=\""<<att_kind<<"\">"<<std::endl;
      if(type==HDF5)
      {
#ifdef HTS_USE_HDF5
          node.h_id = H5Gcreate2(parent.h_id,att_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      }
      return node ;
    }

    void closeFileNode(FileNode const& group)
    {
      for(int i=0;i<group.level-1;++i)
        fout<<"\t";
      fout<<"</"<<group.name<<">"<<std::endl;
      if(type==HDF5)
      {
#ifdef HTS_USE_HDF5
        H5Gclose(group.h_id);
#endif
      }
    }

    void dump(FileNode const& parent,std::vector<Real3>& buffer, int nb_elems_per_line=1)
    {
      switch(type)
      {
        case ASCII :
        {
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"<data format=\"xml\"  type=\"real3\">"<<std::endl;
          int icount = 0 ;
          for (unsigned n=0;n<buffer.size()/nb_elems_per_line;++n)
          {
            for(int k=0;k<nb_elems_per_line;++k)
              fout<<buffer[icount+k].x<<" "<<buffer[icount+k].y<<" "<<buffer[icount+k].z<<" ";
            icount += nb_elems_per_line ;
            fout<<std::endl ;
          }
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        case HDF5 :
        {
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"<data format=\"hdf\"  type=\"real3\">"<<std::endl;
          for(int i=0;i<parent.level+1;++i)
            fout<<"\t";
          fout<<hfile_name<<":"<<parent.path_name<<std::endl ;
#ifdef HTS_USE_HDF5
          hsize_t dim = buffer.size() ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);

          hid_t dataset_id = H5Dcreate2(parent.h_id, "data", m_types.nativeType(Real3()), dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          herr_t status = H5Dwrite(dataset_id,m_types.nativeType(Real3()), H5S_ALL, H5S_ALL, H5P_DEFAULT,buffer.data());
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;
          H5Dclose(dataset_id);
          H5Sclose(dataspace_id);
#endif
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;

        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
      buffer.clear() ;
    }

    template<typename ValueT>
    void dump(FileNode const& parent,std::vector<ValueT>& buffer, int nb_elems_per_line=1)
    {
      switch(type)
      {
        case ASCII :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"xml\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          int icount = 0 ;
          for (size_t n=0;n<buffer.size()/nb_elems_per_line;++n)
          {
            for(int k=0;k<nb_elems_per_line;++k)
              fout<<buffer[icount+k]<<" ";
            icount += nb_elems_per_line ;
            fout<<std::endl ;
          }
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        case HDF5 :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"hdf\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          for(int i=0;i<parent.level+1;++i)
            fout<<"\t";
          fout<<hfile_name<<":"<<parent.path_name<<std::endl ;
#ifdef HTS_USE_HDF5
          hsize_t dim = buffer.size() ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
          hid_t dataset_id = H5Dcreate2(parent.h_id, "data", m_types.nativeType(ValueT()), dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          herr_t status = H5Dwrite(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,buffer.data());
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;

          H5Dclose(dataset_id);
          H5Sclose(dataspace_id);
#endif
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
      buffer.clear() ;
    }

    template<typename ValueT>
    void dump(FileNode const& parent,const std::vector<ValueT>& buffer, int nb_elems_per_line=1)
    {
      switch(type)
      {
        case ASCII :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"xml\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          int icount = 0 ;
          for (int n=0;n<buffer.size()/nb_elems_per_line;++n)
          {
            for(int k=0;k<nb_elems_per_line;++k)
              fout<<buffer[icount+k]<<" ";
            icount += nb_elems_per_line ;
            fout<<std::endl ;
          }
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        case HDF5 :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"hdf\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          for(int i=0;i<parent.level+1;++i)
            fout<<"\t";
          fout<<hfile_name<<":"<<parent.path_name<<std::endl ;
#ifdef HTS_USE_HDF5
          hsize_t dim = buffer.size() ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
          hid_t dataset_id = H5Dcreate2(parent.h_id, "data", m_types.nativeType(ValueT()), dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          herr_t status = H5Dwrite(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,buffer.data());
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;

          H5Dclose(dataset_id);
          H5Sclose(dataspace_id);
#endif
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }

    template<typename ValueT>
    void dump(FileNode const& parent,const ValueT *buffer,const int nb_elem, int nb_elems_per_line=1)
    {
      switch(type)
      {
        case ASCII :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"xml\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          int icount = 0 ;
          for (int n=0;n<nb_elem/nb_elems_per_line;++n)
          {
            for(int k=0;k<nb_elems_per_line;++k)
              fout<<buffer[icount+k]<<" ";
            icount += nb_elems_per_line ;
            fout<<std::endl ;
          }
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        case HDF5 :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"hdf\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          for(int i=0;i<parent.level+1;++i)
            fout<<"\t";
          fout<<hfile_name<<":"<<parent.path_name<<std::endl ;
#ifdef HTS_USE_HDF5
          hsize_t dim = nb_elem;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
          hid_t dataset_id = H5Dcreate2(parent.h_id, "data", m_types.nativeType(ValueT()), dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          herr_t status = H5Dwrite(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,buffer);
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;

          H5Dclose(dataset_id);
          H5Sclose(dataspace_id);
#endif
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }

    template<typename ValueT>
    void dump(FileNode const& parent,const ValueT &val)
    {
      switch(type)
      {
        case ASCII :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"xml\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          fout<<val<<std::endl;

          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        case HDF5 :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"hdf\"  type=\""<<m_types.type(ValueT())<<"\">"<<std::endl;

          for(int i=0;i<parent.level+1;++i)
            fout<<"\t";
          fout<<hfile_name<<":"<<parent.path_name<<std::endl ;
#ifdef HTS_USE_HDF5
          hsize_t dim = 1 ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
          hid_t dataset_id = H5Dcreate2(parent.h_id, "data", m_types.nativeType(ValueT()), dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          herr_t status = H5Dwrite(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,&val);
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;
          H5Dclose(dataset_id);
          H5Sclose(dataspace_id);
#endif
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }


    void dump(FileNode const& parent,const std::string &s)
    {
      switch(type)
      {
        case ASCII :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"xml\"  type=\"string\">"<<std::endl;
          fout<<s<<std::endl;
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        case HDF5 :
        {
          for(int i=0;i<parent.level;++i)
          fout<<"\t";
          fout<<"<data format=\"hdf\"  type=\"string\">"<<std::endl;

          for(int i=0;i<parent.level+1;++i)
            fout<<"\t";
          fout<<hfile_name<<":"<<parent.path_name<<std::endl ;
#ifdef HTS_USE_HDF5
          hsize_t dim = s.size() ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
          hid_t dataset_id = H5Dcreate2(parent.h_id, "data", m_types.nativeType(std::string()),
              dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          herr_t status = H5Dwrite(dataset_id, m_types.nativeType(std::string()), H5S_ALL, H5S_ALL, H5P_DEFAULT,s.c_str());
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;
          H5Dclose(dataset_id);
          H5Sclose(dataspace_id);
#endif
          for(int i=0;i<parent.level;++i)
            fout<<"\t";
          fout<<"</data>"<<std::endl;
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }

    std::ofstream fout ;
  };

  struct Importer : public HDF5Base
  {

    Importer(std::string const& name, std::string const& in_format,int prec)
    : HDF5Base(name)
    {
      if(in_format.compare("ascii")==0)
      {
        format = "xml" ;
        type = ASCII ;
#ifdef HTS_USE_XML2
        doc = xmlParseFile(xfile_name.c_str());
        if(doc == NULL)
        {
          std::cerr<<"Error while parsing XML file : "<<xfile_name<<std::endl ;
        }
#endif
      }
      if(in_format.compare("hdf5")==0)
      {
        format = "hdf" ;
        type = HDF5 ;
#ifdef HTS_USE_HDF5
        hfile = H5Fopen(hfile_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
#endif
      }
    }

    FileNode openFileRootNode(int level=0) {
      FileNode node(level) ;
      switch(type)
      {
        case ASCII:
#ifdef HTS_USE_XML2
          node.x_id =  xmlDocGetRootElement(doc) ;
          //if (xmlStrcmp(node.x_id->name, (const xmlChar *) "form-function"))
          //  std::cerr<<"Error bad root name"<<std::endl;
#endif
          break ;
        case HDF5 :
#ifdef HTS_USE_HDF5
          node.h_id = hfile ;
#endif
          break ;
      }
      return node ;
    }

    void closeFileRootNode(FileNode const& node) {
      switch(type)
      {
        case ASCII:
#ifdef HTS_USE_XML2
          xmlFreeDoc(doc);
#endif
          break ;
        case HDF5:
#ifdef HTS_USE_HDF5
          herr_t status = H5Fclose(node.h_id) ;
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;
#endif
          break ;
      }
    }


    FileNode openFileNode(FileNode const& parent, std::string name)
    {
      FileNode node ;
      node.name = name ;
      node.path_name = parent.path_name+"/"+name ;
      node.level = parent.level + 1 ;
      switch(type)
      {
        case ASCII:
        {
#ifdef HTS_USE_XML2
          xmlNodePtr cur = parent.x_id->xmlChildrenNode;
          while (cur != NULL)
          {
            if ((!xmlStrcmp(cur->name, (const xmlChar *)name.c_str())))
            {
              node.x_id = cur ;
              break ;
            }
            cur = cur->next;
          }
#endif
        }
        break ;
        case HDF5:
        {
#ifdef HTS_USE_HDF5
          /* Open an existing dataset. */
          node.h_id = H5Gopen2(parent.h_id,node.name.c_str(), H5P_DEFAULT);
#endif
        }
        break ;
      }
      return node ;
    }

    void closeFileNode(FileNode const& group)
    {
      switch(type)
      {
        case ASCII :
        {
          //xmlFree(group.x_id) ;
        }
        break ;
        case HDF5 :
        {
#ifdef HTS_USE_HDF5
          H5Gclose(group.h_id);
#endif
        }
        break ;
      }
    }

    template<typename ValueT>
    void read(FileNode const& parent,std::vector<ValueT>& buffer)
    {
      switch(type)
      {
        case ASCII:
          {
#ifdef HTS_USE_XML2
            xmlChar* contenu = xmlNodeGetContent(parent.x_id);
            //printf("%s\n",contenu);
            std::stringstream flux ;
            flux << contenu ;
            for(std::size_t i=0;i<buffer.size();++i)
              flux>>buffer[i];
            xmlFree(contenu);
#endif
          }
          break ;
        case HDF5 :
        {
#ifdef HTS_USE_HDF5
          //TODO: remove
          hsize_t dim = buffer.size() ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);

          /* Open an existing dataset. */
          hid_t dataset_id = H5Dopen2(parent.h_id, "data", H5P_DEFAULT);

          herr_t status = H5Dread(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,buffer.data());
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;
          /* Close the dataset. */
          H5Dclose(dataset_id);

          H5Sclose(dataspace_id);
#endif
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }

    template<typename ValueT>
    void read(FileNode const& parent,ValueT& val)
    {
      switch(type)
      {
        case ASCII:
          {
#ifdef HTS_USE_XML2
            xmlChar* contenu = xmlNodeGetContent(parent.x_id);
            //printf("%s\n",contenu);
            std::stringstream flux ;
            flux<< contenu ;
            flux>>val;
            xmlFree(contenu);
#endif
          }
          break ;
        case HDF5 :
        {
#ifdef HTS_USE_HDF5
          // TODO: remove
          hsize_t dim = 1 ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);

          /* Open an existing dataset. */
          hid_t dataset_id = H5Dopen2(parent.h_id, "data", H5P_DEFAULT);

          herr_t status = H5Dread(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,&val);
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;

          /* Close the dataset. */
          H5Dclose(dataset_id);

          H5Sclose(dataspace_id);
#endif
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }

    template<typename ValueT>
    void read(FileNode const& parent,ValueT *buffer,const int n_elem)
    {
      switch(type)
      {
        case ASCII:
          {
#ifdef HTS_USE_XML2
            xmlChar *contenu = xmlNodeGetContent(parent.x_id);
            //printf("%s\n",contenu);
            std::stringstream flux ;
            flux<<contenu ;
            for(int i=0;i<n_elem;++i)
              flux>>buffer[i];
            xmlFree(contenu);
#endif
          }
          break ;
        case HDF5 :
        {
#ifdef HTS_USE_HDF5
          // TODO: remove
          hsize_t dim = n_elem ;
          hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);

          /* Open an existing dataset. */
          hid_t dataset_id = H5Dopen2(parent.h_id, "data", H5P_DEFAULT);

          herr_t status = H5Dread(dataset_id, m_types.nativeType(ValueT()), H5S_ALL, H5S_ALL, H5P_DEFAULT,buffer);
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;

          /* Close the dataset. */
          H5Dclose(dataset_id);

          H5Sclose(dataspace_id);
#endif
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }


    void read(FileNode const& parent,std::string& val)
    {
      switch(type)
      {
        case ASCII:
          {
#ifdef HTS_USE_XML2
            xmlChar *contenu = xmlNodeGetContent(parent.x_id);
            //printf("%s\n",contenu);
            std::stringstream flux ;
            flux<<contenu ;
            flux>>val;
            xmlFree(contenu);
#endif
          }
          break ;
        case HDF5 :
        {
#ifdef HTS_USE_HDF5
          /* Open an existing dataset. */
          hid_t dataset_id = H5Dopen2(parent.h_id, "data", H5P_DEFAULT);
          // get dataset size
          hid_t dataspace_id = H5Dget_space(dataset_id);
          int ndim = H5Sget_simple_extent_ndims(dataspace_id);
          hsize_t dims[ndim];
          H5Sget_simple_extent_dims(dataspace_id,dims,NULL);
          char string_on_file[dims[0]+1];
          string_on_file[dims[0]] = '\0';
          //cout << "string size=" << dims[0] << endl;
          herr_t status = H5Dread(dataset_id, m_types.nativeType(std::string()), H5S_ALL, H5S_ALL, H5P_DEFAULT,string_on_file);
          if(status) std::cout<<"HDF5 ERROR : "<<status<<std::endl ;
          val = std::string(string_on_file);
          /* Close the dataset. and dataspace */
          H5Sclose(dataspace_id);
          H5Dclose(dataset_id);
#endif
        }
        break ;
        default :
          std::cerr<<"Unknown output format "<<std::endl ;
          break ;
      }
    }

  };


} /* namespace MCGSolver */

#endif /* SRC_MCGS_COMMON_UTILS_HDF5TOOLS_H_ */
