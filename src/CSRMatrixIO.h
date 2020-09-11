/*
 * CSRMatrixIO.h
 *
 *  Created on: 9 nov. 2015
 *      Author: guignont
 */

#ifndef SRC_MCGS_MATRIXVECTOR_IO_CSRMATRIXIO_H_
#define SRC_MCGS_MATRIXVECTOR_IO_CSRMATRIXIO_H_

namespace HartsSolver
{

template <typename ProfileType>
void exportProfile(Exporter &exporter,Exporter::FileNode base_node,const ProfileType &profile)
{
  typedef Exporter::FileNode FileNode;

  FileNode profile_node = exporter.createFileNode(base_node,"profile");
  FileNode profile_type_node = exporter.createFileNode(profile_node,"profile-type");
  exporter.dump(profile_type_node,std::string("CSR"));
  exporter.closeFileNode(profile_type_node);

  FileNode profile_nlin_node = exporter.createFileNode(profile_node,"nlin");
  int nlin = profile.getNRows();
  exporter.dump(profile_nlin_node, nlin);
  exporter.closeFileNode(profile_nlin_node);

  FileNode profile_ncol_node = exporter.createFileNode(profile_node,"ncol");
  exporter.dump(profile_ncol_node,profile.getNCols());
  exporter.closeFileNode(profile_ncol_node);

  FileNode profile_nelem_node = exporter.createFileNode(profile_node,"nelem");
  int nelem = profile.getNElems();
  exporter.dump(profile_nelem_node,nelem);
  exporter.closeFileNode(profile_nelem_node);

  FileNode profile_line_list_node = exporter.createFileNode(profile_node,"line-list");
  exporter.dump(profile_line_list_node,profile.getKCol(),nlin+1);
  exporter.closeFileNode(profile_line_list_node);

  FileNode profile_col_idx_node = exporter.createFileNode(profile_node,"col-idx");
  exporter.dump(profile_col_idx_node,profile.getCols(),nelem);
  exporter.closeFileNode(profile_col_idx_node);

  exporter.closeFileNode(profile_node);
}

template <typename MatrixType>
void exportMatrix(Exporter &exporter,Exporter::FileNode base_node,const MatrixType &matrix,const std::string &node_name)
{
  typedef Exporter::FileNode FileNode;

  FileNode matrix_node = exporter.createFileNode(base_node,node_name);
  FileNode element_blocksize_node = exporter.createFileNode(matrix_node,"element-blocksize");
  int element_blocksize[2];
  element_blocksize[0] = matrix.getBlockSize();
  element_blocksize[1] = matrix.getBlockSize2();
  exporter.dump(element_blocksize_node,element_blocksize,2);
  exporter.closeFileNode(element_blocksize_node);

  FileNode values_node = exporter.createFileNode(matrix_node,"values");
  exporter.dump(values_node,matrix.getAddressScalarData(),
      matrix.getProfile().getNElems()*matrix.getBlockSize()*matrix.getBlockSize2());
  FileNode values_ordering_node = exporter.createFileNode(values_node,"values-ordering");
  exporter.dump(values_ordering_node,std::string("line-order"));
  exporter.closeFileNode(values_ordering_node);
  exporter.closeFileNode(values_node);

  exportProfile(exporter,matrix_node,matrix.getProfile());

  exporter.closeFileNode(matrix_node);
}

template <typename ProfileType>
void importProfile(Importer &importer,Importer::FileNode base_node,ProfileType &profile)
{
  typedef Importer::FileNode FileNode;

  FileNode profile_node = importer.openFileNode(base_node,"profile");
  FileNode profile_type_node = importer.openFileNode(profile_node,"profile-type");
  std::string profile_type;
  importer.read(profile_type_node,profile_type);
  if(profile_type != "CSR"){
    throw BaseException::RunTimeError("CSR profile expected, found :"+profile_type);
  }
  importer.closeFileNode(profile_type_node);

  FileNode profile_nlin_node = importer.openFileNode(profile_node,"nlin");
  int nlin = 0;
  importer.read(profile_nlin_node,nlin);
  importer.closeFileNode(profile_nlin_node);

  FileNode profile_ncol_node = importer.openFileNode(profile_node,"ncol");
  int ncol = 0;
  importer.read(profile_ncol_node,ncol);
  importer.closeFileNode(profile_ncol_node);

  FileNode profile_nelem_node = importer.openFileNode(profile_node,"nelem");
  int nelem = 0;
  importer.read(profile_nelem_node,nelem);
  importer.closeFileNode(profile_nelem_node);

  profile.init(nlin,ncol,nelem);

  FileNode profile_line_list_node = importer.openFileNode(profile_node,"line-list");
  importer.read(profile_line_list_node,profile.getKCol(),nlin+1);
  importer.closeFileNode(profile_line_list_node);

  FileNode profile_col_idx_node = importer.openFileNode(profile_node,"col-idx");
  importer.read(profile_col_idx_node,profile.getCols(),nelem);
  importer.closeFileNode(profile_col_idx_node);

  importer.closeFileNode(profile_node);

  profile.prepare();
}

template <typename MatrixType>
void importMatrix(Importer &importer,Importer::FileNode base_node,MatrixType &matrix,const std::string &node_name=std::string("matrix"))
{
  typedef typename MatrixType::ProfileType ProfileType;
  typedef Importer::FileNode FileNode;

  ProfileType *profile = new ProfileType();

  FileNode matrix_node = importer.openFileNode(base_node,node_name);

  importProfile(importer,matrix_node,*profile);
  matrix.setProfile(profile,true); // get profile ownership
  matrix.allocate();

  FileNode element_blocksize_node = importer.openFileNode(matrix_node,"element-blocksize");
  int element_blocksize[2] = {0,0};
  importer.read(element_blocksize_node,element_blocksize,2);
  importer.closeFileNode(element_blocksize_node);

  if(element_blocksize[0] != MatrixType::m_block_size || element_blocksize[1] != MatrixType::m_block_size2){
    throw BaseException::RunTimeError("Matrix element block size does not correspond to file element block size");
  }

  FileNode values_node = importer.openFileNode(matrix_node,"values");
  // TODO check for line/column data ordering
  importer.read(values_node,matrix.getAddressScalarData(),
      profile->getNElems()*MatrixType::m_block_size*MatrixType::m_block_size2);
  importer.closeFileNode(values_node);
}

}
#endif /* SRC_MCGS_MATRIXVECTOR_IO_CSRMATRIXIO_H_ */
