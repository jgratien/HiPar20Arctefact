/*
 * MatrixGenerator.h
 *
 *  Created on: Mar 15, 2014
 *      Author: gratienj
 */

#ifndef MATRIXGENERATOR_H_
#define MATRIXGENERATOR_H_

#include <iostream>
#include <cassert>
#include <cmath>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

namespace HartsSolver
{

  class MatrixGenerator
  {
  public:

    class BoundaryCondition
    {
    public :
      typedef enum {
        Dirichlet,
        Neunman,
        NbConditionType
      } eConditionType;

      typedef enum {
        Xmin,
        Xmax,
        Ymin,
        Ymax,
        Zmin,
        Zmax,
        NbBoundaryType
      } eBoundaryType;

      BoundaryCondition()
      : m_type(Neunman)
      , m_value(0.)
      {}
      eConditionType m_type ;
      double m_value ;
    };
    MatrixGenerator(){}
    virtual ~MatrixGenerator(){}

    template<typename MatrixT, typename VectorT>
    void buildLaplacian(MatrixT& matrix,
        int nx, int ny,
        VectorT& rhs,
        BoundaryCondition const* bc,
        bool use_perm_file,
        std::string& perm_file)
    {
      typedef typename MatrixT::InfoVectorType InfoVectorType;
      using namespace MatrixVector ;
      int nrows = nx*ny;
      int nnz = 5*(nx-2)*(ny-2) + (nx+ny-4)*8+4*3 ;
      std::cout<<"NROWS : "<<nrows<<std::endl ;
      std::cout<<"NNZ       : "<<nnz<<std::endl ;
      m_permitivity.resize(nrows) ;
      if(use_perm_file)
      {
          std::ifstream fin(perm_file.c_str());
          if(fin.is_open())
              for(int k=0; k<nrows; ++k)
                  fin>>m_permitivity[k];
          fin.close();
      }
      else
          std::fill(m_permitivity.data(),m_permitivity.data()+nrows,1.) ;

      typename MatrixT::ProfileType& profile = matrix.getProfile() ;
      profile.init(nrows,nnz) ;
      matrix.allocate() ;
      rhs.resize(nrows) ;
      rhs.assign(nrows,0.) ;
      InfoVectorType& cols = profile.getCols() ;
      InfoVectorType& kcol = profile.getKCol() ;
      typedef typename MatrixT::MatrixDataType MatrixDataType ;
      typedef Diag<MatrixDataType> OffDiagBlockPatternType ;
      //typedef Diag<MatrixDataType> BlockPatternType ;
      typedef UTri<MatrixDataType> BlockPatternType ;
      //typedef LTri<MatrixDataType> BlockPatternType ;
      typedef ScalarId<typename MatrixT::VectorDataType> Block1DPatternType ;
      typename MatrixT::MatrixDataType* values = matrix.getAddressData() ;
      int irow =0 ;
      int offset = 0 ;
      m_nx = nx ;
      m_ny = ny ;
      {
        int j=0 ;
        {
          int i=0 ;
          double T_m_i = _trans_m_i(i,j) ;
          double T_p_i = _trans_p_i(i,j) ;
          double T_m_j = _trans_m_j(i,j) ;
          double T_p_j = _trans_p_j(i,j) ;

          int row_size = 3 ;
          kcol[irow] = offset ;
          cols[offset] = irow ;
          cols[offset+1] = irow+1 ;
          cols[offset+2] = irow+nx ;
          values[offset] = BlockPatternType::D(T_p_i + T_p_j) ;
          if(bc[BoundaryCondition::Xmin].m_type == BoundaryCondition::Dirichlet)
          {
            values[offset] += BlockPatternType::D(T_m_i);
            rhs[0] += Block1DPatternType::D(T_m_i*bc[BoundaryCondition::Xmin].m_value) ;
          }
          if(bc[BoundaryCondition::Ymin].m_type == BoundaryCondition::Dirichlet)
          {
            values[offset] += BlockPatternType::D(T_m_j);
            rhs[0] += Block1DPatternType::D(T_m_j*bc[BoundaryCondition::Ymin].m_value) ;
          }
          values[offset+1] = OffDiagBlockPatternType::D(-T_p_i) ;
          values[offset+2] = OffDiagBlockPatternType::D(-T_p_j) ;
          offset += row_size ;
          ++irow ;
        }
        for(int i=1;i<nx-1;++i)
        {
          double T_m_i = _trans_m_i(i,j) ;
          double T_p_i = _trans_p_i(i,j) ;
          double T_m_j = _trans_m_j(i,j) ;
          double T_p_j = _trans_p_j(i,j) ;

          int row_size = 4 ;
          kcol[irow] = offset ;
          cols[offset] = irow-1 ;
          cols[offset+1] = irow ;
          cols[offset+2] = irow+1 ;
          cols[offset+3] = irow+nx ;
          values[offset] = OffDiagBlockPatternType::D(-T_m_i) ;
          values[offset+1] = BlockPatternType::D(T_m_i+T_p_i+T_p_j) ;
          if(bc[BoundaryCondition::Ymin].m_type == BoundaryCondition::Dirichlet)
          {
            values[offset+1] += BlockPatternType::D(T_m_j);
            rhs[i] += Block1DPatternType::D(T_m_j*bc[BoundaryCondition::Ymin].m_value) ;
          }
          values[offset+2] = OffDiagBlockPatternType::D(-T_p_i) ;
          values[offset+3] = OffDiagBlockPatternType::D(-T_p_j) ;
          offset += row_size ;
          ++irow ;
        }
        {
            int i=nx-1 ;
            double T_m_i = _trans_m_i(i,j) ;
            double T_p_i = _trans_p_i(i,j) ;
            double T_m_j = _trans_m_j(i,j) ;
            double T_p_j = _trans_p_j(i,j) ;

            int row_size = 3 ;
            kcol[irow] = offset ;
            cols[offset] = irow-1 ;
            cols[offset+1] = irow ;
            cols[offset+2] = irow+nx ;
            values[offset] = OffDiagBlockPatternType::D(-T_m_i) ;
            values[offset+1] = BlockPatternType::D(T_m_i+T_p_j) ;
            values[offset+2] = OffDiagBlockPatternType::D(-T_p_j) ;
            if(bc[BoundaryCondition::Xmax].m_type == BoundaryCondition::Dirichlet)
            {
              values[offset+1] += BlockPatternType::D(T_p_i);
              rhs[nx-1] += Block1DPatternType::D(T_p_i*bc[BoundaryCondition::Xmax].m_value) ;
            }
            if(bc[BoundaryCondition::Ymin].m_type == BoundaryCondition::Dirichlet)
            {
              values[offset+1] += BlockPatternType::D(T_m_j);
              rhs[nx-1] += Block1DPatternType::D(T_m_j*bc[BoundaryCondition::Ymin].m_value) ;
            }
            offset += row_size ;
            ++irow ;
         }
      }
      for(int j=1;j<ny-1;++j)
      {
        {
            int i=0 ;
            double T_m_i = _trans_m_i(i,j) ;
            double T_p_i = _trans_p_i(i,j) ;
            double T_m_j = _trans_m_j(i,j) ;
            double T_p_j = _trans_p_j(i,j) ;

            int row_size = 4 ;
            kcol[irow] = offset ;
            cols[offset] = irow-nx ;
            cols[offset+1] = irow ;
            cols[offset+2] = irow+1 ;
            cols[offset+3] = irow+nx ;
            values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
            values[offset+1] = BlockPatternType::D(T_m_j+T_p_i+T_p_j) ;
            values[offset+2] = OffDiagBlockPatternType::D(-T_p_i) ;
            values[offset+3] = OffDiagBlockPatternType::D(-T_p_j) ;
            if(bc[BoundaryCondition::Xmin].m_type == BoundaryCondition::Dirichlet)
            {
              values[offset+1] += BlockPatternType::D(T_m_i);
              rhs[j*nx] += Block1DPatternType::D(T_m_i*bc[BoundaryCondition::Xmin].m_value) ;
            }
            offset += row_size ;
            ++irow ;
         }
        for(int i=1;i<nx-1;++i)
        {
          int row_size = 5 ;
          double T_m_i = _trans_m_i(i,j) ;
          double T_p_i = _trans_p_i(i,j) ;
          double T_m_j = _trans_m_j(i,j) ;
          double T_p_j = _trans_p_j(i,j) ;

          kcol[irow] = offset ;
          cols[offset] = irow-nx ;
          cols[offset+1] = irow -1;
          cols[offset+2] = irow ;
          cols[offset+3] = irow+1 ;
          cols[offset+4] = irow+nx ;
          values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
          values[offset+1] = OffDiagBlockPatternType::D(-T_m_i) ;
          values[offset+2] = BlockPatternType::D(T_m_j+T_m_i+T_p_i+T_p_j) ;
          values[offset+3] = OffDiagBlockPatternType::D(-T_p_i) ;
          values[offset+4] = OffDiagBlockPatternType::D(-T_p_j) ;
          offset += row_size ;
          ++irow ;
        }
        {
            int i=nx-1 ;
            double T_m_i = _trans_m_i(i,j) ;
            double T_p_i = _trans_p_i(i,j) ;
            double T_m_j = _trans_m_j(i,j) ;
            double T_p_j = _trans_p_j(i,j) ;

            int row_size = 4 ;
            kcol[irow] = offset ;
            cols[offset] = irow-nx ;
            cols[offset+1] = irow-1 ;
            cols[offset+2] = irow ;
            cols[offset+3] = irow+nx ;
            values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
            values[offset+1] = OffDiagBlockPatternType::D(-T_m_i) ;
            values[offset+2] = BlockPatternType::D(T_m_j+T_m_i+T_p_j) ;
            values[offset+3] = OffDiagBlockPatternType::D(-T_p_j) ;

            if(bc[BoundaryCondition::Xmax].m_type == BoundaryCondition::Dirichlet)
            {
              values[offset+2] += BlockPatternType::D(T_p_i);
              rhs[j*nx+nx-1] += Block1DPatternType::D(T_p_i*bc[BoundaryCondition::Xmax].m_value) ;
            }
            offset += row_size ;
            ++irow ;
         }
      }
      {
        int j=ny-1 ;
        {
          int i=0 ;
          double T_m_i = _trans_m_i(i,j) ;
          double T_p_i = _trans_p_i(i,j) ;
          double T_m_j = _trans_m_j(i,j) ;
          double T_p_j = _trans_p_j(i,j) ;

          int row_size = 3 ;
          kcol[irow] = offset ;
          cols[offset] = irow-nx ;
          cols[offset+1] = irow ;
          cols[offset+2] = irow+1 ;
          values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
          values[offset+1] = BlockPatternType::D(T_m_j+T_p_i) ;
          values[offset+2] = OffDiagBlockPatternType::D(-T_p_i) ;
          if(bc[BoundaryCondition::Xmin].m_type == BoundaryCondition::Dirichlet)
          {
            values[offset+1] += BlockPatternType::D(T_m_i);
            rhs[(ny-1)*nx] += Block1DPatternType::D(T_m_i*bc[BoundaryCondition::Xmin].m_value) ;
          }
          if(bc[BoundaryCondition::Ymax].m_type == BoundaryCondition::Dirichlet)
          {
            values[offset+1] += BlockPatternType::D(T_p_j);
            rhs[(ny-1)*nx] += Block1DPatternType::D(T_p_j*bc[BoundaryCondition::Ymax].m_value) ;
          }
          offset += row_size ;
          ++irow ;
        }
        for(int i=1;i<nx-1;++i)
        {
          double T_m_i = _trans_m_i(i,j) ;
          double T_p_i = _trans_p_i(i,j) ;
          double T_m_j = _trans_m_j(i,j) ;
          double T_p_j = _trans_p_j(i,j) ;

          int row_size = 4 ;
          kcol[irow] = offset ;
          cols[offset] = irow-nx ;
          cols[offset+1] = irow-1 ;
          cols[offset+2] = irow ;
          cols[offset+3] = irow+1 ;
          values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
          values[offset+1] = OffDiagBlockPatternType::D(-T_m_i) ;
          values[offset+2] = BlockPatternType::D(T_m_j+T_m_i+T_p_i) ;
          values[offset+3] = OffDiagBlockPatternType::D(-T_p_i) ;

          if(bc[BoundaryCondition::Ymax].m_type == BoundaryCondition::Dirichlet)
          {
            values[offset+2] += BlockPatternType::D(T_p_j);
            rhs[(ny-1)*nx+i] += Block1DPatternType::D(T_p_j*bc[BoundaryCondition::Ymax].m_value) ;
          }
          offset += row_size ;
          ++irow ;
        }
        {
            int i=nx-1 ;
            double T_m_i = _trans_m_i(i,j) ;
            double T_p_i = _trans_p_i(i,j) ;
            double T_m_j = _trans_m_j(i,j) ;
            double T_p_j = _trans_p_j(i,j) ;

            int row_size = 3 ;
            kcol[irow] = offset ;
            cols[offset] = irow-nx ;
            cols[offset+1] = irow -1 ;
            cols[offset+2] = irow ;
            values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
            values[offset+1] = OffDiagBlockPatternType::D(-T_m_i) ;
            values[offset+2] = BlockPatternType::D(T_m_j+T_m_i) ;
            if(bc[BoundaryCondition::Xmax].m_type == BoundaryCondition::Dirichlet)
            {
              values[offset+2] += BlockPatternType::D(T_p_i);
              rhs[ny*nx-1] += Block1DPatternType::D(T_p_i*bc[BoundaryCondition::Ymax].m_value) ;
            }
            if(bc[BoundaryCondition::Ymax].m_type == BoundaryCondition::Dirichlet)
            {
              values[offset+2] += BlockPatternType::D(T_p_j);
              rhs[ny*nx-1] += Block1DPatternType::D(T_p_j*bc[BoundaryCondition::Ymax].m_value) ;
            }
            offset += row_size ;
            ++irow ;
         }
      }
      kcol[irow] = offset ;
      std::cout<<"NROW : "<<irow<<" NNZ : "<<offset<<std::endl ;
    }

    template<typename MatrixT, typename VectorT>
    void buildLaplacian(MatrixT& matrix,
        int nx, int ny, int nz,
        VectorT& rhs,
        BoundaryCondition const* bc,
        bool use_perm_file,
        std::string& perm_file)
    {
      typedef typename MatrixT::InfoVectorType InfoVectorType;
      using namespace MatrixVector ;
      int nrows = nx*ny*nz;
      int nnz = 7*(nx-2)*(ny-2)*(nz-2) + ((nx-2)*(ny-2) +(nx-2)*(nz-2) +(ny-2)*(nz-2))*12 + ((nx-2)+(ny-2)+(nz-2))*20 +8*4 ;
      std::cout<<"NROWS : "<<nrows<<std::endl ;
      std::cout<<"NNZ   : "<<nnz<<std::endl ;
      m_permitivity.resize(nrows) ;
      if(use_perm_file)
      {
          std::ifstream fin(perm_file.c_str());
          if(fin.is_open())
              for(int k=0; k<nrows; ++k)
                  fin>>m_permitivity[k];
          fin.close();
      }
      else
          std::fill(m_permitivity.data(),m_permitivity.data()+nrows,1.) ;

      typename MatrixT::ProfileType& profile = matrix.getProfile() ;
      profile.init(nrows,nnz) ;
      matrix.allocate() ;
      rhs.resize(nrows) ;
      rhs.assign(nrows,0.) ;
      InfoVectorType& cols = profile.getCols() ;
      InfoVectorType& kcol = profile.getKCol() ;
      typedef typename MatrixT::MatrixDataType MatrixDataType ;
      typedef Diag<MatrixDataType> OffDiagBlockPatternType ;
      //typedef Diag<MatrixDataType> BlockPatternType ;
      typedef UTri<MatrixDataType> BlockPatternType ;
      //typedef LTri<MatrixDataType> BlockPatternType ;
      typedef ScalarId<typename MatrixT::VectorDataType> Block1DPatternType ;
      typename MatrixT::MatrixDataType* values = matrix.getAddressData() ;
      m_nx = nx ;
      m_ny = ny ;
      m_nz = nz ;
      m_nxy = nx*ny ;

      int irow =0 ;
      int offset = 0 ;
      for(int k=0;k<nz;++k)
      {
        for(int j=0;j<ny;++j)
        {
          for(int i=0;i<nx;++i)
          {
            double T_m_i = _trans_m_i(i,j,k) ;
            double T_p_i = _trans_p_i(i,j,k) ;
            double T_m_j = _trans_m_j(i,j,k) ;
            double T_p_j = _trans_p_j(i,j,k) ;
            double T_m_k = _trans_m_k(i,j,k) ;
            double T_p_k = _trans_p_k(i,j,k) ;

            kcol[irow] = offset ;
            int diag_offset = offset + diagOffset(i,j,k) ;
            if(isInternal(k-1,nz))
            {
              cols[offset] = irow - m_nxy ;
              values[offset] = OffDiagBlockPatternType::D(-T_m_k) ;
              ++offset ;
            }
            else
            {
              if(bc[BoundaryCondition::Zmin].m_type == BoundaryCondition::Dirichlet)
              {
                values[diag_offset] += BlockPatternType::D(T_m_k);
                rhs[uid(i,j,k)]     += Block1DPatternType::D(T_m_k*bc[BoundaryCondition::Zmin].m_value) ;
              }
            }
            if(isInternal(j-1,ny))
            {
              cols[offset] = irow - m_nx ;
              values[offset] = OffDiagBlockPatternType::D(-T_m_j) ;
              ++offset ;
            }
            else
            {
              if(bc[BoundaryCondition::Ymin].m_type == BoundaryCondition::Dirichlet)
              {
                values[diag_offset] += BlockPatternType::D(T_m_j);
                rhs[uid(i,j,k)]     += Block1DPatternType::D(T_m_j*bc[BoundaryCondition::Ymin].m_value) ;
              }
            }
            if(isInternal(i-1,nx))
            {
              cols[offset] = irow - 1 ;
              values[offset] = OffDiagBlockPatternType::D(-T_m_i) ;
              ++offset ;
            }
            else
            {
              if(bc[BoundaryCondition::Xmin].m_type == BoundaryCondition::Dirichlet)
              {
                values[diag_offset] += BlockPatternType::D(T_m_i);
                rhs[uid(i,j,k)]     += Block1DPatternType::D(T_m_i*bc[BoundaryCondition::Xmin].m_value) ;
              }
            }
            {
              cols[offset] = irow ;
              values[offset] = BlockPatternType::D(T_m_k+T_m_j+T_m_i+T_p_i+T_p_j+T_m_k) ;
              ++offset;
            }
            if(isInternal(i+1,nx))
            {
              cols[offset] = irow + 1 ;
              values[offset] = OffDiagBlockPatternType::D(-T_p_i) ;
              ++offset ;
            }
            else
            {
              if(bc[BoundaryCondition::Xmax].m_type == BoundaryCondition::Dirichlet)
              {
                values[diag_offset] += BlockPatternType::D(T_p_i);
                rhs[uid(i,j,k)]     += Block1DPatternType::D(T_m_i*bc[BoundaryCondition::Xmax].m_value) ;
              }
            }
            if(isInternal(j+1,ny))
            {
              cols[offset] = irow + m_nx ;
              values[offset] = OffDiagBlockPatternType::D(-T_p_j) ;
              ++offset ;
            }
            else
            {
              if(bc[BoundaryCondition::Ymax].m_type == BoundaryCondition::Dirichlet)
              {
                values[diag_offset] += BlockPatternType::D(T_p_j);
                rhs[uid(i,j,k)]     += Block1DPatternType::D(T_m_j*bc[BoundaryCondition::Ymax].m_value) ;
              }
            }
            if(isInternal(k+1,nz))
            {
              cols[offset] = irow + m_nxy ;
              values[offset] = OffDiagBlockPatternType::D(-T_p_j) ;
              ++offset ;
            }
            else
            {
              if(bc[BoundaryCondition::Zmax].m_type == BoundaryCondition::Dirichlet)
              {
                values[diag_offset] += BlockPatternType::D(T_p_k);
                rhs[uid(i,j,k)]     += Block1DPatternType::D(T_m_k*bc[BoundaryCondition::Zmax].m_value) ;
              }
            }
            ++irow ;
          }
        }
      }
      assert(offset==nnz) ;
      kcol[irow] = offset ;
      std::cout<<"NROW : "<<irow<<" NNZ : "<<offset<<std::endl ;
    }



    template<typename MatrixT>
    void buildDenseMatrix(MatrixT& matrix,  int nrows, bool print)
    {
      typename MatrixT::ProfileType& profile = matrix.getProfile() ;
      profile.init(nrows,nrows*nrows) ;
      matrix.allocate() ;
      int* kcol = HartsSolver::kcol(profile) ;
      int* cols = HartsSolver::cols(profile) ;
      typename MatrixT::ValueType* values = HartsSolver::dataPtr(matrix) ;
      int offset = 0 ;
      for(int i=0;i<nrows;++i)
      {
        kcol[i] = offset ;
        for(int j=0;j<nrows;++j)
        {
          cols[offset+j] = j ;
          values[offset+j] = -1./(1+i+j) ;
        }
        values[offset+i] = 1. ;
        offset += nrows ;
      }
      kcol[nrows] = offset ;
      if(print)
      {
        std::cout<<"GENERATE DENSE MATRIX : "<<std::endl ;
        matrix.print() ;
      }
    }

    template<typename MatrixT>
    int readFromFile(MatrixT& matrix,std::string& filename,int nb_file=1)
    {
      using namespace std ;
      if(nb_file==1)
      {
        std::ifstream fin(filename.c_str());
        int nrows, nnz;
        fin>>nrows>>ws>>nnz;
        typename MatrixT::ProfileType& profile = matrix.getProfile() ;
        profile.init(nrows,nnz) ;
        matrix.allocate() ;

        int* kcol = HartsSolver::kcol(profile) ;
        int* cols = HartsSolver::cols(profile) ;
        typename MatrixT::MatrixDataType* values = HartsSolver::dataPtr(matrix) ;
        for(int i=0;i<nrows+1;++i)
          fin>>kcol[i]>>ws ;
        for(int i=0;i<nnz;++i)
          fin>>cols[i]>>ws;
        for(int i=0;i<nnz;++i)
          fin>>values[i]>>ws ;
      }
      else if(nb_file>1)
      {
        int nrows=0 ;
        int nnz=0;
        std::vector<ifstream*> files(nb_file) ;
        std::vector<int> local_nrows(nb_file) ;
        std::vector<int> local_nnz(nb_file) ;
        for(int i=0;i<nb_file;++i)
        {
          stringstream fname ;
          fname<<filename<<"."<<i<<"P";
          std::cout<<"READING FILE : "<<fname.str()<<std::endl ;
          files[i] = new ifstream(fname.str().c_str());
          ifstream& fin = *files[i] ;
          fin>>local_nrows[i]>>ws>>local_nnz[i];
          nrows += local_nrows[i] ;
          nnz += local_nnz[i] ;
        }
        typename MatrixT::ProfileType& profile = matrix.getProfile() ;
        profile.init(nrows,nnz) ;
        matrix.allocate() ;

        int* kcol = HartsSolver::kcol(profile) ;
        int* cols = HartsSolver::cols(profile) ;
        typename MatrixT::ValueType* values = HartsSolver::dataPtr(matrix) ;
        int row_offset = 0 ;
        int offset = 0 ;
        for(int i=0;i<nb_file;++i)
        {
          ifstream& fin = *files[i] ;
          int lnrows = local_nrows[i] ;
          int lnnz = local_nnz[i] ;
          for(int i=0;i<lnrows+1;++i)
          {
            fin>>kcol[row_offset+i]>>ws ;
            kcol[row_offset+i] += offset ;
          }
          for(int i=0;i<lnnz;++i)
            fin>>cols[offset+i]>>ws;
          for(int i=0;i<lnnz;++i)
            fin>>values[offset+i]>>ws ;
          row_offset += lnrows ;
          offset += lnnz ;
        }
      }
      return 0 ;
    }

    template<typename ValueT,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
    int readFromFileB(CSRMatrix<ValueT,1,AllocatorT> & matrix,std::string& filename, bool sum_first_eq=false)
    {
      using namespace std ;
      try
      {
        typedef CSRMatrix<ValueT,1,AllocatorT> MatrixType ;
        ifstream fin(filename.c_str());
        int nrows, nnz, block1d_size;
        fin>>nrows>>ws>>nnz>>ws>>block1d_size;
        int block2d_size = block1d_size*block1d_size ;
        typename MatrixType::ProfileType& profile = matrix.getProfile() ;
        profile.init(nrows,nnz) ;
        matrix.allocate() ;

        int* kcol = HartsSolver::kcol(profile) ;
        int* cols = HartsSolver::cols(profile) ;
        typename MatrixType::MatrixDataType* values = HartsSolver::dataPtr(matrix) ;
        for(int i=0;i<nrows+1;++i)
          fin>>kcol[i]>>ws ;
        for(int i=0;i<nnz;++i)
          fin>>cols[i]>>ws;
        for(int i=0;i<nnz;++i)
        {
          fin>>values[i]>>ws ;
          if(sum_first_eq)
          {
            ValueT tmp ;
            for(int k2=1;k2<block1d_size;++k2)
              fin>>tmp>>ws ;
            for(int k1=1;k1<block1d_size;++k1)
            {
              fin>>tmp>>ws ;
              values[i] += tmp ;
              for(int k2=1;k2<block1d_size;++k2)
                fin>>tmp>>ws ;
            }
          }
          else
          {
            ValueT tmp ;
            for(int k=1;k<block2d_size;++k)
              fin>>tmp>>ws ;
          }
        }
        profile.computeDCol() ;
      }
      catch(std::exception exc)
      {
        cerr<<"Error while importing matrix : "<<exc.what() ;
        return 1 ;
      }
      return 0 ;
    }

    template<typename ValueT, std::size_t N,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
    int readFromFileB(CSRMatrix<ValueT,N,AllocatorT> & matrix,std::string& filename, bool sum_first_eq=false)
    {
      using namespace std ;
      {
        typedef CSRMatrix<ValueT,N> MatrixType ;
        ifstream fin(filename.c_str());
        std::size_t nrows, nnz, block1d_size;
        fin>>nrows>>ws>>nnz>>ws>>block1d_size;
        if(block1d_size==0)
        {
          std::cerr<<" Null Matrix Data Block Size";
          return 1 ;
        }
        if(block1d_size!=N)
          std::cout<<"Warning Matrix Data Block Size "<<block1d_size<<" not equal to "<<N<<std::endl ;
        std::size_t block2d_size = std::min(block1d_size*block1d_size,N*N) ;
        typename MatrixType::ProfileType& profile = matrix.getProfile() ;
        profile.init(nrows,nnz) ;
        matrix.allocate() ;

        int* kcol = HartsSolver::kcol(profile) ;
        int* cols = HartsSolver::cols(profile) ;
        std::size_t NxN = N*N;
        typename MatrixType::MatrixDataType* block_values = HartsSolver::dataPtr(matrix) ;
        typename MatrixType::ValueType* values = block_values->begin() ;
        for(std::size_t i=0;i<nrows+1;++i)
          fin>>kcol[i]>>ws ;
        for(std::size_t i=0;i<nnz;++i)
          fin>>cols[i]>>ws;
        int off=0 ;
        for(std::size_t i=0;i<nnz;++i)
        {
          for(std::size_t j=0;j<block2d_size;++j)
          fin>>values[off+j]>>ws ;
          if(sum_first_eq)
          {
            for(std::size_t ki=1;ki<block1d_size;++ki)
              for(std::size_t kj=0;kj<block1d_size;++kj)
                values[off+kj] += values[off+ki*block1d_size+kj] ;
          }
          off += NxN ;
        }
        profile.computeDCol() ;
      }
      return 0 ;
    }




    template<typename ValueT,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
     int readFromFileXML(CSRMatrix<ValueT,1,AllocatorT>& matrix,
                        typename CSRMatrix<ValueT,1,AllocatorT>::VectorType& rhs,
                        std::string const& filename,
                        std::string const& format,
                        bool sum_first_eq=false,
                        int prec=6)
     {
      typedef CSRMatrix<ValueT,1,AllocatorT> MatrixType ;

      typename MatrixType::ProfileType& profile = matrix.getProfile() ;

      typedef Exporter::FileNode FileNode ;

      Importer importer(filename,format,prec) ;
      std::vector<double>& rbuffer  = importer.rbuffer ;
      std::vector<int>& i32buffer   = importer.i32buffer ;

      int nrows = 0 ;
      int nnz = 0 ;
      int blk_size = 1 ;

      FileNode root_node = importer.openFileRootNode() ;
      //FileNode system_node = importer.openFileNode(root_node,"system") ;
      {
        FileNode matrix_node = importer.openFileNode(root_node,"matrix") ;
        {
          FileNode info_node = importer.openFileNode(matrix_node,"struct-info") ;
          {
            FileNode nrows_node = importer.openFileNode(info_node,"nrows") ;
            i32buffer.resize(1) ;
            importer.read(nrows_node,i32buffer) ;
            importer.closeFileNode(nrows_node) ;
            nrows = i32buffer[0] ;
          }
          {
            FileNode nnz_node = importer.openFileNode(info_node,"nnz") ;
            i32buffer.resize(1) ;
            importer.read(nnz_node,i32buffer) ;
            importer.closeFileNode(nnz_node) ;
            nnz = i32buffer[0] ;
          }
          profile.init(nrows,nnz) ;
          matrix.allocate() ;
          int* kcol = HartsSolver::kcol(profile) ;
          int* cols = HartsSolver::cols(profile) ;
          {
            FileNode blk_size_node = importer.openFileNode(info_node,"blk-size") ;
            i32buffer.resize(1) ;
            importer.read(blk_size_node,i32buffer) ;
            importer.closeFileNode(blk_size_node) ;
            blk_size = i32buffer[0] ;
          }
          {
            FileNode kcol_node = importer.openFileNode(info_node,"kcol") ;
            i32buffer.resize(nrows+1) ;
            importer.read(kcol_node,i32buffer) ;
            importer.closeFileNode(kcol_node) ;
            for(Integer i=0;i<nrows+1;++i)
            {
              kcol[i] = i32buffer[i];
            }
          }
          {
            FileNode cols_node = importer.openFileNode(info_node,"cols") ;
            i32buffer.resize(nnz) ;
            importer.read(cols_node,i32buffer) ;
            importer.closeFileNode(cols_node) ;
            for(Integer i=0;i<nnz;++i)
            {
              cols[i] = i32buffer[i];
            }
          }
          importer.closeFileNode(info_node) ;
          profile.computeDCol() ;

          FileNode data_node = importer.openFileNode(matrix_node,"data") ;
          {
            typename MatrixType::MatrixDataType* values = HartsSolver::dataPtr(matrix) ;

            FileNode values_node = importer.openFileNode(data_node,"values") ;
            rbuffer.resize(nnz*blk_size*blk_size) ;
            importer.read(values_node,rbuffer) ;
            importer.closeFileNode(values_node) ;
            for(int k=0;k<nnz;++k)
            {
              values[k] = rbuffer[k*blk_size*blk_size] ;
               if(sum_first_eq)
               {
                 for(int i=1;i<blk_size;++i)
                   values[k] += rbuffer[k*blk_size*blk_size + i*blk_size] ;
               }
            }
          }
          importer.closeFileNode(data_node) ;
        }
      importer.closeFileNode(matrix_node) ;
    }


      {
        FileNode vector_node = importer.openFileNode(root_node,"vector") ;
        {
          FileNode info_node = importer.openFileNode(vector_node,"struct-info") ;
          {
            FileNode nrows_node = importer.openFileNode(info_node,"nrows") ;
            i32buffer.resize(1) ;
            importer.read(nrows_node,i32buffer) ;
            importer.closeFileNode(nrows_node) ;
            nrows = i32buffer[0] ;
          }
          {
            FileNode blk_size_node = importer.openFileNode(info_node,"blk-size") ;
            i32buffer.resize(1) ;
            importer.read(blk_size_node,i32buffer) ;
            importer.closeFileNode(blk_size_node) ;
            blk_size = i32buffer[0] ;
          }
          importer.closeFileNode(info_node) ;

          rhs.resize(nrows) ;
          FileNode data_node = importer.openFileNode(vector_node,"data") ;
          {
            FileNode values_node = importer.openFileNode(data_node,"values") ;
            rbuffer.resize(nrows*blk_size) ;
            importer.read(values_node,rbuffer) ;
            importer.closeFileNode(values_node) ;
            for(int i=0;i<nrows;++i)
            {
               rhs[i] = rbuffer[i*blk_size] ;
               if(sum_first_eq)
               {
                 for(int j=1;j<blk_size;++j)
                   rhs[i] += rbuffer[i*blk_size+j] ;
               }
            }
          }
          importer.closeFileNode(data_node) ;
        }
      importer.closeFileNode(vector_node) ;
    }

     //importer.closeFileNode(system_node) ;
     importer.closeFileRootNode(root_node) ;
     return 0 ;
    }

    template<typename ValueT, std::size_t N,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
     int readFromFileXML(CSRMatrix<ValueT,N,AllocatorT>& matrix,
                        typename CSRMatrix<ValueT,N,AllocatorT>::VectorType& rhs,
                        std::string const& filename,
                        std::string const& format,
                        bool sum_first_eq=false,
                        int prec=6)
     {
       typedef CSRMatrix<ValueT,N> MatrixType ;

       typename MatrixType::ProfileType& profile = matrix.getProfile() ;

       typedef Exporter::FileNode FileNode ;

       Importer importer(filename,format,prec) ;
       std::vector<double>& rbuffer  = importer.rbuffer ;
       std::vector<int>& i32buffer   = importer.i32buffer ;

       int nrows = 0 ;
       int nnz = 0 ;
       int blk_size = 1 ;

       FileNode root_node = importer.openFileRootNode() ;
       //FileNode system_node = importer.openFileNode(root_node,"system") ;
       {
         FileNode matrix_node = importer.openFileNode(root_node,"matrix") ;
         {
           FileNode info_node = importer.openFileNode(matrix_node,"struct-info") ;
           {
             FileNode nrows_node = importer.openFileNode(info_node,"nrows") ;
             i32buffer.resize(1) ;
             importer.read(nrows_node,i32buffer) ;
             importer.closeFileNode(nrows_node) ;
             nrows = i32buffer[0] ;
           }
           {
             FileNode nnz_node = importer.openFileNode(info_node,"nnz") ;
             i32buffer.resize(1) ;
             importer.read(nnz_node,i32buffer) ;
             importer.closeFileNode(nnz_node) ;
             nnz = i32buffer[0] ;
           }
           profile.init(nrows,nnz) ;
           matrix.allocate() ;
           int* kcol = HartsSolver::kcol(profile) ;
           int* cols = HartsSolver::cols(profile) ;
           {
             FileNode blk_size_node = importer.openFileNode(info_node,"blk-size") ;
             i32buffer.resize(1) ;
             importer.read(blk_size_node,i32buffer) ;
             importer.closeFileNode(blk_size_node) ;
             blk_size = i32buffer[0] ;
           }
           {
             FileNode kcol_node = importer.openFileNode(info_node,"kcol") ;
             i32buffer.resize(nrows+1) ;
             importer.read(kcol_node,i32buffer) ;
             importer.closeFileNode(kcol_node) ;
             for(Integer i=0;i<nrows+1;++i)
             {
               kcol[i] = i32buffer[i];
             }
           }
           {
             FileNode cols_node = importer.openFileNode(info_node,"cols") ;
             i32buffer.resize(nnz) ;
             importer.read(cols_node,i32buffer) ;
             importer.closeFileNode(cols_node) ;
             for(Integer i=0;i<nnz;++i)
             {
               cols[i] = i32buffer[i];
             }
           }
           importer.closeFileNode(info_node) ;
           profile.computeDCol() ;

           FileNode data_node = importer.openFileNode(matrix_node,"data") ;
           {
             typename MatrixType::MatrixDataType* values = HartsSolver::dataPtr(matrix) ;

             FileNode values_node = importer.openFileNode(data_node,"values") ;
             rbuffer.resize(nnz*blk_size*blk_size) ;
             importer.read(values_node,rbuffer) ;
             importer.closeFileNode(values_node) ;
             for(int k=0;k<nnz;++k)
             {
               for(std::size_t i=0;i<N;++i)
                 for(std::size_t j=0;j<N;++j)
                   values[k][i][j] = rbuffer[k*blk_size*blk_size+i*blk_size+j] ;
             }
           }
           importer.closeFileNode(data_node) ;
         }
       importer.closeFileNode(matrix_node) ;
     }


       {
         FileNode vector_node = importer.openFileNode(root_node,"vector") ;
         {
           FileNode info_node = importer.openFileNode(vector_node,"struct-info") ;
           {
             FileNode nrows_node = importer.openFileNode(info_node,"nrows") ;
             i32buffer.resize(1) ;
             importer.read(nrows_node,i32buffer) ;
             importer.closeFileNode(nrows_node) ;
             nrows = i32buffer[0] ;
           }
           {
             FileNode blk_size_node = importer.openFileNode(info_node,"blk-size") ;
             i32buffer.resize(1) ;
             importer.read(blk_size_node,i32buffer) ;
             importer.closeFileNode(blk_size_node) ;
             blk_size = i32buffer[0] ;
           }
           importer.closeFileNode(info_node) ;

           rhs.resize(nrows) ;
           FileNode data_node = importer.openFileNode(vector_node,"data") ;
           {
             FileNode values_node = importer.openFileNode(data_node,"values") ;
             rbuffer.resize(nrows*blk_size) ;
             importer.read(values_node,rbuffer) ;
             importer.closeFileNode(values_node) ;
             for(int i=0;i<nrows;++i)
               for(std::size_t k=0;k<N;++k)
                rhs[i][k] = rbuffer[i*blk_size+k] ;
           }
           importer.closeFileNode(data_node) ;
         }
       importer.closeFileNode(vector_node) ;
     }

      //importer.closeFileNode(system_node) ;
      importer.closeFileRootNode(root_node) ;

      return 0 ;
     }

    template<typename ValueT,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
    int readFromFileIJV(CSRMatrix<ValueT,1,AllocatorT>& matrix,
                        std::string& filename,
                        bool cnumbering=true,
                        int stride=1)
    {
      std::cout<<"STRIDE = "<<stride<<std::endl ;
      typedef typename CSRMatrix<ValueT,1,AllocatorT>::MatrixEntryType MatrixEntryType ;
      std::ifstream fin(filename.c_str()) ;
      int nrows, nnz ;
      int incr = cnumbering?0:-1 ;
      fin>>nrows>>std::ws>>nnz;
      if(stride!=1)
      {
        nrows /= stride ;
        nnz /= stride*stride ;
      }
      std::vector<MatrixEntryType> entries;
      entries.reserve(nnz) ;
      while(!fin.eof())
      {
        int i,j;
        ValueT value ;
        fin>>i>>std::ws>>j>>std::ws>>value;
        if(stride==1)
          entries.push_back(MatrixEntryType(i+incr,j+incr,value)) ;
        else
        {
          if( ((i+incr)%stride==0) && ((j+incr)%stride==0) )
            entries.push_back(MatrixEntryType((i+incr)/stride,(j+incr)/stride,value)) ;
        }
      }
      std::cout<<"NROWS : "<<nrows<<std::endl ;
      setFromTriplets(matrix,nrows,entries) ;
      return 0 ;
    }


    template<typename ValueT,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
    int readFromFileIJV(CSRMatrix<ValueT,1,AllocatorT>& matrix,
                        std::string& filename,
                        bool cnumbering,
                        int dof_type,
                        std::size_t nb_dofs,
                        std::vector<int> const& gdof_lid,
                        std::vector<int> const& gdof_type)
    {
      std::set<int> dof_types ;
      switch(dof_type)
      {
      case 0:
        dof_types.insert(0) ;
        break ;
      case 1:
        dof_types.insert(1) ;
        break ;
      case 2:
        dof_types.insert(2) ;
        break ;
      case 3:
        dof_types.insert(0) ;
        dof_types.insert(1) ;
        break ;
      case 4:
        dof_types.insert(0) ;
        dof_types.insert(2) ;
        break ;
      case 5:
        dof_types.insert(1) ;
        dof_types.insert(2) ;
        break ;
      }
      typedef typename CSRMatrix<ValueT,1,AllocatorT>::MatrixEntryType MatrixEntryType ;
      std::ifstream fin(filename.c_str()) ;
      int nrows, nnz ;
      int incr = cnumbering?0:-1 ;
      fin>>nrows>>std::ws>>nnz;
      std::vector<MatrixEntryType> entries;
      entries.reserve(nnz) ;
      while(!fin.eof())
      {
        int i,j;
        ValueT value ;
        fin>>i>>std::ws>>j>>std::ws>>value;
        if((dof_types.find(gdof_type[i+incr])!=dof_types.end()) &&
           (dof_types.find(gdof_type[j+incr])!=dof_types.end()))
        {
          assert(gdof_lid[i+incr]!=-1) ;
          assert(gdof_lid[j+incr]!=-1) ;
            entries.push_back(MatrixEntryType(gdof_lid[i+incr],gdof_lid[j+incr],value)) ;
        }
      }
     setFromTriplets(matrix,nb_dofs,entries) ;
     return 0 ;
    }

    template<typename ValueT,typename AllocatorT=HartsSolver::MatrixVector::DefaultAllocator>
    int readFromFileMXM(CSRMatrix<ValueT,1,AllocatorT>& matrix,
                       typename CSRMatrix<ValueT,1,AllocatorT>::VectorType& rhs,
                       std::string& filename)
    {
      using namespace std ;
      typedef CSRMatrix<ValueT,1,AllocatorT> MatrixType ;    FILE *fdes;
      size_t size;
      char *line = nullptr;
      fdes = fopen(filename.c_str(),"r");

      if(getline(&line,&size,fdes)==-1){
        perror(filename.c_str());
        throw BaseException::RunTimeError(__FILE__,__LINE__);
      }
      // skip comments
      while (line[0]=='%'){
        std::free(line);
        line = nullptr;
        if(getline(&line,&size,fdes)==-1){
          perror(filename.c_str());
          throw BaseException::RunTimeError(__FILE__,__LINE__);
        }
      }
      // first non comment line is: n m nnz
      int nrows,m,nnz;
      sscanf(line,"%d %d %d",&nrows,&m,&nnz);

      typename MatrixType::ProfileType& profile = matrix.getProfile() ;
      profile.init(nrows,nnz) ;
      matrix.allocate() ;
      int* kcol = HartsSolver::kcol(profile) ;
      int* col = HartsSolver::cols(profile) ;
      typename MatrixType::MatrixDataType* csr_val = HartsSolver::dataPtr(matrix) ;

      std::vector<int> li(nnz);
      std::vector<int> ci(nnz);
      std::vector<ValueT> val(nnz);

      for(int i=0;i<nnz;++i){
        fscanf(fdes,"%d %d %lg\n",&li[i],&ci[i],&val[i]);
        li[i]--;
        ci[i]--;
      }
      fclose(fdes);

      for(int i=0;i<=nrows;++i){
        kcol[i] = 0;
      }
      for(int i=0;i<nnz;++i){
        kcol[1+li[i]]++;
      }
      for(int i=1;i<=nrows;++i){
        kcol[i]+=kcol[i-1];
      }
      std::vector<int> kcolidx(nrows);
      for(int i=0;i<nrows;++i){
        kcolidx[i] = 0;
      }
      for(int i=0;i<nnz;++i){
        int line = li[i];
        int kci = kcol[line]+kcolidx[line];
        kcolidx[line]++;
        csr_val[kci] = val[i];
        col[kci] = ci[i];
      }
      return 0 ;
    }

private :
    int hat(int i,int n)
    {
      return std::max(0,std::min(i,n-1)) ;
    }
    int uid(int i,int j)
    {
      return hat(j,m_ny)*m_nx+hat(i,m_nx) ;
    }
    int uid(int i,int j,int k)
    {
      return hat(k,m_nz)*m_nxy+hat(j,m_ny)*m_nx+hat(i,m_nx) ;
    }
    int isInternal(int i,int n)
    {
      if((i>=0)&&(i<n))
        return 1 ;
      else
        return 0 ;
    }
    int diagOffset(int i,int j,int k)
    {
      return isInternal(k-1,m_nz)+isInternal(j-1,m_ny)+isInternal(i-1,m_nx) ;
    }
    double  _trans_m_i(int i,int j)
    {
      double p1 = m_permitivity[uid(i-1,j)] ;
      double p2 = m_permitivity[uid(i,j)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_p_i(int i,int j)
    {
      double p1 = m_permitivity[uid(i+1,j)] ;
      double p2 = m_permitivity[uid(i,j)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_m_j(int i,int j)
    {
      double p1 = m_permitivity[uid(i,j-1)] ;
      double p2 = m_permitivity[uid(i,j)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_p_j(int i,int j)
    {
      double p1 = m_permitivity[uid(i,j+1)] ;
      double p2 = m_permitivity[uid(i,j)] ;
      return p1*p2/(p1+p2) ;
    }

    double  _trans_m_i(int i,int j,int k)
    {
      double p1 = m_permitivity[uid(i-1,j,k)] ;
      double p2 = m_permitivity[uid(i,j,k)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_p_i(int i,int j,int k)
    {
      double p1 = m_permitivity[uid(i+1,j,k)] ;
      double p2 = m_permitivity[uid(i,j,k)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_m_j(int i,int j,int k)
    {
      double p1 = m_permitivity[uid(i,j-1,k)] ;
      double p2 = m_permitivity[uid(i,j,k)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_p_j(int i,int j,int k)
    {
      double p1 = m_permitivity[uid(i,j+1,k)] ;
      double p2 = m_permitivity[uid(i,j,k)] ;
      return p1*p2/(p1+p2) ;
    }

    double  _trans_m_k(int i,int j,int k)
    {
      double p1 = m_permitivity[uid(i,j,k-1)] ;
      double p2 = m_permitivity[uid(i,j,k)] ;
      return p1*p2/(p1+p2) ;
    }
    double  _trans_p_k(int i,int j,int k)
    {
      double p1 = m_permitivity[uid(i,j,k+1)] ;
      double p2 = m_permitivity[uid(i,j,k)] ;
      return p1*p2/(p1+p2) ;
    }


    int m_nx = 0 ;
    int m_ny = 0 ;
    int m_nz = 0 ;
    int m_nxy = 0 ;

    std::vector<double> m_permitivity;

  };

} /* namespace HartsSolver */

#endif /* MATRIXGENERATOR_H_ */
