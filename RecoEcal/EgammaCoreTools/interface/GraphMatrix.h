#ifndef RecoEcal_EgammaCoreTools_GraphMatrix_h
#define RecoEcal_EgammaCoreTools_GraphMatrix_h

/**
   \file
   Tools for manipulating ECAL Clusters as graphs
   \author Davide Valsecchi, Badder Marzocchi
   \date 05 October 2020
*/

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip> 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
namespace ublas = boost::numeric::ublas;
typedef size_t size_type; 

template<typename T>
class GraphMatrix; 

template<typename T>
std::ostream& operator<<(ostream& s, const GraphMatrix<T>&m);

template<typename T>
class GraphMatrix 
{ 
  friend std::ostream &operator<< <>(ostream& s, const GraphMatrix&m);
  typedef typename ublas::matrix<T>::const_iterator1 const_iterator1;
  typedef typename ublas::matrix<T>::const_iterator2 const_iterator2;
  typedef typename ublas::matrix<T>::iterator1 iterator1;
  typedef typename ublas::matrix<T>::iterator2 iterator2;
  typedef typename ublas::matrix<T>::const_reverse_iterator1 const_reverse_iterator1;
  typedef typename ublas::matrix<T>::const_reverse_iterator2 const_reverse_iterator2;  
  typedef typename ublas::matrix<T>::reverse_iterator1 reverse_iterator1;
  typedef typename ublas::matrix<T>::reverse_iterator2 reverse_iterator2; 

  public:
  
    explicit GraphMatrix();
    explicit GraphMatrix(const size_type r, const size_type c);
    explicit GraphMatrix(const size_type r, const size_type c, const std::vector<std::vector<T>> elements, bool isRow);
    explicit GraphMatrix(const size_type r, const size_type c, const std::vector<std::vector<T>>* elements, bool isRow);
    explicit GraphMatrix(const ublas::matrix<T>& inMatrix);
    ~GraphMatrix();

    //methods
    void Clear();
    void Resize(const size_type r, const size_type c);
    size_type nRows() const { return nRows_; };
    size_type nColumns() const { return nColumns_; };
    T Get(const size_type r, const size_type c) const { return matrix_(r,c); };
    ublas::matrix<T> Get() const { return matrix_; };
    std::vector<T> GetRow(const size_type r);
    std::vector<T> GetColumn(const size_type c); 
    void SetZero(const size_type r, const size_type c) { matrix_.erase_element(r,c); };
    void SetRowZero(size_type r);
    void SetColumnZero(size_type c);
    void Set(const size_type r, const size_type c, const T val) { matrix_.insert_element(r,c,val); };
    void SetRow(size_type r, std::vector<T> row);
    void SetRow(size_type r, std::vector<T>* row);
    void SetColumn(size_type c, std::vector<T> column);
    void SetColumn(size_type c, std::vector<T>* column);
    void SetRows(const size_type r, const size_type c, std::vector<std::vector<T>> rows);
    void SetRows(const size_type r, const size_type c, std::vector<std::vector<T>>* rows);
    void SetColumns(const size_type r, const size_type c, std::vector<std::vector<T>> columns);
    void SetColumns(const size_type r, const size_type c, std::vector<std::vector<T>>* columns);
    int nZeros(int i, bool isRow);
    int nZeros(std::vector<T> elements);
    int nZeros(std::vector<T>* elements);
    bool AllZeros(std::vector<T> elements);
    bool AllZeros(std::vector<T>* elements);
    GraphMatrix Unit(const size_type s) { return GraphMatrix(ublas::identity_matrix(s)); };  
    GraphMatrix Zero(const size_type r, const size_type c) { return GraphMatrix(ublas::zero_matrix(r,c)); };   
    GraphMatrix Transpose() { return GraphMatrix(trans(matrix_)); };
    GraphMatrix ReduceElements(T maxVal, T defaultGood, std::vector<T>& threshold, bool isRow);
    GraphMatrix RemoveDuplicates(T val, bool isRow);
    
    //iterators
    const_iterator1 begin1() const { return matrix_.begin1(); };
    const_iterator2 begin2() const { return matrix_.begin2(); };
    const_iterator1 end1() const { return matrix_.end1(); };
    const_iterator2 end2() const { return matrix_.end2(); }; 
    iterator1 begin1() { return matrix_.begin1(); };
    iterator2 begin2() { return matrix_.begin2(); };
    iterator1 end1() { return matrix_.end1(); };
    iterator2 end2() { return matrix_.end2(); };
    const_reverse_iterator1 rbegin1() const { return matrix_.rbegin1(); };
    const_reverse_iterator2 rbegin2() const { return matrix_.rbegin2(); };
    const_reverse_iterator1 rend1() const { return matrix_.rend1(); };
    const_reverse_iterator2 rend2() const { return matrix_.rend2(); }; 
    reverse_iterator1 rbegin1() { return matrix_.rbegin1(); };
    reverse_iterator2 rbegin2() { return matrix_.rbegin2(); };
    reverse_iterator1 rend1() { return matrix_.rend1(); };
    reverse_iterator2 rend2() { return matrix_.rend2(); };
    
    
    //operators
    GraphMatrix& operator=(const GraphMatrix<T>& other);
    bool operator==(const GraphMatrix<T>& other);
    bool operator!=(const GraphMatrix<T>& other);
    GraphMatrix operator+(const GraphMatrix<T>& other) { return GraphMatrix(matrix_+other.Get()); };
    GraphMatrix operator-(const GraphMatrix<T>& other) { return GraphMatrix(matrix_-other.Get()); };
    GraphMatrix operator*(const GraphMatrix<T>& other) { return GraphMatrix(prod(matrix_,other.Get())); };
    GraphMatrix operator*(const T& other) { return GraphMatrix(matrix_*=other); };
    GraphMatrix operator/(const T& other) { return GraphMatrix(matrix_/=other); };
    std::vector<T> operator*(const std::vector<T>& other) { return ToVector(prod(matrix_,ToVector(other))); };
    
  private:
  
    size_type nRows_;
    size_type nColumns_;
    ublas::matrix<T> matrix_;

    //methods
    ublas::vector<T> ToVector(const std::vector<T> vector); 
    std::vector<T> ToVector(const ublas::vector<T> vector); 
    std::vector<T> ToVector(const ublas::matrix_row<ublas::matrix<T>> row);
    std::vector<T> ToVector(const ublas::matrix_column<ublas::matrix<T>> column); 
};

template<typename T> GraphMatrix<T>::GraphMatrix()
{
   nRows_ = 0;
   nColumns_ = 0; 
   matrix_ = ublas::zero_matrix<T> ();
}

template<typename T> GraphMatrix<T>::GraphMatrix(const size_type r, const size_type c)
{
   nRows_ = r;
   nColumns_ = c; 
   matrix_ = ublas::zero_matrix<T> (nRows_,nColumns_);
}

template<typename T> GraphMatrix<T>::GraphMatrix(const size_type r, const size_type c, const std::vector<std::vector<T>> elements, bool isRow)
{
   matrix_.resize(r,c);
   nRows_ = r;
   nColumns_ = c; 

   if(isRow)
   {
      std::vector<std::vector<T>> rows_ = elements;
      rows_.resize(nRows_);
      for(typename std::vector<std::vector<T>>::iterator it = rows_.begin() ; it != rows_.end(); ++it)
          (*it).resize(nColumns_);

      typename std::vector<std::vector<T>>::iterator it;
      typename std::vector<T>::iterator jt;
      for(it = rows_.begin() ; it != rows_.end(); ++it)
          for(jt = (*it).begin(); jt != (*it).end(); jt++)
              matrix_(it-rows_.begin(),jt-(*it).begin()) = *jt;
   }else{
      std::vector<std::vector<T>> columns_ = elements;
      columns_.resize(nColumns_);
      for(typename std::vector<std::vector<T>>::iterator it = columns_.begin() ; it != columns_.end(); ++it)
          (*it).resize(nRows_);

      typename std::vector<std::vector<T>>::iterator it;
      typename std::vector<T>::iterator jt;
      for(it = columns_.begin() ; it != columns_.end(); ++it)
          for(jt = (*it).begin(); jt != (*it).end(); jt++)
              matrix_(jt-(*it).begin(),it-columns_.begin()) = *jt;
   } 
}

template<typename T> GraphMatrix<T>::GraphMatrix(const size_type r, const size_type c, const std::vector<std::vector<T>>* elements, bool isRow)
{
   matrix_.resize(r,c);
   nRows_ = r;
   nColumns_ = c; 

   if(isRow)
   {
      std::vector<std::vector<T>> rows_ = *elements;
      rows_.resize(nRows_);
      for(typename std::vector<std::vector<T>>::iterator it = rows_.begin() ; it != rows_.end(); ++it)
          (*it).resize(nColumns_);

      typename std::vector<std::vector<T>>::iterator it;
      typename std::vector<T>::iterator jt;
      for(it = rows_.begin() ; it != rows_.end(); ++it)
          for(jt = (*it).begin(); jt != (*it).end(); jt++)
              matrix_(it-rows_.begin(),jt-(*it).begin()) = *jt;
   }else{
      std::vector<std::vector<T>> columns_ = *elements;
      columns_.resize(nColumns_);
      for(typename std::vector<std::vector<T>>::iterator it = columns_.begin() ; it != columns_.end(); ++it)
          (*it).resize(nRows_);

      typename std::vector<std::vector<T>>::iterator it;
      typename std::vector<T>::iterator jt;
      for(it = columns_.begin() ; it != columns_.end(); ++it)
          for(jt = (*it).begin(); jt != (*it).end(); jt++)
              matrix_(jt-(*it).begin(),it-columns_.begin()) = *jt;
   } 
}

template<typename T> GraphMatrix<T>::GraphMatrix(const ublas::matrix<T>& inMatrix)
{
   nRows_ = inMatrix.size1();
   nColumns_ = inMatrix.size2();
   matrix_ = inMatrix;
}

template<typename T> GraphMatrix<T>::~GraphMatrix()
{
    nRows_ = 0; 
    nColumns_ = 0;
    matrix_.clear();
}

//operators
template<typename T> std::ostream& operator<<(std::ostream &s, const GraphMatrix<T>&m)
{
   s << "GraphMatrix: " << std::endl;
   for(auto is = m.begin1(); is!=m.end1(); is++){
       for(auto ic = is.begin(); ic!=is.end(); ic++){
            s << std::setprecision(3) << *ic << " ";
        }
        s << endl;
   }
   s << endl;
   return s;
}

template<typename T> GraphMatrix<T>& GraphMatrix<T>::operator=(const GraphMatrix<T>& other) 
{ 
   if(this == &other) return *this;
   this->Resize(other.nRows_,other.nColumns_);
   matrix_ = other.Get();
   return *this;
}

template<typename T> bool GraphMatrix<T>::operator==(const GraphMatrix<T>& other)  
{
   if(matrix_.size1() != other.nRows() || matrix_.size2() != other.nColumns()) return false;
   else{
      for(auto is = matrix_.begin1() ; is != matrix_.end1(); ++is)
          for(auto ic = is.begin(); ic != is.end(); ic++)
              if(matrix_(is-matrix_.begin1(),ic-is.begin())!=other.Get(is-matrix_.begin1(),ic-is.begin())) return false;
   }
   return true;
}

template<typename T> bool GraphMatrix<T>::operator!=(const GraphMatrix<T>& other) { 

   if(matrix_.size1() != other.nRows() || matrix_.size2() != other.nColumns()) return true;
   else{
      for(auto is = matrix_.begin1() ; is != matrix_.end1(); ++is)
          for(auto ic = is.begin(); ic != is.end(); ic++)
              if(matrix_(is-matrix_.begin1(),ic-is.begin())!=other.Get(is-matrix_.begin1(),ic-is.begin())) return true;
   }
   return false;
}


//other methods
template<typename T> void GraphMatrix<T>::Clear()
{
    nRows_ = 0; 
    nColumns_ = 0;
    matrix_.clear();
}

template<typename T> void GraphMatrix<T>::Resize(const size_type r, const size_type c)
{
   nRows_ = r;
   nColumns_ = c; 
   matrix_.resize(nRows_,nColumns_);
}

template<typename T> std::vector<T> GraphMatrix<T>::GetRow(const size_type r) 
{ 
   size_type row_ = r;
   if(row_>=nRows_) row_ = nRows_-1;  
   return ToVector(row(matrix_,row_)); 

}

template<typename T> std::vector<T> GraphMatrix<T>::GetColumn(const size_type c) 
{ 
   size_type column_ = c;
   if(column_>=nColumns_) column_ = nColumns_-1;  
   return ToVector(column(matrix_,column_)); 
}; 

template<typename T> void GraphMatrix<T>::SetRowZero(size_type r)
{
   std::vector<T> row_;
   row_.resize(nColumns_);
   if(r>=nRows_){ 
      Resize(nRows_+1,nColumns_);
      r = nRows_-1;
   }

   for(typename std::vector<T>::iterator it = row_.begin() ; it != row_.end(); ++it)
       matrix_.erase_element(r,it-row_.begin());
   
   nRows_ = matrix_.size1();
   nColumns_ = matrix_.size2();
}

template<typename T> void GraphMatrix<T>::SetRow(size_type r, std::vector<T> row)
{
   std::vector<T> row_ = row;
   row_.resize(nColumns_);
   if(r>=nRows_){ 
      Resize(nRows_+1,nColumns_);
      r = nRows_-1;
   }

   for(typename std::vector<T>::iterator it = row_.begin() ; it != row_.end(); ++it)
       matrix_.insert_element(r,it-row_.begin(),*it);

   nRows_ = matrix_.size1();
   nColumns_ = matrix_.size2();
}

template<typename T> void GraphMatrix<T>::SetRow(size_type r, std::vector<T>* row)
{
   std::vector<T> row_ = *row;
   row_.resize(nColumns_);
   if(r>=nRows_){ 
      Resize(nRows_+1,nColumns_);
      r = nRows_-1;
   }

   for(typename std::vector<T>::iterator it = row_.begin() ; it != row_.end(); ++it)
       matrix_.insert_element(r,it-row_.begin(),*it);

   nRows_ = matrix_.size1();
   nColumns_ = matrix_.size2();
}

template<typename T> void GraphMatrix<T>::SetColumnZero(size_type c)
{
   std::vector<T> column_;
   column_.resize(nRows_);
   if(c>=nColumns_){ 
      Resize(nRows_,nColumns_+1);
      c = nColumns_-1;
   }

   for(typename std::vector<T>::iterator it = column_.begin() ; it != column_.end(); ++it)
       matrix_.erase_element(it-column_.begin(),c);

   nRows_ = matrix_.size1();
   nColumns_ = matrix_.size2();
}

template<typename T> void GraphMatrix<T>::SetColumn(size_type c, std::vector<T> column)
{
   std::vector<T> column_ = column;
   column_.resize(nRows_);
   if(c>=nColumns_){ 
      Resize(nRows_,nColumns_+1);
      c = nColumns_-1;
   }

   for(typename std::vector<T>::iterator it = column_.begin() ; it != column_.end(); ++it)
       matrix_.insert_element(it-column_.begin(),c,*it);

   nRows_ = matrix_.size1();
   nColumns_ = matrix_.size2();
}

template<typename T> void GraphMatrix<T>::SetColumn(size_type c, std::vector<T>* column)
{
   std::vector<T> column_ = *column;
   column_.resize(nRows_);
   if(c>=nColumns_){ 
      Resize(nRows_,nColumns_+1);
      c = nColumns_-1;
   }

   for(typename std::vector<T>::iterator it = column_.begin() ; it != column_.end(); ++it)
       matrix_.insert_element(it-column_.begin(),c,*it);

   nRows_ = matrix_.size1();
   nColumns_ = matrix_.size2();
}

template<typename T> int GraphMatrix<T>::nZeros(int i, bool isRow)
{
   std::vector<T> vector_;
   if(isRow) vector_ = this->GetRow(i);
   else vector_ = this->GetColumn(i);

   int n = std::count(vector_.begin(), vector_.end(), T(0));
   return n;
}

template<typename T> int GraphMatrix<T>::nZeros(std::vector<T> elements)
{
   int n = std::count(elements.begin(), elements.end(), T(0));
   return n;
}

template<typename T> int GraphMatrix<T>::nZeros(std::vector<T>* elements)
{
   int n = std::count(elements->begin(), elements->end(), T(0));
   return n;
}

template<typename T> bool GraphMatrix<T>::AllZeros(std::vector<T> elements)
{
   bool isZeros = false;
   if(nZeros(&elements) == (int)elements.size()) isZeros = true;
   return isZeros;
}

template<typename T> bool GraphMatrix<T>::AllZeros(std::vector<T>* elements)
{
   bool isZeros = false;
   if(nZeros(elements) == (int)elements->size()) isZeros = true;
   return isZeros;
}

template<typename T> GraphMatrix<T> GraphMatrix<T>::ReduceElements(T maxVal, T defaultGood, std::vector<T>& thresholds, bool isRow)
{
   GraphMatrix<T> m_(nRows_,nColumns_);
   
   size_type n;
   if(isRow) n = nRows_;
   else n = nColumns_;
   thresholds.resize(n);
   
   for(size_type i=0; i<n; i++)
   {
       std::vector<T> vector_;
       if(isRow) vector_ = this->GetRow(i);
       else vector_ = this->GetColumn(i);
   
       T threshold = thresholds[i];
       if(AllZeros(&vector_)){
          if(isRow) m_.SetRow(i,vector_);
          else m_.SetColumn(i,vector_);
       }else{
          std::replace_if(vector_.begin(), vector_.end(), [threshold](T &v){ return v<=threshold; }, T(0));
          std::replace_if(vector_.begin(), vector_.end(), [threshold](T &v){ return v>threshold; }, defaultGood);
       
          if(isRow) m_.SetRow(i,vector_);
          else m_.SetColumn(i,vector_);
       } 
   }
   
   return m_;
}

template<typename T> GraphMatrix<T> GraphMatrix<T>::RemoveDuplicates(T val, bool isRow)
{
   GraphMatrix<T> m_(nRows_,nColumns_);
   
   size_type n;
   if(isRow) n = nRows_;
   else n = nColumns_;
   
   for(size_type i=0; i<n; i++)
   {
       std::vector<T> vector_;
       if(isRow) vector_ = this->GetRow(i);
       else vector_ = this->GetColumn(i);

       if(AllZeros(&vector_)){
          if(isRow) m_.SetRow(i,vector_);
          else m_.SetColumn(i,vector_);
       }else{ 
          typename std::vector<T>::iterator it = std::find(vector_.begin(), vector_.end(),val);

          T val = *it;
          int index = it-vector_.begin();
   
          vector_ = std::vector<T>(vector_.size(),T(0));
          if(it!=vector_.end()) vector_[index] = val;
       
          if(isRow) m_.SetRow(i,vector_);
          else m_.SetColumn(i,vector_);
       }
   }
   
   return m_;
}

//Private methods 
template<typename T> ublas::vector<T> GraphMatrix<T>::ToVector(const std::vector<T> vector)
{
   ublas::vector<T> v(vector.size());
   std::copy(vector.begin(),vector.end(),v.begin());
   return v;
} 

template<typename T> std::vector<T> GraphMatrix<T>::ToVector(const ublas::vector<T> vector)
{
   std::vector<T> v(vector.size());
   std::copy(vector.begin(),vector.end(),v.begin());
   return v;
} 


template<typename T> std::vector<T> GraphMatrix<T>::ToVector(const ublas::matrix_row<ublas::matrix<T>> row)
{
   std::vector<T> v(row.size());
   std::copy(row.begin(),row.end(),v.begin());
   return v;
} 

template<typename T> std::vector<T> GraphMatrix<T>::ToVector(const ublas::matrix_column<ublas::matrix<T>> column)
{
   std::vector<T> v(column.size());
   std::copy(column.begin(),column.end(),v.begin());
   return v;
} 


#endif
