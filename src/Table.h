#ifndef _TABLE_H_
#define _TABLE_H_

class Table {
  public:
  double *x, *y;
  Table(const int &dim){
    x = new double [dim];
    y = new double [dim];
  }
  
  virtual ~Table(){
    delete x;
    delete y;
  }
};
#endif
