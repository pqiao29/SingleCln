#include <Rcpp.h>
using namespace Rcpp;

#include <map>
#include <vector>
using namespace std; 

map<int, vector<int> > grid_ring(int n, int r = 1, bool center = true){
  map<int, vector<int> > ret;
  
  int row_n, col_n;
  for(int a = 1; a != n*n + 1; ++a){
    
    if(center) ret[a].push_back(a);
    
    // get row_n and col_n
    if ( a % n == 0) {
      row_n = n;
      col_n = a/n;
    } else {
      row_n = a % n;
      col_n = a/n + 1;
    }
    
    // get space_left
    vector<int> space_left;
    space_left.push_back(row_n - 1);
    space_left.push_back(n - row_n);
    space_left.push_back(col_n - 1);
    space_left.push_back(n - col_n);
    
    // Get neighbours:
    // furtherest horizontal or vertical neighbours
    if(space_left[0] >= r) ret[a].push_back(a - r);
    if(space_left[1] >= r) ret[a].push_back(a + r);
    if(space_left[2] >= r) ret[a].push_back(a - r*n);
    if(space_left[3] >= r) ret[a].push_back(a + r*n);
    
    // up and down edges
    if(r >= 2){
      for(int h = 1; h != r; ++h){
        if(space_left[2] >= h){
          if(space_left[0] >= r) ret[a].push_back(a - h*n - r);
          if(space_left[1] >= r) ret[a].push_back(a - h*n + r);
        }
        if(space_left[3] >= h){
          if(space_left[0] >= r) ret[a].push_back(a + h*n - r);
          if(space_left[1] >= r) ret[a].push_back(a + h*n + r);
        }
      }
    }
    
    // left and right edges
    for(int v = 1; v != r + 1; ++v){
      if(space_left[2] >= r){
        if(space_left[0] >= v) ret[a].push_back(a - r*n - v);
        if(space_left[1] >= v) ret[a].push_back(a - r*n + v);
      }
      if(space_left[3] >= r){
        if(space_left[0] >= v) ret[a].push_back(a + r*n - v);
        if(space_left[1] >= v) ret[a].push_back(a + r*n + v);
      }
    }
  }
  
  return ret;
}
